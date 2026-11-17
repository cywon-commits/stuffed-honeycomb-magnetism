"""Reuse verified cell-star eigenvectors via exact point-group maps, then directly recheck."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,time,json,resource,numpy as np,scipy.linalg as la
from ed import Model
from cyclic_model import CyclicModel
from sectors import point_permutations
from observables import full_samples,total_spin_check
from crossing_transitions import pack
ROOT=pathlib.Path(__file__).resolve().parent;DATA=ROOT.parent/'crossing'
def main():
 start=time.perf_counter();result=json.loads((ROOT/'weak_moments.json').read_text());z=np.load(ROOT/'weak_seeds.npz');B=z['B'];E=float(z['E0']);m=Model(primitive=False,jp=.2375);m.sector(z['q']);poles=[]
 for r in [0,1]:
  name=f'p0.237500000_n01_r{r}_e0';old=np.load(DATA/(name+'.npz'));s=CyclicModel(jp=.2375,cell_n=(0,1),angular=r);oldchars=np.exp(-2j*np.pi*(m.gs@s.q));newchars=np.exp(-2j*np.pi*(m.gs@m.q));lookup={tuple(p):i for i,p in enumerate(m.perms)}
  for index,p in enumerate(point_permutations()):
   inv=np.argsort(p)
   if all(abs(oldchars[lookup[tuple(inv[t[p]])]]-newchars[i])<1e-8 for i,t in enumerate(m.perms)):break
  else:raise RuntimeError('No star map')
  lo=np.zeros(1<<m.split,np.uint32);hi=np.zeros(1<<(m.n-m.split),np.uint32)
  for table,offset in [(lo,0),(hi,m.split)]:
   for mask in range(1,len(table)):
    bit=mask&-mask;table[mask]=table[mask^bit]|(1<<int(inv[offset+bit.bit_length()-1]))
  pre=lo[m.masks&((1<<m.split)-1)]|hi[m.masks>>m.split];v=full_samples(s,old['vector'],pre)*np.sqrt(m.lengths.astype(float));norm=la.norm(v);assert abs(norm-1)<1e-10;v/=norm;energy=float(old['energy']);res=float(la.norm(m.mv(v)-energy*v));spin=total_spin_check(m,v);a=v.conj()@B;verified=res<1e-8 and spin['S2_eigen_residual']<1e-5;assert verified
  poles.append({'energy':energy,'omega':energy-E,'H_residual':res,'spin':spin,'verified':verified,'amplitudes':pack(a),'residue_matrix':pack(np.outer(a.conj(),a)),'source_checkpoint':name,'point_group_map_index':index,'status':'SYMMETRY_MAPPED_AND_FRESHLY_VERIFIED'})
  np.savez_compressed(ROOT/f'weak_pole_{r}.npz',vector=v,energy=energy,q=m.q)
 result.update(poles=sorted(poles,key=lambda x:x['energy']),search_start='Exact spatial-star map of previously verified minima in both compatible reflection characters; each vector rechecked in target cell sector.',ncv=None,arpack_returned_all=None,H_calls=len(poles),S2_calls=len(poles),operator_calls=m.calls,solve_seconds=time.perf_counter()-start,total_wall_seconds=time.perf_counter()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,scope='Two verified symmetry-character minima; partial spectral matrix. Higher eigenpairs not enumerated.')
 sums=[np.zeros((3,3),complex) for _ in range(3)]
 for d in poles:
  R=np.array(d['residue_matrix']['real'])+1j*np.array(d['residue_matrix']['imag'])
  for power in range(3):sums[power]+=d['omega']**power*R
 result['partial_spectral_moment_sums']={f'M{p}':pack(sums[p]) for p in range(3)};result['moment_remainders']={f'M{p}':pack(np.array(result['moments']['1'][f'M{p}']['real'])+1j*np.array(result['moments']['1'][f'M{p}']['imag'])-sums[p]) for p in range(3)}
 (ROOT/'symmetry_weak_spectrum.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=='__main__':main()
