import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,time,resource,numpy as np,scipy.linalg as la
from ed import Model
from cyclic_model import CyclicModel
from observables import convert_basis,spin_seed,cell_q,total_spin_check,full_samples,coefficients
from sectors import point_permutations
from search import pack
ROOT=pathlib.Path(__file__).resolve().parent;OLD=ROOT.parent/'static'
def main():
 start=time.perf_counter();j=.2;z=np.load(OLD/'j0.2_n12_r0.npz');E=float(z['energy']);s=CyclicModel(jp=j,cell_n=(1,2),angular=0);source=Model(primitive=False,jp=j);source.sector(z['q']);v=convert_basis(s,source,z['vector']);del s;sourcecheck={'H_residual':float(la.norm(source.mv(v)-E*v)),'spin':total_spin_check(source,v)};assert sourcecheck['H_residual']<1e-8
 moments={}
 for sign in [1,-1]:
  k=np.array([0,sign/3]);m=Model(primitive=False,jp=j);m.sector(source.q+cell_q(k));B=np.column_stack([spin_seed(source,m,v,k,ch) for ch in ['n','m','H']]);HB=np.column_stack([m.mv(B[:,a])-E*B[:,a] for a in range(3)]);moments[str(sign)]={f'M{p}':pack(x) for p,x in enumerate([B.conj().T@B,B.conj().T@HB,HB.conj().T@HB])}
  if sign==1:target=m;forwardB=B
  else:del m,B,HB
 del source,v
 poles=[]
 for r in [0,1]:
  name=f'j0.2_n01_r{r}';z=np.load(OLD/(name+'.npz'));s=CyclicModel(jp=j,cell_n=(0,1),angular=r);oldchars=np.exp(-2j*np.pi*(target.gs@s.q));newchars=np.exp(-2j*np.pi*(target.gs@target.q));lookup={tuple(t):i for i,t in enumerate(target.perms)}
  for pg,p in enumerate(point_permutations()):
   inv=np.argsort(p)
   if all(abs(oldchars[lookup[tuple(inv[t[p]])]]-newchars[i])<1e-8 for i,t in enumerate(target.perms)):break
  else:raise RuntimeError('No map')
  lo=np.zeros(1<<target.split,np.uint32);hi=np.zeros(1<<(target.n-target.split),np.uint32)
  for table,offset in [(lo,0),(hi,target.split)]:
   for mask in range(1,len(table)):
    bit=mask&-mask;table[mask]=table[mask^bit]|(1<<int(inv[offset+bit.bit_length()-1]))
  pre=lo[target.masks&((1<<target.split)-1)]|hi[target.masks>>target.split];u=full_samples(s,z['vector'],pre)*np.sqrt(target.lengths.astype(float));assert abs(la.norm(u)-1)<1e-9;en=float(z['energy']);res=float(la.norm(target.mv(u)-en*u));spin=total_spin_check(target,u);a=u.conj()@forwardB;assert res<1e-8 and spin['S2_eigen_residual']<1e-5
  poles.append({'name':name,'Jp':j,'reflection_character':r,'reflection_parity':(-1)**r,'cell_n':[1,0],'energy':en,'omega':en-E,'H_residual':res,'spin':spin,'verified':True,'amplitudes':pack(a),'residue_matrix':pack(np.outer(a.conj(),a)),'source_point_group_map':pg});np.savez_compressed(ROOT/(name+'_target.npz'),vector=u,energy=en,q=target.q)
 result={'Jp':j,'E0':E,'source_cell_n':[1,2],'target_cell_n':[1,0],'probe_cell_n':[0,1],'channels':['n','m','h'],'source_check':sourcecheck,'moments':moments,'poles':poles,'wall_seconds':time.perf_counter()-start,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
 refs=[.33177055287103,.0069218759870593,.0075353314103994];result['hub_benchmark_errors']=[moments['1'][f'M{p}']['real'][2][2]-refs[p] for p in range(3)];assert max(abs(x) for x in result['hub_benchmark_errors'])<1e-10
 (ROOT/'matched_j0.2.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=='__main__':main()
