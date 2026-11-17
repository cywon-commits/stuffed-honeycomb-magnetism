import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,time,resource,numpy as np,scipy.linalg as la
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
from ed import Model
from cyclic_model import CyclicModel
from observables import convert_basis,spin_seed,cell_q,total_spin_check
from crossing_transitions import pack
ROOT=pathlib.Path(__file__).resolve().parent;DATA=ROOT.parent/'crossing'
class Limit(Exception):pass
def main():
 start=time.perf_counter();j=.2375;z=np.load(DATA/'p0.237500000_n12_r0_e0.npz');E=float(z['energy']);source=Model(primitive=False,jp=j);source.sector(z['q']);c=CyclicModel(jp=j,cell_n=(1,2),angular=0);v=convert_basis(c,source,z['vector']);del c
 channels=['n','m','H'];moments={};target=None;B=None
 for sign in [1,-1]:
  k=np.array([0,sign/3]);t=Model(primitive=False,jp=j);t.sector(source.q+cell_q(k));b=np.column_stack([spin_seed(source,t,v,k,ch) for ch in channels]);hb=np.column_stack([t.mv(b[:,i])-E*b[:,i] for i in range(3)]);moments[str(sign)]={'M0':pack(b.conj().T@b),'M1':pack(b.conj().T@hb),'M2':pack(hb.conj().T@hb)}
  if sign==1:target=t;B=b;HB=hb
  else:del t,b,hb
 # Spin 1/2 makes diagonal longitudinal static moments equal to one third scalar moments.
 hf=moments['1']['M1']['real'][2][2];hr=moments['-1']['M1']['real'][2][2];m0=moments['1']['M0']['real'][2][2]
 result={'Jp':j,'source_cell_n':[1,2],'probe_cell_n':[0,1],'target_cell_n':[1,0],'source_M':.5,'E0':E,'dimension':target.nb,'channels':['n','m','h'],'moments':moments,'reference_checks':{'hub_M0_reference':.33103929,'hub_M0_error':m0-.33103929,'hub_sym_M1_reference':.01006870,'hub_sym_M1_error':(hf+hr)/2-.01006870,'hub_forward_M1':hf,'hub_reverse_M1':hr},'elastic_ground_manifold_weight':0.,'elastic_reason':'Target momentum differs from both ground corner partners.','moment_seconds':time.perf_counter()-start}
 (ROOT/'weak_moments.json').write_text(json.dumps(result,indent=2));np.savez_compressed(ROOT/'weak_seeds.npz',B=B,HB=HB,E0=E,q=target.q)
 if '--moments-only' in __import__('sys').argv:return
 del source,v;calls=0;t0=time.perf_counter()
 def mv(x):
  nonlocal calls
  if time.perf_counter()-t0>2400:raise Limit()
  calls+=1
  if calls%25==0:print(json.dumps({'H_calls':calls,'seconds':time.perf_counter()-t0}),flush=True)
  return target.mv(x)
 A=LinearOperator((target.nb,target.nb),matvec=mv,dtype=complex);rng=np.random.default_rng(2026);roots=[];converged=False
 try:evals,U=eigsh(A,k=6,which='SA',ncv=40,tol=1e-11,maxiter=1000,v0=rng.normal(size=target.nb)+1j*rng.normal(size=target.nb));converged=True
 except ArpackNoConvergence as e:evals,U=e.eigenvalues,e.eigenvectors
 except Limit:evals=[];U=None
 for index,e in enumerate(evals):
  u=U[:,index];res=float(la.norm(target.mv(u)-e*u));spin=total_spin_check(target,u);a=u.conj()@B;R=np.outer(a.conj(),a);verified=res<1e-8 and spin['S2_eigen_residual']<1e-5
  roots.append({'energy':float(e),'omega':float(e-E),'H_residual':res,'spin':spin,'verified':verified,'amplitudes':pack(a),'residue_matrix':pack(R)})
  if verified:np.savez_compressed(ROOT/f'weak_pole_{index}.npz',vector=u,energy=e,q=target.q)
 result.update(poles=sorted(roots,key=lambda d:d['energy']),arpack_returned_all=converged,H_calls=calls,solve_seconds=time.perf_counter()-t0,total_wall_seconds=time.perf_counter()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,scope='Lowest six target eigenpairs; partial spectral matrix. Exact operator moments include all higher states. No broadened lines or damping inferred.')
 sums=[np.zeros((3,3),complex) for _ in range(3)]
 for d in roots:
  if d['verified']:
   r=np.array(d['residue_matrix']['real'])+1j*np.array(d['residue_matrix']['imag'])
   for p in range(3):sums[p]+=d['omega']**p*r
 result['partial_spectral_moment_sums']={f'M{p}':pack(sums[p]) for p in range(3)}
 result['moment_remainders']={f'M{p}':pack(np.array(moments['1'][f'M{p}']['real'])+1j*np.array(moments['1'][f'M{p}']['imag'])-sums[p]) for p in range(3)}
 (ROOT/'weak_spectrum.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=='__main__':main()
