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
 start=time.perf_counter();j=.2375;result=json.loads((ROOT/'weak_moments.json').read_text());moments=result['moments'];E=result['E0'];z=np.load(ROOT/'weak_seeds.npz');B=z['B'];target=Model(primitive=False,jp=j);target.sector(z['q']);result['search_start']='generic complex linear combination of n,m,h seeds';result['ncv']=64
 calls=0;t0=time.perf_counter()
 def mv(x):
  nonlocal calls
  if time.perf_counter()-t0>2400:raise Limit()
  calls+=1
  if calls%25==0:print(json.dumps({'H_calls':calls,'seconds':time.perf_counter()-t0}),flush=True)
  return target.mv(x)
 A=LinearOperator((target.nb,target.nb),matvec=mv,dtype=complex);rng=np.random.default_rng(2026);roots=[];converged=False
 try:evals,U=eigsh(A,k=6,which='SA',ncv=64,tol=1e-11,maxiter=1000,v0=B@np.array([.17+.13j,.23-.11j,1.]));converged=True
 except ArpackNoConvergence as e:evals,U=e.eigenvalues,e.eigenvectors
 except Limit:evals=[];U=None
 for index,e in enumerate(evals):
  u=U[:,index];res=float(la.norm(target.mv(u)-e*u));spin=total_spin_check(target,u);a=u.conj()@B;R=np.outer(a.conj(),a);verified=res<1e-8 and spin['S2_eigen_residual']<1e-5
  roots.append({'energy':float(e),'omega':float(e-E),'H_residual':res,'spin':spin,'verified':verified,'amplitudes':pack(a),'residue_matrix':pack(R)})
  if verified:np.savez_compressed(ROOT/f'weak_pole_{index}.npz',vector=u,energy=e,q=target.q)
 result.update(poles=sorted(roots,key=lambda d:d['energy']),arpack_returned_all=converged,H_calls=calls,solve_seconds=time.perf_counter()-t0,total_wall_seconds=time.perf_counter()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,scope='Six low eigenpairs searched with a probe-generated starting vector (only S=1/2,3/2 directly accessible); partial spectral matrix. Exact operator moments include all higher states. No broadened lines or damping inferred.')
 sums=[np.zeros((3,3),complex) for _ in range(3)]
 for d in roots:
  if d['verified']:
   r=np.array(d['residue_matrix']['real'])+1j*np.array(d['residue_matrix']['imag'])
   for p in range(3):sums[p]+=d['omega']**p*r
 result['partial_spectral_moment_sums']={f'M{p}':pack(sums[p]) for p in range(3)}
 result['moment_remainders']={f'M{p}':pack(np.array(moments['1'][f'M{p}']['real'])+1j*np.array(moments['1'][f'M{p}']['imag'])-sums[p]) for p in range(3)}
 (ROOT/'weak_spectrum.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=='__main__':main()
