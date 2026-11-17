"""Memory-bounded restarted real-Hermitian Lanczos, with explicit final residual."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np, scipy.linalg as la
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
import argparse,pathlib,json,time,resource
from ed import Model
from cyclic_model import CyclicModel
ROOT=pathlib.Path(__file__).resolve().parent
class TimeLimit(Exception):pass

def solve(jp,n1,n2,ncv=32,seconds=900,L=3,up=14,seed=2026,cyclic=False,angular=0):
 t=time.perf_counter();m=CyclicModel(L,up,jp,(n1,n2),angular) if cyclic else Model(L,up,primitive=False,jp=jp);real_phase=np.max(abs(m.phase.imag))<1e-13 if cyclic else (n1==0 and n2==0);k=np.array([n1,n2])/L;q=[(k[0]+k[1])/3,(2*k[0]-k[1])/3];m.sector(q)
 setup=time.perf_counter()-t;history=[];calls=0;start=time.perf_counter()
 def mv(x):
  nonlocal calls
  if time.perf_counter()-start>seconds:raise TimeLimit()
  z=x.astype(complex) if real_phase else x[:m.nb]+1j*x[m.nb:];y=m.mv(z);calls+=1
  if calls%50==0:
   record={'matvec_calls':calls,'elapsed_seconds':time.perf_counter()-start,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024};history.append(record);print(json.dumps(record),flush=True)
  return y.real.copy() if real_phase else np.concatenate([y.real,y.imag])
 dim=m.nb if real_phase else 2*m.nb
 op=LinearOperator((dim,dim),matvec=mv,dtype=float)
 rng=np.random.default_rng(seed);v0=rng.standard_normal(dim);v0/=la.norm(v0)
 result={'Jp':jp,'cyclic':cyclic,'angular':angular,'cyclic_order':m.order if cyclic else 1,'real_solver':bool(real_phase),'L':L,'up':up,'cell_n':[n1,n2],'cell_k':k.tolist(),'q_argument':q,'dimension':m.nb,'setup_seconds':setup,'ncv':ncv,'arpack_tol':1e-11,'residual_tolerance':1e-8,'converged':False,'seed':seed,'iteration_count':None,'iteration_count_note':'SciPy does not expose restart iterations; matvec calls are measured.'}
 v=None
 try:
  e,vs=eigsh(op,k=1,which='SA',ncv=ncv,tol=1e-11,maxiter=2000,v0=v0)
  v=vs[:,0].astype(complex) if real_phase else vs[:m.nb,0]+1j*vs[m.nb:,0]
 except ArpackNoConvergence as exc:
  result['status']='ARPACK_NO_CONVERGENCE'
  if len(exc.eigenvalues):v=exc.eigenvectors[:,0].astype(complex) if real_phase else exc.eigenvectors[:m.nb,0]+1j*exc.eigenvectors[m.nb:,0]
 except TimeLimit:result['status']='WALL_TIME_LIMIT_NO_VERIFIED_EIGENPAIR'
 if v is not None:
  v/=la.norm(v);hv=m.mv(v);energy=np.vdot(v,hv).real;res=la.norm(hv-energy*v);result.update(energy=float(energy),residual=float(res),converged=bool(res<1e-8));result['status']='CONVERGED' if result['converged'] else 'UNCONVERGED'
 result.update(matvec_calls=calls,direct_residual_matvecs=1 if v is not None else 0,matvec_seconds=m.mv_seconds,solve_seconds=time.perf_counter()-start,history=history,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
 return result,v

def main():
 p=argparse.ArgumentParser();p.add_argument('--cyclic',action='store_true');p.add_argument('--angular',type=int,default=0);p.add_argument('--jp',type=float,required=True);p.add_argument('--n',type=int,nargs=2,default=[0,0]);p.add_argument('--seconds',type=float,default=900);p.add_argument('--ncv',type=int,default=64);a=p.parse_args()
 name=f'j{a.jp:g}_n{a.n[0]}{a.n[1]}'+(f'_r{a.angular}' if a.cyclic else '');r,v=solve(a.jp,*a.n,ncv=a.ncv,seconds=a.seconds,cyclic=a.cyclic,angular=a.angular)
 if v is not None:np.savez_compressed(ROOT/(name+'.npz'),vector=v,energy=r['energy'],q=r['q_argument'],cyclic=a.cyclic,angular=a.angular,cell_n=a.n,jp=a.jp)
 (ROOT/(name+'.json')).write_text(json.dumps(r,indent=2));print(json.dumps(r),flush=True)
if __name__=='__main__':main()
