import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,time,resource,argparse,numpy as np,scipy.linalg as la
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
from ed import Model
from cyclic_model import CyclicModel
from observables import convert_basis,inversion_vector,total_spin_check
from sectors import point_permutations
ROOT=pathlib.Path(__file__).resolve().parent
class Limit(Exception):pass
def pack(x):return {'real':np.asarray(x).real.tolist(),'imag':np.asarray(x).imag.tolist()}
def run(r,count=2):
 j=.2375;sgn=(-1)**r
 from direct_inputs import inputs
 m,projected,q,E0=inputs(j,r);Q=[q]
 for idx in range(1,count+1):
  name=f'j{j:g}_r{r}_root{idx}';path=ROOT/(name+'.json')
  if path.exists() and (ROOT/(name+'.npz')).exists() and json.loads(path.read_text()).get('verified'):
   Q.append(np.load(ROOT/(name+'.npz'))['vector']);continue
  start=time.perf_counter();calls=0
  def mv(v):
   nonlocal calls
   if time.perf_counter()-start>1200:raise Limit()
   y=m.mv(v)
   for q in Q:y+=64*q*np.vdot(q,v)
   calls+=1
   if calls%25==0:print(json.dumps({'name':name,'H_calls':calls,'seconds':time.perf_counter()-start}),flush=True)
   return y
  A=LinearOperator((m.nb,m.nb),matvec=mv,dtype=complex);v=projected@np.array([.13+.17j,.29-.11j,1.])
  for q in Q:v-=q*np.vdot(q,v)
  v/=la.norm(v);result={'name':name,'Jp':j,'reflection_character':r,'reflection_parity':sgn,'cell_n':[1,0],'dimension':m.nb,'ncv':64,'verified':False,'previous_physical_states_deflated':len(Q),'shift':64}
  try:
   e,U=eigsh(A,k=1,which='SA',ncv=64,tol=2e-12,maxiter=250,v0=v);u=U[:,0];energy=float(np.vdot(u,m.mv(u)).real);res=float(la.norm(m.mv(u)-energy*u));spin=total_spin_check(m,u);a=u.conj()@projected;orth=max(abs(np.vdot(q,u)) for q in Q)
   result.update(energy=energy,omega=energy-E0,H_residual=res,spin=spin,orthogonality_to_deflated=orth,amplitudes=pack(a),residue_matrix=pack(np.outer(a.conj(),a)),verified=bool(res<1e-8 and spin['S2_eigen_residual']<1e-5 and orth<1e-7))
   if result['verified']:Q.append(u);np.savez_compressed(ROOT/(name+'.npz'),vector=u,energy=energy,q=m.q)
  except (Limit,ArpackNoConvergence) as ex:result['status']=type(ex).__name__+'_NO_VERIFIED_PAIR'
  result.update(H_calls=calls,wall_seconds=time.perf_counter()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024);path.write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
  if not result['verified']:break
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('r',type=int);p.add_argument('--count',type=int,default=2);a=p.parse_args();run(a.r,a.count)
