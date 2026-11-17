import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,numpy as np,scipy.linalg as la
from scipy.sparse.linalg import LinearOperator,eigsh
from cyclic_model import CyclicModel
ROOT=pathlib.Path(__file__).resolve().parent
m=CyclicModel(up=2,jp=.2375,cell_n=(1,0),angular=0);H=np.column_stack([m.mv(np.eye(m.nb,dtype=complex)[:,i]) for i in range(m.nb)]);exact=la.eigvalsh(H);Q=[];rows=[]
for k in range(3):
 def mv(v):
  y=m.mv(v)
  for q in Q:y+=64*q*np.vdot(q,v)
  return y
 A=LinearOperator(H.shape,matvec=mv,dtype=complex);rng=np.random.default_rng(2026+k);v=rng.normal(size=m.nb)+1j*rng.normal(size=m.nb)
 for q in Q:v-=q*np.vdot(q,v)
 e,U=eigsh(A,k=1,which='SA',ncv=min(20,m.nb),tol=1e-13,v0=v);u=U[:,0];energy=np.vdot(u,H@u).real;res=la.norm(H@u-energy*u);orth=max([abs(np.vdot(q,u)) for q in Q],default=0.);assert abs(energy-exact[k])<1e-10 and res<1e-10 and orth<1e-10;Q.append(u);rows.append({'root':k,'energy_error':float(abs(energy-exact[k])),'residual':float(res),'orthogonality':float(orth)})
for t in [.01,.04,.05,.1]:
 x=np.linspace(0,20.2,100001);p=x*(20.2+t-x)/(20.2*t);assert min(p)>-1e-10 and min(p[x>=t])>1-1e-10
(ROOT/'small_validation.json').write_text(json.dumps({'status':'PASS','complex_physical_deflation_tests':rows,'quadratic_tail_majorant_on_bounded_spectrum':True},indent=2));print('PASS',rows)
