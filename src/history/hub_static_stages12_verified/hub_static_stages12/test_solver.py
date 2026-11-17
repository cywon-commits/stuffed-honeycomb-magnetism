import pathlib,json,numpy as np,scipy.linalg as la
from solve_sector import solve
from cyclic_model import CyclicModel
from ed import Model
out=[]
for cyclic,ns,ell in [(False,(0,0),0),(False,(1,0),0),(True,(0,0),0),(True,(0,0),3)]:
 r,v=solve(.2,*ns,L=2,up=6,ncv=16,seconds=30,cyclic=cyclic,angular=ell)
 m=CyclicModel(2,6,.2,ns,ell) if cyclic else Model(2,6,False,.2)
 k=np.array(ns)/2;m.sector([(k[0]+k[1])/3,(2*k[0]-k[1])/3]);H=np.column_stack([m.mv(np.eye(m.nb,dtype=complex)[:,i]) for i in range(m.nb)])
 error=abs(r['energy']-la.eigvalsh(H)[0]);assert r['converged'] and r['residual']<1e-8 and error<1e-10
 out.append({'cyclic':cyclic,'cell_n':ns,'angular':ell,'energy_error':float(error),'residual':r['residual'],'real_solver':r['real_solver']})
r,v=solve(.2,0,0,L=2,up=6,seconds=-1,cyclic=True);assert not r['converged'] and 'energy' not in r and v is None
out.append({'time_limit_never_labeled_converged':True})
pathlib.Path(__file__).with_name('solver_validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));print('RESTARTED SOLVER TESTS PASSED')
