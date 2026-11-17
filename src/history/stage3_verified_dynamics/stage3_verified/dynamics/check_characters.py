import pathlib,json,numpy as np,scipy.linalg as la
from cyclic_model import CyclicModel
from observables import inversion_vector
from sectors import point_permutations
ROOT=pathlib.Path(__file__).resolve().parent;DATA=ROOT.parent/'crossing';out=[]
for j in [.3875,.3921875]:
 for r in [0,3]:
  name=f'p{j:.9f}_n00_r{r}_e0';v=np.load(DATA/(name+'.npz'))['vector'];m=CyclicModel(jp=j,cell_n=(0,0),angular=r);P=point_permutations();sign=1 if r==0 else -1;rr=float(la.norm(inversion_vector(m,v,P[1])-sign*v));fr=float(la.norm(inversion_vector(m,v,P[6])-sign*v));assert max(rr,fr)<1e-8;out.append({'name':name,'C6_character':r,'reflection_parity':sign,'rotation_eigen_residual':rr,'reflection_eigen_residual':fr})
(ROOT/'point_group_validation.json').write_text(json.dumps(out,indent=2));print(out)
