import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,numpy as np,json,time
from cyclic_model import CyclicModel
from observables import total_spin_check
ROOT=pathlib.Path(__file__).resolve().parent
from cross_scan import OLD
plan=[((0,0),r) for r in [0,1,2,3]]+[((0,1),r) for r in [0,1]]+[((1,2),r) for r in [0,1,2]]
out=[]
for n,r in plan:
 m=CyclicModel(jp=.2,cell_n=n,angular=r)
 for jp in [.2,.5]:
  name=f'j{jp:g}_n{n[0]}{n[1]}_r{r}';z=np.load(OLD/(name+'.npz'));v=z['vector'];spin=total_spin_check(m,v);record={'Jp':jp,'cell_n':n,'character':r,'energy':float(z['energy']),**spin};out.append(record);print(json.dumps(record),flush=True)
 (ROOT/'endpoint_spin_profiles.json').write_text(json.dumps(out,indent=2))
 del m
