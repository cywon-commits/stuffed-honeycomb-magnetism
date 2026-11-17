import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,scipy.linalg as la,json
from experiment import ROOT,Model,gram
from little_model import LittleModel
from observables import convert_basis
z=np.load(ROOT/'j1_G_seeds.npz');full=Model();full.sector(z['q']);B=z['B'];records=[];total=np.zeros((6,6),complex)
for a,p in [(0,1),(0,-1),(1,None),(2,None)]:
 m=LittleModel(z['q'],a,p);C=m.project(full,B);M=gram(C);total+=M;index=int(np.argmax(np.diag(M).real));v=C[:,index];v/=la.norm(v);back=convert_basis(m,full,v);hv=convert_basis(m,full,m.mv(v));err=la.norm(full.mv(back)-hv);norm=la.norm(back);assert err<1e-10 and abs(norm-1)<1e-10
 r={'angular':a,'reflection':p,'dimension':m.nb,'projected_M0':{'real':M.real.tolist(),'imag':M.imag.tolist()},'channel_weights':np.diag(M).real.tolist(),'H_intertwining_error':float(err),'norm':float(norm)};records.append(r);print(json.dumps(r),flush=True)
err=la.norm(total-z['M0']);assert err<1e-8;(ROOT/'little_group_validation.json').write_text(json.dumps({'status':'PASS','completeness_error':float(err),'blocks':records},indent=2))
