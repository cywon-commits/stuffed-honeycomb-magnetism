"""Verify all available saved eigenvectors using a freshly constructed operator."""
from ed import Model
import pathlib,numpy as np,scipy.linalg as la,json
root=pathlib.Path(__file__).resolve().parent;m=Model();out=[]
for name in ['Gamma_primitive','K_primitive','cell11_fold1','cell11_fold2','cell11_fold3','other']:
 path=root/(name+'.npz')
 if not path.exists():continue
 z=np.load(path);r=json.loads((root/(name+'.json')).read_text());m.sector(z['q']);v=z['vector'];assert abs(la.norm(v)-1)<1e-12 and np.max(abs(z['q']-r['q']))<1e-14
 hv=m.mv(v);e=float(np.vdot(v,hv).real);res=float(la.norm(hv-e*v));assert res<1e-8 and abs(e-r['energy'])<1e-12
 out.append({'name':name,'energy':e,'residual':res,'status':'VERIFIED FROM SAVED VECTOR'})
assert out
(root/'checkpoint_verification.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
