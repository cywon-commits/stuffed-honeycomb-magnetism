"""Rebuild matrix-free operator and verify a saved Ritz vector independently of Lanczos history."""
import argparse,numpy as np,json,scipy.linalg as la
from ed import Model
p=argparse.ArgumentParser();p.add_argument('checkpoint');p.add_argument('--tol',type=float,default=1e-8);a=p.parse_args()
z=np.load(a.checkpoint);m=Model();m.sector(z['q']);v=z['vector'];assert len(v)==m.nb
assert abs(la.norm(v)-1)<1e-12
hv=m.mv(v);e=np.vdot(v,hv).real;res=la.norm(hv-e*v)
print(json.dumps({'energy':float(e),'residual':float(res),'stored_energy':float(z['energy']),'converged':bool(res<a.tol)}))
if res>=a.tol:raise SystemExit(1)
