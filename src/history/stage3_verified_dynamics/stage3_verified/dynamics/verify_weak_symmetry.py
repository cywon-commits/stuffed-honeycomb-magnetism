import pathlib,json,numpy as np,scipy.linalg as la
from ed import Model
from observables import inversion_vector
from sectors import point_permutations
ROOT=pathlib.Path(__file__).resolve().parent
m=Model(primitive=False,jp=.2375);z=np.load(ROOT/'weak_pole_0.npz');m.sector(z['q']);chars=np.exp(-2j*np.pi*(m.gs@m.q));lookup={tuple(p):i for i,p in enumerate(m.perms)}
for index,F in enumerate(point_permutations()[6:]):
 inv=np.argsort(F)
 if all(abs(chars[lookup[tuple(F[t[inv]])]]-chars[j])<1e-8 for j,t in enumerate(m.perms)):break
else:raise RuntimeError('No compatible reflection')
out=[]
for file in sorted(ROOT.glob('weak_pole_*.npz')):
 z=np.load(file);u=z['vector'];fu=inversion_vector(m,u,F);ov=np.vdot(u,fu);sgn=1 if ov.real>=0 else -1;res=float(la.norm(fu-sgn*u));out.append({'file':file.name,'reflection_index_in_pg6':index,'reflection_parity':sgn if res<1e-6 else None,'reflection_eigen_residual':res,'energy':float(z['energy']),'cell_n':[1,0]})
(ROOT/'weak_symmetry.json').write_text(json.dumps(out,indent=2));print(out)
