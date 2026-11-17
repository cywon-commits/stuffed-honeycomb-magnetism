import numpy as np, scipy.sparse as sp, scipy.linalg as la, json
from itertools import combinations
from ed import Model,geometry,translations,lanczos

def independent_bonds(L):
 c,_=geometry(L);N=len(c);out=[]
 # Euclidean triangular nearest-neighbor criterion in oblique 60-degree basis.
 for i in range(N):
  for j in range(i+1,N):
   ds=[c[j]-c[i]+np.array([L*x+2*L*y,L*x-L*y]) for x in range(-2,3) for y in range(-2,3)]
   if any(u*u+v*v+u*v==1 for u,v in ds):out.append((i,j,1.))
 return out

def full_h(n,up,bonds):
 masks=sorted(sum(1<<p for p in ps) for ps in combinations(range(n),up));idx={m:i for i,m in enumerate(masks)};rows=[];cols=[];vals=[]
 for i,m in enumerate(masks):
  d=0.
  for a,b,j in bonds:
   anti=((m>>a)^(m>>b))&1;d+=j*(.25-.5*anti)
   if anti:rows.append(idx[m^(1<<a)^(1<<b)]);cols.append(i);vals.append(j/2)
  rows.append(i);cols.append(i);vals.append(d)
 return masks,idx,sp.coo_matrix((vals,(rows,cols)),shape=(len(masks),len(masks))).tocsr()

def run():
 results=[]
 # Extract only the uploaded geometry definition, without executing its imports or solver.
 import ast,pathlib
 tree=ast.parse(pathlib.Path(__file__).with_name('original_tri27_matrix_free_k_ed.py').read_text())
 start=next(i for i,n in enumerate(tree.body) if isinstance(n,ast.FunctionDef) and n.name=='sid')
 stop=next(i for i,n in enumerate(tree.body) if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='trans_lo' for t in n.targets))
 body=[n for n in tree.body[start:stop] if not isinstance(n,(ast.Expr,ast.Assert))]
 ns={'np':np,'J':1.,'Jp':1.}
 exec(compile(ast.Module(body=body,type_ignores=[]),'<uploaded geometry only>','exec'),ns)
 assert sorted(map(tuple,ns['bonds']))==sorted(map(tuple,geometry(3)[1]))
 results.append({'uploaded_geometry_exact_match':True})
 # Three-site quotient retains periodic edge multiplicities: each triangle edge has weight 3J.
 tiny=Model(1,2)
 _,_,Htiny=full_h(3,2,[(0,1,3.),(0,2,3.),(1,2,3.)])
 assert np.max(np.abs(la.eigvalsh(Htiny.toarray())-[-2.25,-2.25,2.25]))<1e-12
 es=[]
 for q in [(0,0),(1/3,2/3),(2/3,1/3)]:
  tiny.sector(q);es.append(tiny.mv(np.ones(1,complex))[0].real)
 assert np.max(np.abs(np.sort(es)-[-2.25,-2.25,2.25]))<1e-12
 results.append({'three_site_analytic_spectrum_error':float(np.max(np.abs(np.sort(es)-[-2.25,-2.25,2.25])))})
 polarized=Model(3,27);polarized.sector([0,0]);assert abs(polarized.mv(np.ones(1,complex))[0]-20.25)<1e-12
 results.append({'27_site_fully_polarized_energy':20.25})
 for L in [2,3]:
  c,b=geometry(L);actual={(min(a,z),max(a,z)) for a,z,j in b};ind={(a,z) for a,z,j in independent_bonds(L)}
  assert len(actual)==9*L*L and actual==ind
  degree=np.bincount(np.array(list(actual)).ravel(),minlength=len(c));assert np.all(degree==6)
  results.append({'geometry_L':L,'sites':len(c),'bonds':len(actual),'degree':6,'independent_geometry_match':True})
 for primitive,jp in [(True,1.),(False,.73)]:
  m=Model(2,6,primitive,jp);masks,idx,H=full_h(12,6,independent_bonds(2) if primitive else m.bonds)
  # Enumerate all characters of the translation group, retaining distinct phase lists.
  qs={}
  for a in range(6):
   for b in range(6):
    q=np.array([a/6,b/6]);ph=np.exp(-2j*np.pi*(m.gs@q))
    if all(abs(np.exp(2j*np.pi*np.dot(q,t))-1)<1e-9 for t in [(2,2),(4,-2)]):qs.setdefault(tuple(np.round(ph,9)),q)
  spectrum=[]
  for q in qs.values():
   op=m.sector(q);U=np.zeros((len(masks),m.nb),complex)
   for j,r in enumerate(m.masks):
    for g,perm in enumerate(m.perms):
     t=sum(1<<int(perm[p]) for p in range(12) if int(r)>>p&1)
     U[idx[t],j]=m.phase[g]/np.sqrt(float(m.lengths[j]))
   assert la.norm(U.conj().T@U-np.eye(m.nb))<1e-12
   # Derivation: row H_ij = sqrt(Li/Lj) phase(g) for t=Tg rep_j.
   exact=U.conj().T@(H@U);reduced=np.column_stack([op@np.eye(m.nb,dtype=complex)[:,i] for i in range(m.nb)])
   err=la.norm(reduced-exact);assert err<1e-11,(q,err)
   spectrum.extend(la.eigvalsh(reduced))
   results.append({'primitive':primitive,'Jp':jp,'q':q.tolist(),'dimension':m.nb,'projection_error_frobenius':err})
  assert len(spectrum)==len(masks)
  err=np.max(np.abs(np.sort(spectrum)-la.eigvalsh(H.toarray())))
  assert err<1e-10
  results.append({'primitive':primitive,'complete_spectrum_error':err})
 m=Model(2,6);m.sector([0,0]);r=lanczos(m,maxiter=80,tol=1e-10,check_every=5);assert r['converged']
 Hred=np.column_stack([m.mv(np.eye(m.nb,dtype=complex)[:,i]) for i in range(m.nb)])
 assert abs(r['energy']-la.eigvalsh(Hred)[0])<1e-10
 results.append({'lanczos_small_energy':r['energy'],'residual':r['residual']})
 short=lanczos(m,maxiter=3,tol=1e-12);assert not short['converged'] and short['residual']>1e-12
 rejected=False
 try:m.sector([.25,0])
 except ValueError:rejected=True
 assert rejected
 results.append({'nonconvergence_label_test':True,'incompatible_momentum_rejected':True})
 __import__('pathlib').Path(__file__).with_name('tests.json').write_text(json.dumps(results,indent=2));print('ALL TESTS PASSED')
if __name__=='__main__':run()
