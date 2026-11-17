import numpy as np,scipy.linalg as la,json,pathlib
from cyclic_model import CyclicModel
from ed import Model
from test_observables import projector
from observables import convert_basis
out=[]
for L,up,ns in [(2,6,[(0,0)]),(3,1,[(0,0),(0,1),(1,2)])]:
 for celln in ns:
  cm=Model(L,up,False,.2);k=np.array(celln)/L;q=[(k[0]+k[1])/3,(2*k[0]-k[1])/3];cm.sector(q);Uc,_,_=projector(cm);Hc=np.column_stack([cm.mv(np.eye(cm.nb,dtype=complex)[:,i]) for i in range(cm.nb)])
  values=[];dims=[];errs=[];order=CyclicModel(L,up,.2,celln,0).order
  for ell in range(order):
   m=CyclicModel(L,up,.2,celln,ell);U,_,_=projector(m)
   H=np.column_stack([m.mv(np.eye(m.nb,dtype=complex)[:,i]) for i in range(m.nb)]) if m.nb else np.zeros((0,0),complex)
   P=Uc.conj().T@U
   if m.nb:
    vv=np.arange(1,m.nb+1,dtype=complex);vv/=la.norm(vv);assert la.norm(convert_basis(m,cm,vv)-P@vv)<1e-12
   err=la.norm(H-P.conj().T@Hc@P);assert err<1e-10
   values.extend(la.eigvalsh(H));dims.append(m.nb);errs.append(err)
  err=np.max(abs(np.sort(values)-la.eigvalsh(Hc)));assert err<1e-10 and sum(dims)==cm.nb
  out.append({'L':L,'up':up,'cell_n':celln,'cyclic_order':order,'dimensions':dims,'spectrum_error':float(err),'projection_error':max(errs)})
pathlib.Path(__file__).with_name('cyclic_validation.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2));print('CYCLIC LITTLE-GROUP TESTS PASSED')
