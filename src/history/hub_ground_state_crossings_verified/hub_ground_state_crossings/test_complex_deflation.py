import pathlib,json,numpy as np,scipy.linalg as la
import cross_scan as cs
ROOT=pathlib.Path(__file__).resolve().parent;cs.ROOT=ROOT/'small_tests';cs.ROOT.mkdir(exist_ok=True)
n=(0,1);r=0;j=.319;m,_,_=cs.get_model(j,n,r,L=3,up=2);H=np.column_stack([m.mv(np.eye(m.nb,dtype=complex)[:,i]) for i in range(m.nb)]);exact=la.eigvalsh(H);records=[]
for idx in range(3):
 d=cs.solve_one(j,n,r,index=idx,ncv=24,seconds=60,force=True,L=3,up=2);assert d['verified'];error=abs(d['energy']-exact[idx]);assert error<1e-10;records.append({'physical_root_index':idx,'eigenvalue_error':float(error),'original_H_residual':d['residual'],'orthogonality':d['orthogonality_to_lower_roots']})
(ROOT/'complex_deflation_validation.json').write_text(json.dumps(records,indent=2));print('COMPLEX PHYSICAL-ROOT DEFLATION TESTS PASSED')
