import pathlib,json,numpy as np,scipy.linalg as la
import cross_scan as cs
ROOT=pathlib.Path(__file__).resolve().parent;cs.ROOT=ROOT/'small_tests';cs.ROOT.mkdir(exist_ok=True)
records=[]
for n,r in [((0,0),0),((0,0),3)]:
 m,_,_=cs.get_model(.371,n,r,L=2,up=6);H=np.column_stack([m.mv(np.eye(m.nb,dtype=complex)[:,i]) for i in range(m.nb)]);exact=la.eigvalsh(H)
 for idx in range(3):
  rr=cs.solve_one(.371,n,r,index=idx,ncv=24,seconds=60,force=True,L=2,up=6);assert rr['verified'];error=abs(rr['energy']-exact[idx]);assert error<1e-10;records.append({'character':r,'root':idx,'energy_error':error,'residual':rr['residual'],'S2_residual':rr['spin_check']['S2_eigen_residual']})
 # Check affine Rayleigh upper lines and eigenvalue chord lower bounds against exact spectra.
 values=[]
 for j in [.2,.35,.5]:
  mm,_,_=cs.get_model(j,n,r,L=2,up=6);HH=np.column_stack([mm.mv(np.eye(mm.nb,dtype=complex)[:,i]) for i in range(mm.nb)]);e,U=la.eigh(HH);values.append(e[0])
  if j==.2:seed=U[:,0];E=e[0];derivative=cs.slope(mm,seed)
 assert values[1] >=(values[0]+values[2])/2-1e-11
 assert values[1] <= E+.15*derivative+1e-11
cs.ROOT=ROOT;cs.certify_reflection()
(ROOT/'crossing_validation.json').write_text(json.dumps({'deflation_tests':records,'concavity_and_variational_line_tests':True,'corner_reflection_equivalence_test':True},indent=2));print('CROSSING TESTS PASSED')
