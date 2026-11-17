"""Check the correct conjugate source/target pair, not just +/- probes of one source."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,numpy as np,scipy.linalg as la
from cyclic_model import CyclicModel
from observables import site_weights,coefficients
ROOT=pathlib.Path(__file__).resolve().parent
out=[]
for n,partner,r in [((1,2),(2,1),0),((1,0),(2,0),0),((1,0),(2,0),1)]:
 a=CyclicModel(jp=.2375,cell_n=n,angular=r);b=CyclicModel(jp=.2375,cell_n=partner,angular=r);assert np.array_equal(a.masks,b.masks) and np.array_equal(a.lengths,b.lengths) and np.array_equal(a.rawbasis,b.rawbasis);phase_error=float(np.max(abs(a.phase.conj()-b.phase)));rng=np.random.default_rng(2026);v=rng.normal(size=a.nb)+1j*rng.normal(size=a.nb);v/=la.norm(v);err=float(la.norm(b.mv(v.conj())-a.mv(v).conj()));assert phase_error<1e-12 and err<1e-10;out.append({'cell_n':n,'conjugate_cell_n':partner,'reflection_character':r,'phase_conjugation_error':phase_error,'H_conjugation_error':err});del a,b,v
operr=max(float(np.max(abs(site_weights(3,[0,-1/3],coefficients([0,-1/3],ch))-site_weights(3,[0,1/3],coefficients([0,1/3],ch)).conj()))) for ch in ['n','m','H','Phi']);assert operr<1e-12
(ROOT/'partner_source_validation.json').write_text(json.dumps({'status':'PASS','blocks':out,'operator_conjugation_error':operr,'identity':'The conjugate source at -k0 probed at -k has target -(k0+k) and moment/residue matrices conjugate to the original forward matrices. This is distinct from reversing only the probe around the original k0 source. Real bond couplings make the identity valid at both couplings.'},indent=2));print(out)
