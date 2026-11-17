import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,hashlib,json,numpy as np,scipy.linalg as la
from experiment import ROOT,geometry_weights,pack
from ensemble import unfold,fourier
R=ROOT/'reference/ABC_Dimer_Chirality_Dynamics_2026-11-16';manifest=[]
for line in (R/'SHA256SUMS').read_text().splitlines():
 h,p=line.split('  ',1);actual=hashlib.sha256((R/p).read_bytes()).hexdigest();manifest.append({'file':p,'match':actual==h});assert h==actual
raw=np.load(ROOT/'inputs/raw_observables_j1.npz');ground=np.load(ROOT/'inputs/isotropic_ground_manifold.npz');rows=[]
swap18=np.r_[np.arange(6,12),np.arange(6),np.arange(12,18)];sg18=np.r_[np.ones(12),-np.ones(6)]
for p in sorted(ROOT.glob('j1*_moments.json')):
 r=json.loads(p.read_text());name=r['name'];Q=np.array(r['Q']);n=np.rint(np.mod([Q[0]+Q[1],2*Q[0]-Q[1]],1)*3).astype(int)%3;W=geometry_weights(Q)
 if r['source_partner']:
  nn=(-n)%3;G=raw['selected_operator_Gram_q'][3*nn[0]+nn[1]];G=G[np.ix_(swap18,swap18)].conj()*sg18[:,None]*sg18[None,:]
 else:G=raw['selected_operator_Gram_q'][3*n[0]+n[1]]
 expected=W.conj().T@G@W;got=np.array(r['Mraw']['real'])+1j*np.array(r['Mraw']['imag']);err=la.norm(got-expected);assert err<1e-6;rows.append({'experiment':name,'selected_static_M0_error':float(err)})
for prefix,j in [('j09',.9),('j07',.7)]:
 for p in sorted(ROOT.glob(prefix+'*_moments.json')):
  r=json.loads(p.read_text());old=np.load(ROOT/f'inputs/raw_observables_j{j:g}.npz');n=np.rint(np.asarray(r['Q'])*3).astype(int)%3;expected=old['selected_operator_Gram_q'][3*n[0]+n[1]];got=np.array(r['Mraw']['real'])+1j*np.array(r['Mraw']['imag']);err=float(la.norm(got-expected));assert err<1e-8;rows.append({'experiment':r['name'],'selected_static_M0_error':err})
# All-q local trace sums retain overlapping order and both D circular senses.
C=raw['selected_full_raw_Gram'];rho=raw['selected_rho_ABC'];from local_geometry import chi,P
local=[]
for mu in range(6):
 ch=float(((-1)**mu*np.trace(rho[mu]@chi)).real);prob=float(np.trace(rho[mu]@P).real);delta=C[mu,mu]-C[54+mu,54+mu];sym=(C[mu,mu]+C[54+mu,54+mu])/2;ch2=C[108+mu,108+mu]
 errors=[abs(delta-3*np.sqrt(3)*ch),abs(sym-9*prob/8),abs(ch2-3*prob/16)];assert max(errors)<1e-10;local.append({'mu':mu,'chi_CCW':ch,'D_circular_integrated_asymmetry':float(delta.real),'P_doublet':prob,'errors':[float(x) for x in errors]})
result={'manifest':manifest,'reference_M0':rows,'local_sum_rules_from_all_q_archived_Gram_including_zero_frequency':local,'notes':'Supplied 27-spin reference specifies M0; no independent 27-spin M1/M2 reference numbers are present. New M1/M2 use actual H actions and are checked by Krylov moment reproduction and independent small-cluster spectral sums.'};(ROOT/'input_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps({'manifest_matches':len(manifest),'M0':rows,'maximum_local_sum_rule_error':max(max(r['errors']) for r in local)}))
