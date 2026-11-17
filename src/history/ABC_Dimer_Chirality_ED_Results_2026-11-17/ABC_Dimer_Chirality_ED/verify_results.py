"""Check measured data, authoritative subspace residues and partial-weight positivity."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,numpy as np,scipy.linalg as la
from experiment import ROOT,dump,geometry_weights
from group_poles import unpack
from ensemble import unfold,fourier

def run():
 responses=json.loads((ROOT/'responses.json').read_text());groups=json.loads((ROOT/'degenerate_pole_groups.json').read_text())['groups'];checks={'moment_matrix_min_eigenvalue':1.,'residue_min_eigenvalue':1.,'remainder_moment_min_eigenvalue':1.,'max_pole_H_residual':0.,'max_pole_S2_residual':0.,'maximum_forbidden_spin_weight':0.,'reference_comparisons':[],'quadrature_moment_errors':{}}
 queue_checks=[]
 for qfile in ['queue_records.json','extension_records.json']:
  if not (ROOT/qfile).exists():continue
  for task in json.loads((ROOT/qfile).read_text()):
   allowed=False
   if task['returncode']!=0 and task['task']=='j09_q_solve' and (ROOT/'finite_solver_restart.json').exists():allowed=json.loads((ROOT/'finite_solver_restart.json').read_text())['replacement_returncode']==0
   assert task['returncode']==0 or allowed,(qfile,task)
   if task['returncode']!=0:queue_checks.append({'task':task['task'],'status':'ABANDONED_NO_POLES; successful exact reflection replacements checked separately'})
 checks['explicit_solver_restarts']=queue_checks
 for f in ROOT.glob('*_moments.json'):
  d=json.loads(f.read_text())
  if 'validation' in d:
   src=d['validation']['source'];assert src['H_residual']<1e-8 and src['spin']['S2_eigen_residual']<1e-7
 for g in groups:
  if g['verified_eigenpair'] and not g['scalar_accessible']:checks['maximum_forbidden_spin_weight']=max(checks['maximum_forbidden_spin_weight'],g['residue_trace'])
 for name,r in responses.items():
  mm=r.get('moments',r);N=len(r['channel_names']);partial=[np.zeros((N,N),complex) for _ in range(3)]
  for p in range(3):
   M=unpack(mm[f'M{p}']);assert la.norm(M-M.conj().T)<1e-7;checks['moment_matrix_min_eigenvalue']=min(checks['moment_matrix_min_eigenvalue'],float(la.eigvalsh((M+M.conj().T)/2)[0]))
  assert la.norm(unpack(mm['Mraw'])-unpack(mm['M0_elastic'])-unpack(mm['M0']))<1e-7
  for pole in r['poles']:
   R=unpack(pole['residue']);assert la.norm(R-R.conj().T)<1e-8;checks['residue_min_eigenvalue']=min(checks['residue_min_eigenvalue'],float(la.eigvalsh(R)[0]));checks['max_pole_H_residual']=max(checks['max_pole_H_residual'],pole['max_direct_H_residual']);checks['max_pole_S2_residual']=max(checks['max_pole_S2_residual'],pole['max_S2_residual'])
   for p in range(3):partial[p]+=pole['omega']**p*R
  for p in range(3):
   rem=unpack(mm[f'M{p}'])-partial[p];checks['remainder_moment_min_eigenvalue']=min(checks['remainder_moment_min_eigenvalue'],float(la.eigvalsh((rem+rem.conj().T)/2)[0]))
 raw=np.load(ROOT/'inputs/raw_observables_j1.npz');oldg=np.load(ROOT/'inputs/isotropic_ground_manifold.npz')
 for name,Q in [('isotropic_G_pair_average',[0,0]),('isotropic_Q_pair_average',[8/9,7/9]),('isotropic_K_pair_average',[1/3,2/3])]:
  if name not in responses:continue
  for kind,a,off in [('D',0,0),('chi',4,108)]:
   for key,C in [('Mraw',raw['ensemble_full_connected']),('M0_elastic',oldg['known_ground_connected']),('M0',oldg['excited_connected'])]:
    data=unfold(fourier(C[off:off+54,off:off+54]),kind)[0];entry=min(data,key=lambda x:la.norm(x['Q']-Q));M=unpack(responses[name][key])[a:a+2,a:a+2];err=float(la.norm(M-entry['matrix']));assert err<1e-6;checks['reference_comparisons'].append({'response':name,'operator':kind,'matrix':key,'error':err})
 for f in ROOT.glob('*_quadrature.json'):
  q=json.loads(f.read_text());checks['quadrature_moment_errors'][q['name']]=q['moment_errors'];assert max(q['moment_errors'])<1e-6
 if 'isotropic_G_pair_average' in responses:
  cross=max([la.norm(unpack(p['residue'])[:4,4:]) for p in responses['isotropic_G_pair_average']['poles']],default=0);checks['Gamma_forbidden_D_chi_cross_max']=float(cross);assert cross<1e-8
 if 'isotropic_Q_pair_average' in responses:checks['finite_Q_allowed_complex_D_chi_cross_max']=float(max([np.max(abs(unpack(p['residue'])[:4,4:])) for p in responses['isotropic_Q_pair_average']['poles']],default=0))
 assert checks['maximum_forbidden_spin_weight']<1e-14 and checks['max_pole_H_residual']<1e-8 and checks['max_pole_S2_residual']<1e-7
 assert min(checks['moment_matrix_min_eigenvalue'],checks['residue_min_eigenvalue'],checks['remainder_moment_min_eigenvalue'])>-1e-7
 checks['status']='PASS';checks['scope']='Directly converged partial poles; unresolved remainders retained. Quadrature nodes are not exported as verified poles. Numerical floating point checks, not certified interval arithmetic.';dump(ROOT/'validation.json',checks);print(json.dumps(checks,indent=2))
if __name__=='__main__':run()
