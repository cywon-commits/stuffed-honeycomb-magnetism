import pathlib,json,numpy as np
ROOT=pathlib.Path(__file__).resolve().parent
def read(n):return json.loads((ROOT/n).read_text())
def mat(x):return np.array(x['real'])+1j*np.array(x['imag'])
def run():
 c=[read(f'crossing_{j:.7f}.json') for j in [.3875,.3921875]];w=read('weak_spectrum.json');q=read('weak_quadrature.json');ref=read('moment_reference_validation.json');small=read('small_dynamics_validation.json')
 assert ref['maximum_error']<1e-8 and small['status']=='PASS'
 for d in c:
  assert d['omega']>0 and mat(d['residue_matrix'])[0,0].real>1
  assert d['checks']['total_z_conservation_seed_error']<1e-12 and d['checks']['m_h_opposite_branch_amplitude_max']<1e-8
  for k,v in d['checks'].items():
   if isinstance(v,dict):assert v['H_residual']<1e-8 and v['spin']['S2_eigen_residual']<1e-5
 for d in w['poles']:
  assert d['verified'] and d['H_residual']<1e-8 and d['spin']['S2_eigen_residual']<1e-5 and d['spin']['nearest_half_integer_spin'] in [.5,1.5]
  assert np.linalg.eigvalsh(mat(d['residue_matrix']))[0]>-1e-10
 remainder=[]
 for p in range(3):
  R=mat(w['moment_remainders'][f'M{p}']);remainder.append(float(np.linalg.eigvalsh((R+R.conj().T)/2)[0]));assert remainder[-1]>-1e-8
 assert max(q['moment_errors'])<1e-7 and q['orthogonality_error']<1e-9
 d={'status':'PASS','verified_crossing_poles':2,'verified_weak_poles':len(w['poles']),'reference_moment_max_error':ref['maximum_error'],'quadrature_M0_M1_M2_matrix_norm_errors':q['moment_errors'],'positive_partial_moment_remainder_minimum_eigenvalues_M0_M1_M2':remainder,'interpretation':'Quadrature nodes remain approximate, irrespective of moment closure. No unconverged Ritz value is promoted to a physical pole.'}
 (ROOT/'final_validation.json').write_text(json.dumps(d,indent=2));print(d)
if __name__=='__main__':run()
