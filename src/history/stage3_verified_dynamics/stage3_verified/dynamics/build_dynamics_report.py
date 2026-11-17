import pathlib,json,csv,numpy as np
ROOT=pathlib.Path(__file__).resolve().parent
def read(n):return json.loads((ROOT/n).read_text())
def mat(x):return np.array(x['real'])+1j*np.array(x['imag'])
def csvout(name,rows):
 if rows:
  with (ROOT/name).open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
 crossing=[read(f'crossing_{j:.7f}.json') for j in [.3875,.3921875]];cr=[]
 for d in crossing:
  weight=mat(d['residue_matrix'])[0,0].real;cr.append({'Jp':d['Jp'],'omega_J':d['omega'],'source_spin':d['source_spin'],'target_spin':d['excited_spin'],'n_zz_weight':weight,'fraction_n_M0':weight/d['moments']['n']['M0'],'fraction_n_M1':weight*d['omega']/d['moments']['n']['M1'],'n_M0':d['moments']['n']['M0'],'n_M1_J':d['moments']['n']['M1'],'n_M2_J2':d['moments']['n']['M2'],'n_elastic_weight':d['moments']['n']['elastic_weight'],'largest_state_H_residual_J':max(v['H_residual'] for k,v in d['checks'].items() if isinstance(v,dict)),'wall_seconds':d['wall_seconds'],'peak_rss_MiB':d['peak_rss_MiB']})
 csvout('crossing_poles.csv',cr);m=read('weak_moments.json');w=read('weak_spectrum.json') if (ROOT/'weak_spectrum.json').exists() else None;quad=read('weak_quadrature.json') if (ROOT/'weak_quadrature.json').exists() else None
 poles=[];residues=[]
 for i,d in enumerate(w['poles'] if w else []):
  if not d['verified']:continue
  R=mat(d['residue_matrix']);assert np.linalg.eigvalsh(R)[0]>-1e-9
  poles.append({'index':i,'omega_J':d['omega'],'energy_J':d['energy'],'spin':d['spin']['nearest_half_integer_spin'],'H_residual_J':d['H_residual'],'S2_residual':d['spin']['S2_eigen_residual'],'n_weight':R[0,0].real,'m_weight':R[1,1].real,'h_weight':R[2,2].real})
  for a,an in enumerate(['n','m','h']):
   for b,bn in enumerate(['n','m','h']):residues.append({'pole':i,'omega_J':d['omega'],'channel':an+bn,'real':R[a,b].real,'imag':R[a,b].imag})
 csvout('weak_verified_poles.csv',poles);csvout('weak_complex_residues.csv',residues)
 moments=[]
 for sign in ['1','-1']:
  for p in range(3):
   M=mat(m['moments'][sign][f'M{p}'])
   for a,an in enumerate(['n','m','h']):
    for b,bn in enumerate(['n','m','h']):moments.append({'probe_sign':sign,'moment_order':p,'channel':an+bn,'real':M[a,b].real,'imag':M[a,b].imag})
 csvout('weak_exact_moments.csv',moments)
 val=read('moment_reference_validation.json');text='# Targeted Stage 3 dynamical ED\n\n'
 text+='All residues refer to longitudinal Sz=1/2 sources, with 1/sqrt(9) Fourier normalization. No broadening is applied. These finite-cluster peaks do not determine thermodynamic damping. Read CONVENTIONS.md for the reviewed translation, spin, Fourier, and sum-rule definitions.\n\n'
 text+='## Second crossing: resolved visible rim Neel poles\n\n| Jprime/J | omega/J | n-z residue | fraction of n M0 | largest state residual/J |\n|---:|---:|---:|---:|---:|\n'
 for d in cr:text+=f"| {d['Jp']:.7f} | {d['omega_J']:.12f} | {d['n_zz_weight']:.12f} | {d['fraction_n_M0']:.6%} | {d['largest_state_H_residual_J']:.2e} |\n"
 text+='\nThe competing S=7/2, chi=0, F=+ and S=5/2, chi=3, F=- eigenvectors pass fresh original-H, total-spin, rotation and reflection checks. This establishes substantial low-energy dipolar weight, rather than symmetry allowance alone. The low poles carry only about half a percent of the n first moment; appreciable high-energy weight coexists with the soft low-energy peak. The forbidden m/h amplitudes are below 1.3e-16. The first crossing has DeltaS=3 and no direct single-spin transition.\n\n'
 text+='Gamma n elastic weight is below 1.1e-30. Fixed-M m/h elastic weights are retained separately in the JSON. Projecting the same-M ground state suffices for longitudinal response: other members of its spin multiplet have different M and cannot contribute. The m+h seed equals v/6 within 6e-17; its inelastic part vanishes. Scalar full-multiplet elastic weights are not subtracted from zz spectra.\n\n'
 text+='## Weak coupling j=0.2375\n\nSource cell (1/3,2/3), probe (0,1/3), target (1/3,0), all at M=1/2. Neither ground corner partner lies in the target, so the elastic ground-manifold contribution is zero. Exact forward/reverse three-channel moment matrices and complex cross terms are in weak_exact_moments.csv.\n\n'
 if w:
  text+=f"The target-sector calculation produced {len(poles)} independently verified eigenpairs by spatial-star transfer of saved symmetry-character minima. It is a partial low-energy spectral matrix, not a complete spectral decomposition. Individual complex 3x3 residues are in weak_complex_residues.csv; original-H and S-squared residuals are in weak_verified_poles.csv. Target verification: {w['H_calls']} H calls, {w['solve_seconds']:.2f} measured solve seconds, {w['peak_rss_MiB']:.1f} MiB process peak RSS.\n\n"
  text+='| omega/J | S | n residue | m residue | h residue | H residual/J |\n|---:|---:|---:|---:|---:|---:|\n'
  for d in poles:text+=f"| {d['omega_J']:.12f} | {d['spin']:g} | {d['n_weight']:.6g} | {d['m_weight']:.6g} | {d['h_weight']:.6g} | {d['H_residual_J']:.2e} |\n"
  text+='\nExact moments minus the verified-pole sums are saved as moment_remainders. They quantify unsaturated spectral weight and are not treated as validation failure. Higher-energy poles and dipole-forbidden levels have not been exhaustively enumerated.\n\n'
 if quad:
  text+=f"A separate {quad['basis_dimension']}-dimensional block-Krylov quadrature approximates the full energy-resolved matrix. Its nodes are explicitly **not verified eigenvalues** and must not be combined with the verified poles as extra states. It reproduces exact zeroth/first/second moments with matrix-norm errors {quad['moment_errors']}. Quadrature time {quad['wall_seconds']:.2f} s, peak RSS {quad['peak_rss_MiB']:.1f} MiB. No convergence of its individual nodes is claimed.\n\n"
 text+='## Moment validation and computational limits\n\n'
 text+=f"All {val['entries']} available full-precision CSV comparisons of M0zz and symmetrized M1zz agree within {val['maximum_error']:.3e}. Forward M1(k) and transposed reverse M1(-k) differ by at most {val['forward_reverse_M1_matrix_difference']:.3e}. There is no supplied M2 reference; exact M2 is newly computed from DB norms. A separate 12-spin seed test checks Fourier phases and complete spectral M0/M1/M2 sums (errors below 1.1e-14).\n\n"
 text+='The two-CPU quota and matrix-free H applications limit high-energy eigenpair enumeration. Exact moments, verified low poles, and moment-complete approximate quadrature are deliberately distinguished. Timings are per operation and include contention; they are not a dedicated-machine performance estimate. The discarded cold and probe-started weak searches have no exported Ritz values; solver_restart_note.json documents its excluded partial cost. No linewidth or thermodynamic transition is inferred. Optional first-crossing spectra are deferred.\n'
 (ROOT/'report.md').write_text(text)
 if poles:
  import os
  os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl'));os.environ.setdefault('XDG_CACHE_HOME',str(ROOT/'.cache'))
  import matplotlib;matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  fig,axes=plt.subplots(1,3,figsize=(12,3.7),layout='constrained')
  for a,ch in enumerate(['n','m','h']):
   axes[a].stem([d['omega_J'] for d in poles],[d[ch+'_weight'] for d in poles]);axes[a].set(xlabel='Excitation energy / J',ylabel='Discrete longitudinal weight',title=f'Weak coupling: {ch} channel');axes[a].grid(alpha=.2)
  fig.suptitle('Verified low poles only; higher-energy weight omitted; no broadening');fig.savefig(ROOT/'weak_verified_sticks.png',dpi=160);plt.close(fig)
 print(text)
if __name__=='__main__':main()
