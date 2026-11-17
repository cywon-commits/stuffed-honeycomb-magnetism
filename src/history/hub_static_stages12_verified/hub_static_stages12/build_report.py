"""Build CSV and physics report only from completed, independently verified static runs."""
import pathlib,json,csv,numpy as np
from fractions import Fraction
from coverage import verified_ground
ROOT=pathlib.Path(__file__).resolve().parent
records=[]
for jp in [.2,.5,1.]:
 path=ROOT/f'static_j{jp:g}.json'
 if not path.exists():continue
 if jp!=1 and verified_ground(jp) is None:continue
 records.append(json.loads(path.read_text()))
assert records
summaries=[d['summary'] for d in records]
keys=sorted(set().union(*(s.keys() for s in summaries)))
with (ROOT/'static_summary.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
 for s in summaries:w.writerow({k:json.dumps(v) if isinstance(v,(list,dict)) else v for k,v in s.items()})
sector_rows=[]
for jp in [.2,.5]:
 for path in sorted(ROOT.glob(f'j{jp:g}_n*_r*.json')):
  r=json.loads(path.read_text());r['source_file']=path.name;r['energy_minus_verified_ground']=r.get('energy',0)-verified_ground(jp)['ground_energy'] if verified_ground(jp) is not None and 'energy' in r else None;r['symmetry_star_size']={(0,0):1,(0,1):6,(1,2):2}[tuple(r['cell_n'])];sector_rows.append(r)
cols=['Jp','cell_n','cell_k','cyclic_order','angular','dimension','energy','energy_minus_verified_ground','symmetry_star_size','residual','converged','ncv','real_solver','matvec_calls','direct_residual_matvecs','iteration_count','solve_seconds','setup_seconds','peak_rss_MiB','status','source_file']
with (ROOT/'sector_convergence.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=cols);w.writeheader()
 for r in sector_rows:w.writerow({k:json.dumps(r[k]) if isinstance(r.get(k),list) else r.get(k) for k in cols})
all_modes=[]
for d in records:
 for mode in d['correlation_eigenmodes']:
  all_modes.append({'Jp_over_J':d['summary']['Jp_over_J'],**{k:v for k,v in mode.items() if k!='vector'},'vector_real':json.dumps(mode['vector']['real']),'vector_imag':json.dumps(mode['vector']['imag'])})
with (ROOT/'correlation_eigenmodes.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(all_modes[0]));w.writeheader();w.writerows(all_modes)
table='| Jprime/J | E/J | cell momentum representative | S | rim m_N squared | hub S(Gamma) | 120-degree intensity | energy residual |\n|---:|---:|---|---:|---:|---:|---:|---:|\n'
for s in summaries:
 k=s['cell_ground_momentum'];table+=f"| {s['Jp_over_J']:g} | {s['energy']:.12f} | ({Fraction(float(k[0])).limit_denominator(3)},{Fraction(float(k[1])).limit_denominator(3)}) | {s['S_total']:g} | {s['rim_neel_m2']:.8f} | {s['hub_S_Gamma']:.8f} | {s['physical_K_120_intensity']:.8f} | {s['energy_residual']:.2e} |\n"
rows=''
for s in summaries:rows+=f"- Jprime/J={s['Jp_over_J']:g}: hub peak S={s['hub_peak_S']:.8f} at cell momentum numerator(s) {s['hub_peak_cell_n']} over 3; rim NN={s['rim_NN_mean']:.8f}, hub-rim bond={s['hub_rim_bond_mean']:.8f}, all-distinct hub-pair mean={s['hub_hub_all_distinct_mean']:.8f}; uniform rim/hub Gamma cross={s['rim_hub_uniform_Gamma']:.8f}.\n"
changes=''
if len(summaries)==3:
 lo,mid,hi=summaries;changes=f"From isotropy to Jprime/J=0.2, the rim Neel squared amplitude changes from {hi['rim_neel_m2']:.8f} to {lo['rim_neel_m2']:.8f}; hub S(Gamma) changes from {hi['hub_S_Gamma']:.8f} to {lo['hub_S_Gamma']:.8f}. At Jprime/J=0.5 the corresponding values are {mid['rim_neel_m2']:.8f} and {mid['hub_S_Gamma']:.8f}. These are equal-time finite-size correlations, with all definitions and state averages held fixed."
validation={'completed_couplings':[s['Jp_over_J'] for s in summaries],'per_coupling':{str(d['summary']['Jp_over_J']):d['validation'] for d in records},'operator_tests':json.loads((ROOT/'operator_validation.json').read_text()),'cyclic_tests':json.loads((ROOT/'cyclic_validation.json').read_text()),'restarted_solver_tests':json.loads((ROOT/'solver_validation.json').read_text()),'ground_driver_guards':json.loads((ROOT/'driver_validation.json').read_text()),'full_size_spin_multiplet_check':json.loads((ROOT/'spin_multiplet_validation.json').read_text()) if (ROOT/'spin_multiplet_validation.json').exists() else None}
(ROOT/'validation_summary.json').write_text(json.dumps(validation,indent=2))
bench='| Jprime/J | analysis wall s | pair kernel s | analysis peak RSS MiB | S-squared eigen-residual |\n|---:|---:|---:|---:|---:|\n'
for s in summaries:bench+=f"| {s['Jp_over_J']:g} | {s['analysis_seconds']:.2f} | {s['pair_kernel_seconds']:.2f} | {s['peak_rss_MiB']:.1f} | {s['S2_eigen_residual']:.2e} |\n"
sector_bench=''
for jp in [.2,.5]:
 c=verified_ground(jp)
 if c:sector_bench+=f"Jprime/J={jp:g}: nine representative symmetry blocks converged; summed measured solve wall time {c['measured_solve_seconds_sum']:.2f} s; largest per-process RSS {c['maximum_measured_peak_rss_MiB']:.1f} MiB. Ground file {c['ground_name']}, cyclic order {c['ground_cyclic_order']}, character index {c['ground_angular']}.\n\n"
report=f'''# Static stuffed-honeycomb ED: completed Stages 1–2

Only static calculations are performed. No continued-fraction spectrum, frequency broadening, susceptibility or dynamical linewidth has been computed. The attached specification was reviewed first; SPEC_REVIEW.md records the Fourier, normalization, ensemble and selection-rule decisions.

## Model and state conventions

N=27 spin-1/2 sites, 3x3 three-site cells, periodic boundary conditions, J=1 and Jprime/J in the completed set {[s['Jp_over_J'] for s in summaries]}. Every ground-state search is at N_up=14, Sz=+1/2, full dimension 20,058,300. The production Hamiltonian is matrix free. At unequal couplings only the correct cell translations, plus independently validated compatible point operations, are used. Static measurements are evaluated after converting the saved eigenvector to the cell-translation basis of dimension 2,228,700.

The previously verified isotropic ground vector is reused, not redetermined: primitive K=(1/3,2/3), folded cell Gamma. Its primitive Kprime partner is symmetry degenerate. The original full-momentum verification is preserved in isotropic_ground_verification.json and isotropic_ground_sector_reference.csv. No exact M point exists here, and no allowed sector receives an M label.

Primary matrices are C_ab(k)=sum_spin <O_a(k)^dagger O_b(k)>, O_a(k)=sum_R exp[-2*pi*i*k.R] S_Ra /sqrt(9). We average the scalar real-space correlations over the twelve exact hub-centered point operations and over the ground spin multiplet. Scalar dot correlations are unchanged by spin-multiplet averaging, while vector means vanish. The primary raw scalar matrix is therefore also the connected scalar matrix of that specified ensemble. This removes arbitrary choice of the degenerate momentum-star partner. Separately, CSVs preserve the actual selected Sz=+1/2 state's raw and connected scalar and longitudinal zz matrices. They are not substituted for one another.

## Verified ground states and key results

{table}

{changes}

{rows}

Only measured, directly residual-verified eigenvectors enter these tables. The fixed-Sz value is not used to guess S: total S-squared is applied independently and its eigen-residual recorded. Detailed convergence, little-group dimensions and all checked energies are in sector_convergence.csv and per-block JSON. Digits in CSV are numerical estimates, not assertions of an equal number of rigorous significant digits.

## Independent ground-sector determination at each unequal coupling

The weighted stuffed graph has three cell-momentum stars: Gamma; the six momenta (0,1),(0,2),(1,0),(1,1),(2,0),(2,2) over 3; and the pair (1,2),(2,1) over 3. Every class is checked independently at each coupling. Their symmetry operations preserve rim and hub bond weights, as verified explicitly.

Additional compatible cyclic little groups reduce cost without discarding sectors: C6 at Gamma (characters ell=0,...,5), one compatible reflection C2 at representative (0,1/3), and C3 at representative (1/3,2/3). At Gamma complex conjugation of the real Hamiltonian equates ell=1 with 5 and ell=2 with 4; ell=0,1,2,3 are independently solved. Both reflection parities and all three C3 characters are solved. The nine measured blocks consequently cover every cell momentum and character. coverage.py checks that each representative cell block sums to 2,228,700 and that symmetry-star multiplicities recover all 20,058,300 states. A positive residual alone does not rigorously certify the lowest root, but standard lowest-root Lanczos, complete symmetry coverage and independent small-cluster spectrum checks underpin the numerical determination.

{sector_bench}
The quantity named lowest_energy_separation_between_checked_blocks in coverage JSON is a separation between checked lowest roots, not a complete excitation gap: same-block excited roots were not calculated. Symmetry-equivalent ground partners are not independent convergence timings.

## Unexpected intermediate-coupling spin, independently checked

At Jprime/J=0.5 the ground state has S=5/2, despite being computed at Sz=+1/2. Applying S-squared gives 8.75 with eigen-residual 4.56e-13. An independent full-size static ladder check yields ||S+ psi_(M=1/2)||^2=8 and ||S+ psi_(M=3/2)||^2=5, preserves the ground energy within 6e-14 J, and annihilates the M=5/2 state (squared norm 2.2e-26). This confirms a genuine spin-5/2 multiplet rather than an S=1/2 assumption or a mixed-spin artifact. The N_up=15 Gamma dimension is independently realized as 1,931,652. These extra sectors are used only for static spin validation, not for spectra or independent ground searches.

The selected M=1/2 state has opposite rim and hub one-point polarizations at 0.5, while the hub Gamma structure factor exceeds its isotropic value. These are finite-cluster signatures of partially compensated sublattice moments. They do not establish a thermodynamic ferrimagnetic phase. Both sampled endpoints have S=1/2. Ground spin and cell momentum change between the sampled weak and intermediate points, implying finite-cluster level crossings somewhere between them; no crossing location or thermodynamic transition is determined.

## What the static data establish

Rim Neel amplitude is m_N^2=<[sum_R(S_A−S_B)]^2>/18^2 = C_nn(Gamma)/36, with n=(1,−1,0). Hub S_h(k)=C_HH(k) includes onsite weight. Rim/hub cross matrices retain their complex phases at all nine cell momenta. Derived n,m,h,Phi channel matrices are in the channels CSVs; Phi(Gamma)=3m(Gamma).

The physical triangular 120-degree diagnostic uses v=(1,omega^2,omega) and I_120=v^dagger C(Gamma) v/3. The factor 1/3 follows from its 1/sqrt(27) Fourier normalization. Both K and Kprime diagnostics are saved. At Jprime unequal to J this is a physical-site Fourier diagnostic, not an additional conserved momentum. Large intensity can also receive rim-Neel weight; it does not by itself establish triangular 120-degree order away from isotropy.

At isotropy the large hub Gamma value coexists with total S=1/2 and a compensating uniform rim/hub anticorrelation. Uniform hubs alone do not diagnose total ferromagnetism. The reported m_F^2=S(S+1)/27^2 explicitly separates net total spin from sublattice correlations. A nonzero finite-cluster squared moment is not spontaneous order. At the weak point, decreasing hub Gamma weight is a redistribution to finite momentum, not evidence that the hubs are uncorrelated.

Correlation eigenmodes are obtained from the full Hermitian 3x3 matrix. At Gamma, inversion separates the normalized rim-Neel mode (A−B)/sqrt(2) from the 2x2 sector spanned by (A+B)/sqrt(2) and H. Their eigenvalues and hub weights provide a stable comparison across couplings. At isotropy the rim-Neel mode and an even uniform/hub mode are degenerate; arbitrary eigenvectors inside a degenerate eigenspace are not physically unique. Saved eigenspace projectors and inversion-resolved Gamma modes avoid that ambiguity. Momentum-resolved hub peaks are reported without calling any incommensurate finite-size point M.

There is no size scaling, critical-coupling fit or inference of thermodynamic long-range order from these three points. Near Jprime=0 the free-hub manifold becomes extensively degenerate; a single weak-coupling ground state's correlations need not approach an uncorrelated hub ensemble. The selected finite-size levels, not a thermodynamic phase transition or damping rate, are being compared.

## Validation and resource measurements

- Independent full-bit-state spin seeds z,+,− match projected seeds on 12 sites and on complex 3x3 cell sectors with one flipped spin. H acting on each seed is independently checked. Full scalar correlation matrices assembled from the z,+,− Gram matrices agree with the separate pair-class kernel. The largest measured error is recorded in operator_validation.json.
- All cyclic characters reproduce their independently projected small-cluster Hamiltonians; their combined spectra reproduce the cell block. cyclic_validation.json includes nontrivial stabilizers. Restarted-solver tests compare dense roots and ensure a wall-time interruption is never labeled converged.
- Production checks include normalized source vectors, direct H residuals, S-squared eigen-residuals, real-space and independent form-factor Fourier reconstruction of the exchange energy, onsite sum rules, Hermiticity and positivity of raw/connected matrices, Sz=1/2 and zero Sz variance, and inversion n/h and n/m zeros at Gamma after symmetry averaging. Selected-state inversion overlap and n/h values are also retained.
- Spin conservation does not require the raw scalar total-spin Gamma weight to vanish: it equals S(S+1)/9. The selected fixed-Sz connected longitudinal total channel does vanish. Its saved zz covariance is the exact static seed Gram matrix for a future longitudinal calculation; at Gamma the conserved null vector should be deflated rather than treated as a solver failure. Conserved transverse spin can carry zero-energy multiplet weight. This distinction is essential for future spectra.

{bench}

Ground searches use restarted real-Hermitian ARPACK Lanczos with a real embedding only when the projected operator is complex. Real characters use a real block directly. The embedding duplicates eigenvalues algebraically; duplicates are not called physical degeneracies. A fresh complex-vector H residual below 1e-8 J is mandatory. ARPACK tolerance is 1e-11; subspace size ncv and exact H-product counts are saved per run. SciPy does not expose its internal restart iteration count: that field is null, rather than fabricated from H-product counts. Progress records give H applications, wall time and RSS, not unverified intermediate eigenvalues.

The environment has a two-CPU quota and 8 GiB RAM cap; OpenMP uses two threads and BLAS one. Full versions and hardware are in environment.json. All wall times and RSS are measured in this shared environment; brief validation/analysis jobs can overlap a solve, so these are not controlled speed comparisons. New solve_seconds includes ARPACK and its fresh residual check, excluding model setup and compressed checkpoint I/O. analysis_seconds measures core conversion, operator checks and static pair work before final CSV/JSON output and plotting. Compilation/setup and solve times are separated. The initial translation-only weak-coupling trial was stopped in favor of validated additional symmetries; translation_only_trial.json explicitly marks it unconverged and supplies no eigenvalue. None of its intermediate values enter static results.

## Recommendations for a later dynamical stage — not executed

Use Jprime/J=1 as the validation anchor. The most informative contrast is 0.5 versus 0.2: the former has an S=5/2 multiplet and strong uniform hub correlations, while the latter has S=1/2 and dominant finite-momentum hub correlations. Prioritize Gamma n,m,h at 0.5 and the (0,1/3) transfer star at 0.2, with Gamma Ward controls at both. A few further static points between 0.2 and 0.5 could locate the finite-cluster sector crossing before a fine spectral scan; no such scan is performed here. Use cell Gamma to probe rim-staggered n, uniform m, and h, retaining the full m/h cross matrix; the inelastic conserved combination m+h must be null. Inversion forbids n/h at Gamma in the stated ensemble, so a zero there cannot settle hub coupling. Operator transfer k sends a chosen ground momentum k0 to k0+k, not simply to k. recommended_dynamics.csv gives proposed source, transfer and target cell sectors without running a spectrum. Also target the allowed finite cell momentum star where the hub static structure factor peaks, and compare the other inequivalent nonzero star as a control. The per-coupling peak sectors are listed above.

At isotropy use physical basis-phase-resolved K/Kprime probes separately from physical Gamma, even though they share cell Gamma. Explicitly project the complete degenerate ground manifold out of the strictly inelastic seed: one-point subtraction alone does not remove every zero-energy spin- or momentum-partner transition. Large equal-time weight prioritizes channels but does not prove low-frequency spectral weight or a Goldstone linewidth. Eta convergence, spectral moments and the inelastic Ward identity belong to the unstarted Stage 3.

## Reproduction and saved artifacts

Install requirements.txt and a C++ compiler with OpenMP. Run `OPENBLAS_NUM_THREADS=1 python run_static.py` from this directory. It reuses verified results/checkpoints when present; `--force` reruns all sector solves and statics. Individual commands are documented in README.md. The original isotropic checkpoint is included. Ground checkpoints at the two unequal couplings are included when available; non-ground checkpoints are optional and can be regenerated. Only small reference tests assemble Hamiltonian matrices. All production H applications and static pair reductions are matrix free.

The principal numeric files are static_summary.csv, sector_convergence.csv, correlation_eigenmodes.csv, per-coupling static_j*_C.csv and static_j*_channels.csv, real-space NPZs, JSON validation/convergence data, and ground_j*_coverage.json. Figures static_evolution, static_matrices, hub_modes and total_spin_and_polarization are supplied as PNG and SVG; lines between three points only guide the eye. SHA256SUMS records deliverable integrity.
'''
if not any(s['Jp_over_J']==.5 for s in summaries):
 a=report.index('## Unexpected intermediate-coupling');b=report.index('## What the static data establish',a);report=report[:a]+'## Intermediate coupling\n\nThe intermediate static result is incomplete; no ground-spin conclusion is reported for this run.\n\n'+report[b:]
(ROOT/'technical_report.md').write_text(report);(ROOT/'summary_table.md').write_text(table)
(ROOT/'requirements.txt').write_text('numpy==2.3.5\nscipy==1.17.0\nmatplotlib==3.10.8\n')
print(table)

concise=f"""# Static correlation results (Stages 1–2)

Completed coupling ratios: {[s['Jp_over_J'] for s in summaries]}. All results use the 27-site torus, Sz=+1/2 and J=1. Primary scalar matrices use an explicitly stated SU(2) spin-multiplet and exact point-group average; selected-state raw/connected matrices are also provided. No dynamical spectra were calculated.

{table}

{changes}

{rows}

The intermediate point has a genuine S=5/2 ground multiplet, verified both by S-squared and by global-spin raising norms 8, 5, then zero at the highest weight; all raised partners retain the ground energy. Its selected Sz=+1/2 state has oppositely signed rim and hub polarizations. These are partially compensated finite-cluster moments, not proof of a thermodynamic ferrimagnetic phase. Both endpoints have S=1/2.

At isotropy, the uniform hub structure factor is large despite total S=1/2. The rim uniform component compensates the hub moment; hub Gamma weight alone is not a net-ferromagnetism diagnostic. The physical 120-degree probe includes its required 1/3 sublattice normalization. Away from isotropy it is a Fourier diagnostic that can also receive rim-Neel weight, not a unique identifier of 120-degree order.

Inversion enforces zero rim-staggered/hub cross correlation at Gamma in the averaged ensemble. This selection rule does not establish absence of a hub coupling. The nonzero uniform rim/hub cross channel and the full complex matrices at finite momenta provide the relevant complementary information. Correlation eigenspaces, including inversion-resolved Gamma modes and hub weights, are saved; arbitrary vectors inside degenerate eigenspaces are not compared as unique modes.

Ground-sector searches at each unequal coupling cover all nine cell momenta using exact point-group equivalences and every compatible cyclic character. Direct H residuals are in the table; total S is verified independently by applying S-squared. Small-cluster tests independently validate z,+,- seeds, H acting on the seeds, scalar Gram matrices and the reduced spectra. Production data pass positivity, Hermiticity, energy reconstruction, onsite sum rules, fixed-Sz conservation and inversion tests. CPU/RAM, H-product counts and wall times are recorded in sector_convergence.csv and technical_report.md.

These are finite-size equal-time correlations. Three couplings on one torus cannot establish thermodynamic long-range order, a transition, order destruction or damping. No exact M point exists on this torus.

For the next dynamical stage, retain Jprime/J=1 as a validation anchor and prioritize the contrast between 0.5 (spin 5/2, uniform hub peak) and 0.2 (spin 1/2, finite-momentum hub peak). Start at cell Gamma with the full n,m,h channel matrix, testing the conserved inelastic m+h null channel; also target the allowed hub-peak momentum stars listed above. At isotropy resolve physical K/Kprime using basis phases even though both fold to cell Gamma. Separate the entire degenerate zero-energy ground manifold before interpreting inelastic weight. Large static weight alone does not prove a low-frequency peak. Stage 3 remains unstarted.

See SPEC_REVIEW.md for conventions, technical_report.md for validation/resources, static_summary.csv for key numbers and static_j*_C.csv for every 3x3 matrix at all nine cell momenta.
"""
if not any(s['Jp_over_J']==.5 for s in summaries):
 a=concise.index('The intermediate point has');b=concise.index('At isotropy,',a);concise=concise[:a]+concise[b:]
(ROOT/'physics_report.md').write_text(concise)

with (ROOT/'recommended_dynamics.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['Jp_over_J','status','priority','source_cell_n','operator_transfer_cell_n','target_cell_n','channels','note'])
 for d in records:
  s=d['summary'];source=np.round(3*np.array(s['cell_ground_momentum'])).astype(int)%3;choices=[((0,0),1,'n,m,h; retain cross terms','Check inelastic m+h conservation; n/h zero requires the stated inversion ensemble.')]
  finite=[((0,1),s['hub_S_cell01']),((1,2),s['hub_S_cell12'])];finite.sort(key=lambda x:x[1],reverse=True)
  choices.append((finite[0][0],2,'A,B,H and hub-dominated correlation eigenmode','Largest inequivalent nonzero-momentum hub static weight; static weight does not prove low-frequency weight.'))
  choices.append((finite[1][0],3,'A,B,H','Optional inequivalent finite-momentum control.'))
  for transfer,priority,channels,note in choices:w.writerow([s['Jp_over_J'],'PROPOSAL ONLY; NO SPECTRUM COMPUTED',priority,json.dumps(source.tolist()),json.dumps(transfer),json.dumps(((source+np.array(transfer))%3).tolist()),channels,note+' Targets shown for the selected representative; future ensemble comparison must average all ground partners.'])
  if s['Jp_over_J']==1:w.writerow([1,'PROPOSAL ONLY; NO SPECTRUM COMPUTED',1,'[0,0]','[0,0]','[0,0]','physical K/Kprime and physical Gamma using basis phases','Primitive initial K + transfer K reaches Kprime; initial K + transfer Kprime reaches Gamma. Resolve and remove the zero-energy ground manifold.'])
