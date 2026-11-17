"""Write an evidence-based finite-cluster crossing report; no spectral calculation."""
import pathlib,json,csv,numpy as np
from fractions import Fraction
from export_results import build
ROOT=pathlib.Path(__file__).resolve().parent

def qtext(q):return '('+','.join(str(Fraction(n,3)) for n in q)+')'
def main():
 ground,sectors,obs,intervals=build();coarse=[d for d in ground if d['Jp_over_J'] in [.275,.35,.425]]
 table='| Jprime/J | E0/J | cell momentum | S | rotation character | residual/J |\n|---:|---:|---|---:|---:|---:|\n'
 for d in coarse:table+=f"| {d['Jp_over_J']:g} | {d['energy_J']:.12f} | {qtext(d['cell_n'])} | {d['spin']:g} | {d['character']} | {d['direct_residual_J']:.2e} |\n"
 cross_table='| Crossing bracket Jprime/J | Width | Left quantum numbers | Right quantum numbers | DeltaE at left / J | DeltaE at right / J |\n|---|---:|---|---|---:|---:|\n'
 selection=[]
 if intervals:
  for i,c in enumerate(intervals['crossings'],1):
   A=c['left_quantum_numbers'];B=c['right_quantum_numbers'];la=f'{qtext(A[0])}, chi={A[1]}, S={A[2]:g}, F={A[3]}';lb=f'{qtext(B[0])}, chi={B[1]}, S={B[2]:g}, F={B[3]}';ga=c.get('gap_right_minus_left_at_a');gb=c.get('gap_right_minus_left_at_b')
   cross_table+=f"| [{c['a']:.9f}, {c['b']:.9f}] | {c['width']:.7f} | {la} | {lb} | {ga if ga is not None else 'not computed'} | {gb if gb is not None else 'not computed'} |\n"
   rank1=abs(A[2]-B[2])<=1 and A[2]+B[2]>=1;bothGamma=A[0]==B[0]==[0,0];parity_ok=A[3] is not None and B[3] is not None;n_allowed=rank1 and bothGamma and (B[1]-A[1])%6==3 and (not parity_ok or B[3]==-A[3]);mh_allowed=rank1 and bothGamma and A[1]==B[1] and (not parity_ok or B[3]==A[3]);selection.append({'crossing':i,'bracket':[c['a'],c['b']],'delta_S':abs(A[2]-B[2]),'single_spin_rank_one_allowed':rank1,'required_transfer_cell_n':[(b-a)%3 for a,b in zip(A[0],B[0])],'n_Gamma_symmetry_allowed':n_allowed,'m_h_Gamma_symmetry_allowed':mh_allowed,'scope':'Analytic selection rules only; transition matrix elements and spectral weights have NOT been calculated.'})
 (ROOT/'selection_rules.json').write_text(json.dumps(selection,indent=2))
 cost=json.loads((ROOT/'computational_cost.json').read_text());proofs=json.loads((ROOT/'interval_exclusions.json').read_text()) if (ROOT/'interval_exclusions.json').exists() else []
 unexcluded=[d for d in proofs if not d['excluded_other_families']];(ROOT/'remaining_interval_ambiguities.json').write_text(json.dumps(unexcluded,indent=2))
 st='| Jprime/J | Role | S | rim m_N squared | hub S(Gamma) | hub S(0,1/3) | hub S(1/3,2/3) |\n|---:|---|---:|---:|---:|---:|---:|\n'
 for d in obs:
  if d['role'].startswith('GROUND'):st+=f"| {d['Jp']:.9g} | {d['role']} | {d['spin']:g} | {d['rim_neel_m2']:.7f} | {d['hub_S_Gamma']:.7f} | {d['hub_S_cell01']:.7f} | {d['hub_S_cell12']:.7f} |\n"
 # Three actual, verified sample points spanning the observed sectors.
 choices=[.2375,.35,.425];recs=[]
 for j in choices:
  point=next((d for d in ground if abs(d['Jp_over_J']-j)<1e-10),None)
  if point is None:continue
  if j==.2375:channels='hub at transfer (0,1/3), with Gamma n,m,h controls';reason='Weak, finite-momentum hub correlations; below the first observed crossing.'
  elif j==.35:channels='Gamma n, and separate m/h channels';reason='Interior S=7/2 state, below the higher-coupling crossing.'
  else:channels='Gamma n, and separate m/h channels';reason='Interior S=5/2 state, above the higher-coupling crossing.'
  recs.append({'Jp_over_J':j,'spin':point['spin'],'source_cell_n':point['cell_n'],'source_name':point['name'],'channels':channels,'reason':reason,'status':'RECOMMENDATION ONLY; STAGE 3 NOT STARTED'})
 with (ROOT/'stage3_recommendations.csv').open('w') as f:
  if recs:
   writer=csv.DictWriter(f,fieldnames=list(recs[0]));writer.writeheader()
   for r in recs:writer.writerow({k:json.dumps(v) if isinstance(v,list) else v for k,v in r.items()})
 status='Refinement complete to the requested bracket width.' if intervals and intervals['all_bracket_widths_met'] else 'Refinement is incomplete; only completed samples and any existing brackets are reported.'
 rules=''
 for r in selection:
  if not r['single_spin_rank_one_allowed']:rules+=f"Crossing {r['crossing']} has DeltaS={r['delta_S']:g}: the two crossing levels cannot couple through a single-spin vector probe, at any momentum. Closing that interbranch gap therefore does not by itself imply a dipolar soft pole. Other low-energy states are not excluded by this rule.\n\n"
  elif r['n_Gamma_symmetry_allowed']:rules+=f"Crossing {r['crossing']} is symmetry-allowed in the rim-staggered n(Gamma) channel (DeltaS=1 and opposite rotation/reflection characters). The direct pair transition is forbidden for uniform m/h(Gamma). Allowed does not mean its spectral weight has been measured.\n\n"
 contrasts=list(csv.DictReader((ROOT/'same_coupling_branch_contrasts.csv').open())) if (ROOT/'same_coupling_branch_contrasts.csv').exists() else []
 contrast_text=''
 for c in contrasts:
  contrast_text+=f"At the common coupling {float(c['common_Jp']):g}, comparing the left branch with its competing right branch changes m_N^2 from {float(c['rim_neel_m2_left']):.7f} to {float(c['rim_neel_m2_right']):.7f}, and S_h(Gamma) from {float(c['hub_S_Gamma_left']):.7f} to {float(c['hub_S_Gamma_right']):.7f}. This isolates branch dependence from changing the coupling; the right branch is an excited eigenstate at that point.\n\n"
 contrast_text+='The saved subsystem moments satisfy 324*m_N^2=2*(<SA^2>+<SB^2>)-<SR^2> and S*(S+1)=<SR^2>+<SH^2>+2*<SR dot SH>. Subsystem moments are expectation values, not conserved subsystem-spin labels. The first branch contrast is dominated by hub polarization and increased rim collective spin; the second has a larger rim-staggered reduction and stronger antialignment of rim and hub moments. These are finite-cluster correlations.'
 text=f'''# 27-spin stuffed-honeycomb ground-state level crossings

{status} All ground-state searches use J=1, N=27, N_up=14, Sz=+1/2. These are finite-cluster level crossings, not thermodynamic phase transitions. Stage 3 has not started.

## Requested initial points

{table}

The S=7/2 Gamma ground state is an additional intermediate sector, absent from the supplied endpoint labels. Every half-integer total-spin multiplet has an Sz=+1/2 member, so the fixed-Sz ground search includes higher-spin ground states. S is measured, not inferred from Sz.

`cell_sector_comparison.csv` gives the three inequivalent momentum-star minima (or explicitly labeled bounds at refinement points); the star multiplicities are 1,6,2. At the three initial points all cell-sector minima are identified. `character_evidence.csv` distinguishes measured eigenvalues from lower bounds, and `corner_character_equivalence.json` records the exact reflection relating the two complex C3 characters. Cell momenta are n/3 in the three-site convention. No momentum is labeled M; there is no exact triangular M point on this torus.

## Refined observed ground-sector changes

{cross_table}

DeltaE means E(right-hand branch) minus E(left-hand branch), evaluated at the same coupling. Its sign reverses across a verified bracket. These are competing-level splittings, not a complete first excitation gap. Same-block excited roots not explicitly calculated, and exact spin/momentum ground partners, must be distinguished. A linear-interpolation estimate in `crossing_brackets.csv` is only an estimate; the measured bracket is authoritative.

The target width is approximately 0.005 in Jprime/J. Eigensolver residuals are reported separately from coupling-bracket uncertainty. `ground_samples.csv` lists all classified samples; `all_levels.csv` retains the individual physical roots and Casimir residuals. Rotation character chi labels C6 at Gamma (eigenvalue exp[2*pi*i*chi/6]); F is the hub-centered reflection parity where it is defined. Reflection and total-spin labels prevent mistaken continuation through a quantum-number change.

## Static checks around the crossings

{st}

All scalar correlations use the same normalization as the previous static study: m_N^2=C_nn(Gamma)/36 and S_h(k)=C_HH(k), with spin-multiplet and exact point-group averaging. The observable columns are transfer momenta, not the initial ground momentum. Competing non-ground eigenstates at a common coupling are explicitly labeled in the full static CSV. Their contrasts are not mislabeled ground-state measurements. Real-space matrices and full sublattice matrices are saved for each measured state.

The first observed sector change accompanies a strong redistribution from finite-momentum hub weight to a nearly uniform hub subsystem. The S=7/2 point at 0.275 has S_h(Gamma)=2.7047294, while the verified weak endpoint has 0.0951725. Rim Neel correlations remain appreciable on both sides. The higher-coupling change reduces total spin to 5/2 while uniform hub correlations remain strong. These finite-size observations do not establish spontaneous order, a phase boundary, a linewidth or damping.

{contrast_text}

## Validation, branch tracking and remaining limits

Hamiltonian-vector products are matrix free in the verified translation/cyclic-little-group bases. Cold generic pseudorandom starts (seed 2026) are used for lowest roots; they do not lock the solver to the endpoint spin. Direct original-H residuals below 1e-8 J and direct S-squared eigen-residuals below 1e-5 are required; actual errors are in the numeric tables. The new spin-7/2 state was independently verified by global raising norms 15,12,7,0 and equal-energy partners. That auxiliary static SU(2) check uses other M components only for verification; all ground searches stay at Sz=+1/2.

Additional roots are obtained by physical-vector projector deflation and checked for orthogonality, original-H residual and total spin. The shift exceeds the full spectral-width bound, so previously found roots are raised above the physical spectrum. Complex real-embedding duplicates are not counted as new physical states. Small-cluster tests independently verify deflation, spectra, concavity and affine Rayleigh bounds. `branch_symmetry_and_overlaps.json`, when present, records same-quantum-number overlaps across couplings; the figures separate spin/spatial labels rather than connecting unlike states as one adiabatic branch.

For a fixed spatial symmetry block H(j)=H0+jV, its lowest energy is concave. Chords between lower endpoint energy estimates provide numerical lower bounds; each verified normalized state's E(j0)+(j-j0)<V> is a variational upper line. These bounds exclude higher characters without assigning them fictitious eigenvalues. Exclusions are conditional on the numerically verified lowest-root identification, not formal interval-arithmetic certificates. `interval_exclusions.json` records continuous spatial-family checks. {len(unexcluded)} stable sampled intervals retain at least one unexcluded alternative spatial family; their ranges are listed in `remaining_interval_ambiguities.json`. `bracket_third_family_exclusions.json` additionally excludes third spatial families throughout both crossing brackets using those conditional bounds (minimum margins approximately 0.000205 J and 0.02685 J). Same-family spin islands narrower than the audit spacing are not rigorously excluded. The result is a targeted crossing analysis, not an exhaustive analytic phase diagram.

Measured new work: {cost['new_verified_eigenpairs']} verified eigenpairs, {cost['measured_H_calls_sum']} H applications, {cost['measured_solve_wall_seconds_sum']/60:.2f} summed solve minutes plus {cost['measured_setup_seconds_sum']:.2f} setup seconds, process peak RSS {cost['process_peak_rss_MiB']:.1f} MiB. Hardware quota: two CPUs, 8 GiB RAM; two OpenMP threads and one BLAS thread. These sums exclude reused endpoint solves, uninstrumented endpoint profiling, plotting, and separately recorded static analyses. Cached bases make RSS a process high-water measure, not an isolated per-block allocation. No unconverged Ritz value enters a ground-state or crossing table. `computational_cost.json` gives exact cost definitions and any incomplete-pair counts.

## Three representative Stage 3 couplings

Recommend Jprime/J = **0.2375, 0.35, 0.425**, provided the corresponding verified rows remain present in the delivered table. These are interior representative samples, rather than arbitrary states mixed at an exact crossing. They cover the weak low-spin regime, the spin-7/2 intermediate regime, and the spin-5/2 regime; 0.35 and 0.425 lie on opposite sides of the higher-coupling crossing. `stage3_recommendations.csv` records source sectors and proposed channels.

{rules}

At 0.2375 prioritize the allowed (0,1/3) transfer star for hub response, with Gamma controls. At 0.35 and 0.425 prioritize separate Gamma n and m/h probes, keeping cross terms and testing the conserved inelastic m+h null channel. For every calculation the final cell momentum is k0+k_transfer. Project the complete degenerate ground manifold out before interpreting inelastic weight. These are recommendations only; no transition matrix element, susceptibility or dynamical spectrum was calculated here.
'''
 (ROOT/'crossing_report.md').write_text(text);(ROOT/'initial_points_table.md').write_text(table);print(table);print(cross_table)
if __name__=='__main__':main()
