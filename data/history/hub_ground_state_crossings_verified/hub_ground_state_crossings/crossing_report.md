# 27-spin stuffed-honeycomb ground-state level crossings

Refinement complete to the requested bracket width. All ground-state searches use J=1, N=27, N_up=14, Sz=+1/2. These are finite-cluster level crossings, not thermodynamic phase transitions. Stage 3 has not started.

## Requested initial points

| Jprime/J | E0/J | cell momentum | S | rotation character | residual/J |
|---:|---:|---|---:|---:|---:|
| 0.275 | -10.290051743981 | (0,0) | 3.5 | 0 | 1.78e-14 |
| 0.35 | -10.472017603235 | (0,0) | 3.5 | 0 | 5.00e-11 |
| 0.425 | -10.732527324629 | (0,0) | 2.5 | 3 | 2.78e-13 |


The S=7/2 Gamma ground state is an additional intermediate sector, absent from the supplied endpoint labels. Every half-integer total-spin multiplet has an Sz=+1/2 member, so the fixed-Sz ground search includes higher-spin ground states. S is measured, not inferred from Sz.

`cell_sector_comparison.csv` gives the three inequivalent momentum-star minima (or explicitly labeled bounds at refinement points); the star multiplicities are 1,6,2. At the three initial points all cell-sector minima are identified. `character_evidence.csv` distinguishes measured eigenvalues from lower bounds, and `corner_character_equivalence.json` records the exact reflection relating the two complex C3 characters. Cell momenta are n/3 in the three-site convention. No momentum is labeled M; there is no exact triangular M point on this torus.

## Refined observed ground-sector changes

| Crossing bracket Jprime/J | Width | Left quantum numbers | Right quantum numbers | DeltaE at left / J | DeltaE at right / J |
|---|---:|---|---|---:|---:|
| [0.256250000, 0.260937500] | 0.0046875 | (1/3,2/3), chi=0, S=0.5, F=1 | (0,0), chi=0, S=3.5, F=1 | 0.0015922604671203544 | -0.0034764877275517136 |
| [0.387500000, 0.392187500] | 0.0046875 | (0,0), chi=0, S=3.5, F=1 | (0,0), chi=3, S=2.5, F=-1 | 0.005756116930829336 | -0.003820478183513387 |


DeltaE means E(right-hand branch) minus E(left-hand branch), evaluated at the same coupling. Its sign reverses across a verified bracket. These are competing-level splittings, not a complete first excitation gap. Same-block excited roots not explicitly calculated, and exact spin/momentum ground partners, must be distinguished. A linear-interpolation estimate in `crossing_brackets.csv` is only an estimate; the measured bracket is authoritative.

The target width is approximately 0.005 in Jprime/J. Eigensolver residuals are reported separately from coupling-bracket uncertainty. `ground_samples.csv` lists all classified samples; `all_levels.csv` retains the individual physical roots and Casimir residuals. Rotation character chi labels C6 at Gamma (eigenvalue exp[2*pi*i*chi/6]); F is the hub-centered reflection parity where it is defined. Reflection and total-spin labels prevent mistaken continuation through a quantum-number change.

## Static checks around the crossings

| Jprime/J | Role | S | rim m_N squared | hub S(Gamma) | hub S(0,1/3) | hub S(1/3,2/3) |
|---:|---|---:|---:|---:|---:|---:|
| 0.2 | GROUND_REUSED_ENDPOINT | 0.5 | 0.2453222 | 0.0951725 | 0.9953117 | 0.3414788 |
| 0.2375 | GROUND | 0.5 | 0.2439606 | 0.1008686 | 0.9931179 | 0.3452121 |
| 0.25625 | GROUND | 0.5 | 0.2431458 | 0.1042740 | 0.9918225 | 0.3473954 |
| 0.2609375 | GROUND | 3.5 | 0.2391022 | 2.7084050 | 0.5067597 | 0.5005183 |
| 0.275 | GROUND | 3.5 | 0.2389333 | 2.7047294 | 0.5073565 | 0.5005658 |
| 0.35 | GROUND | 3.5 | 0.2379818 | 2.6838011 | 0.5107525 | 0.5008419 |
| 0.3875 | GROUND | 3.5 | 0.2374801 | 2.6726211 | 0.5125655 | 0.5009929 |
| 0.3921875 | GROUND | 2.5 | 0.2244896 | 2.5532664 | 0.5258202 | 0.5209062 |
| 0.425 | GROUND | 2.5 | 0.2238387 | 2.5551533 | 0.5278251 | 0.5139480 |
| 0.5 | GROUND_REUSED_ENDPOINT | 2.5 | 0.2225708 | 2.5372596 | 0.5326034 | 0.5085599 |


All scalar correlations use the same normalization as the previous static study: m_N^2=C_nn(Gamma)/36 and S_h(k)=C_HH(k), with spin-multiplet and exact point-group averaging. The observable columns are transfer momenta, not the initial ground momentum. Competing non-ground eigenstates at a common coupling are explicitly labeled in the full static CSV. Their contrasts are not mislabeled ground-state measurements. Real-space matrices and full sublattice matrices are saved for each measured state.

The first observed sector change accompanies a strong redistribution from finite-momentum hub weight to a nearly uniform hub subsystem. The S=7/2 point at 0.275 has S_h(Gamma)=2.7047294, while the verified weak endpoint has 0.0951725. Rim Neel correlations remain appreciable on both sides. The higher-coupling change reduces total spin to 5/2 while uniform hub correlations remain strong. These finite-size observations do not establish spontaneous order, a phase boundary, a linewidth or damping.

At the common coupling 0.25625, comparing the left branch with its competing right branch changes m_N^2 from 0.2431458 to 0.2391577, and S_h(Gamma) from 0.1042740 to 2.7096104. This isolates branch dependence from changing the coupling; the right branch is an excited eigenstate at that point.

At the common coupling 0.3875, comparing the left branch with its competing right branch changes m_N^2 from 0.2374801 to 0.2245959, and S_h(Gamma) from 2.6726211 to 2.5516495. This isolates branch dependence from changing the coupling; the right branch is an excited eigenstate at that point.

The saved subsystem moments satisfy 324*m_N^2=2*(<SA^2>+<SB^2>)-<SR^2> and S*(S+1)=<SR^2>+<SH^2>+2*<SR dot SH>. Subsystem moments are expectation values, not conserved subsystem-spin labels. The first branch contrast is dominated by hub polarization and increased rim collective spin; the second has a larger rim-staggered reduction and stronger antialignment of rim and hub moments. These are finite-cluster correlations.

## Validation, branch tracking and remaining limits

Hamiltonian-vector products are matrix free in the verified translation/cyclic-little-group bases. Cold generic pseudorandom starts (seed 2026) are used for lowest roots; they do not lock the solver to the endpoint spin. Direct original-H residuals below 1e-8 J and direct S-squared eigen-residuals below 1e-5 are required; actual errors are in the numeric tables. The new spin-7/2 state was independently verified by global raising norms 15,12,7,0 and equal-energy partners. That auxiliary static SU(2) check uses other M components only for verification; all ground searches stay at Sz=+1/2.

Additional roots are obtained by physical-vector projector deflation and checked for orthogonality, original-H residual and total spin. The shift exceeds the full spectral-width bound, so previously found roots are raised above the physical spectrum. Complex real-embedding duplicates are not counted as new physical states. Small-cluster tests independently verify deflation, spectra, concavity and affine Rayleigh bounds. `branch_symmetry_and_overlaps.json`, when present, records same-quantum-number overlaps across couplings; the figures separate spin/spatial labels rather than connecting unlike states as one adiabatic branch.

For a fixed spatial symmetry block H(j)=H0+jV, its lowest energy is concave. Chords between lower endpoint energy estimates provide numerical lower bounds; each verified normalized state's E(j0)+(j-j0)<V> is a variational upper line. These bounds exclude higher characters without assigning them fictitious eigenvalues. Exclusions are conditional on the numerically verified lowest-root identification, not formal interval-arithmetic certificates. `interval_exclusions.json` records continuous spatial-family checks. 3 stable sampled intervals retain at least one unexcluded alternative spatial family; their ranges are listed in `remaining_interval_ambiguities.json`. `bracket_third_family_exclusions.json` additionally excludes third spatial families throughout both crossing brackets using those conditional bounds (minimum margins approximately 0.000205 J and 0.02685 J). Same-family spin islands narrower than the audit spacing are not rigorously excluded. The result is a targeted crossing analysis, not an exhaustive analytic phase diagram.

Measured new work: 54 verified eigenpairs, 9924 H applications, 120.06 summed solve minutes plus 70.11 setup seconds, process peak RSS 2302.8 MiB. Hardware quota: two CPUs, 8 GiB RAM; two OpenMP threads and one BLAS thread. These sums exclude reused endpoint solves, uninstrumented endpoint profiling, plotting, and separately recorded static analyses. Cached bases make RSS a process high-water measure, not an isolated per-block allocation. No unconverged Ritz value enters a ground-state or crossing table. `computational_cost.json` gives exact cost definitions and any incomplete-pair counts.

## Three representative Stage 3 couplings

Recommend Jprime/J = **0.2375, 0.35, 0.425**, provided the corresponding verified rows remain present in the delivered table. These are interior representative samples, rather than arbitrary states mixed at an exact crossing. They cover the weak low-spin regime, the spin-7/2 intermediate regime, and the spin-5/2 regime; 0.35 and 0.425 lie on opposite sides of the higher-coupling crossing. `stage3_recommendations.csv` records source sectors and proposed channels.

Crossing 1 has DeltaS=3: the two crossing levels cannot couple through a single-spin vector probe, at any momentum. Closing that interbranch gap therefore does not by itself imply a dipolar soft pole. Other low-energy states are not excluded by this rule.

Crossing 2 is symmetry-allowed in the rim-staggered n(Gamma) channel (DeltaS=1 and opposite rotation/reflection characters). The direct pair transition is forbidden for uniform m/h(Gamma). Allowed does not mean its spectral weight has been measured.



At 0.2375 prioritize the allowed (0,1/3) transfer star for hub response, with Gamma controls. At 0.35 and 0.425 prioritize separate Gamma n and m/h probes, keeping cross terms and testing the conserved inelastic m+h null channel. For every calculation the final cell momentum is k0+k_transfer. Project the complete degenerate ground manifold out before interpreting inelastic weight. These are recommendations only; no transition matrix element, susceptibility or dynamical spectrum was calculated here.
