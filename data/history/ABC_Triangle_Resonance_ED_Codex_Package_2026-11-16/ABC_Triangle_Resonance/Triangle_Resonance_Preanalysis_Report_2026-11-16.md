# ABC triangle correlations from verified 27-spin ED — 2026-11-16

## Question and result

Does an isotropic J'=J limit admit a useful local picture `(AB,C) <-> (BC,A) <-> (CA,B)` where an elementary ABC triangle occupies the spin-1/2 doublet sector with alternative spin-pair singlets? **Yes as a local basis picture**, but the present two-point ED results do **not** independently establish coherent resonating-bond dynamics or a separate quantum-disordered phase.

The quantitative result is that the probability that an elementary triangle has total spin 1/2 rises from 0.769963 at J'/J=0.2 to 0.873481 at J'/J=1, and the local singlet probabilities on AB and AC/BC become equal as symmetry is restored. Across two known 27-spin finite-size total-spin crossings, the triangle doublet probability also jumps.

## Geometry and source

Hamiltonian `H = J sum_AB S_A·S_B + J' sum_A-C,B-C S_i·S_j` on 27 spin-1/2 sites, 3x3 three-site cell, periodic. Exact geometry: 27 AB edges, 54 AC/BC edges, 54 elementary ABC triangles. Observables derive from `G_pointgroup` (the previously verified point-group averaged SU(2)-scalar spin correlation matrix) in `hub_static_stages12_verified.zip` (j=0.2,0.5,1) and `hub_ground_state_crossings_verified.zip` (additional 8 ground values, 2 competitor branches). Correlation matrices are reused, not re-diagonalized. The code reconstructs the Hamiltonian energy from 27 and 54 bond counts and validates no face double-counting.

## Exact formulas

Define `B_ij=S_i·S_j`, `omega=exp(2pi i/3)` and `S_triangle=S_A+S_B+S_C`. For three spin-1/2's,

`S_triangle² = 9/4 + 2(B_AB+B_BC+B_CA)`.

`P_1/2 = (15/4 - S_triangle²)/3 = 1/2 - (2/3)(B_AB+B_BC+B_CA)`.

`P_singlet(ij) = 1/4 - <B_ij>`.

`D_triangle = B_AB + omega B_BC + omega² B_CA`.

`chi_triangle = S_A · (S_B x S_C)` with consistent orientation.

Exact matrix identities, independently checked by an 8x8 Pauli-matrix construction:

`chi²=(3/16)P_1/2`;

`(D†D+DD†)/2=(9/8)P_1/2`;

`D†D - DD†=3 sqrt(3) chi`.

Therefore `<chi²>` and the *symmetrized* local dimer norm add no independent information beyond the local doublet probability. A nonzero `<chi>` (not determined by two-point scalar correlations) is required to separate `<D†D>` from `<DD†>` in a non-TR-invariant state.

## ED-derived numerical results (ground states)

| J'/J | P_1/2 | <B_AB> | <B_AC> = <B_BC> | |<D>| | <chi²> |
|---:|---:|---:|---:|---:|---:|
| 0.2000000 | 0.76996325 | -0.37033550 | -0.01730469 | 0.35303081 | 0.14436811 |
| 0.2375000 | 0.77401699 | -0.36863095 | -0.02119727 | 0.34743369 | 0.14512819 |
| 0.2562500 | 0.77606798 | -0.36762220 | -0.02323988 | 0.34438232 | 0.14551275 |
| 0.2609375 | 0.79620800 | -0.35713360 | -0.04358920 | 0.31354440 | 0.14928900 |
| 0.2750000 | 0.79657887 | -0.35692999 | -0.04396916 | 0.31296083 | 0.14935854 |
| 0.3500000 | 0.79830822 | -0.35575491 | -0.04585371 | 0.30990120 | 0.14968279 |
| 0.3875000 | 0.79903952 | -0.35511453 | -0.04672238 | 0.30839215 | 0.14981991 |
| 0.3921875 | 0.82991759 | -0.32546230 | -0.08470704 | 0.24075526 | 0.15560955 |
| 0.4250000 | 0.83077045 | -0.32458235 | -0.08578666 | 0.23879569 | 0.15576946 |
| 0.5000000 | 0.83199015 | -0.32302122 | -0.08748200 | 0.23553922 | 0.15599815 |
| 1.0000000 | 0.87348082 | -0.18674041 | -0.18674041 | ~0 | 0.16377765 |

The finite-size crossing intervals are [0.25625,0.2609375] and [0.3875,0.3921875]. The competitor state at 0.25625 has `P_1/2 = 0.79608069` (versus 0.77606798 for ground), and the competitor at 0.3875 has `P_1/2=0.82975537` (versus 0.79903952 ground), showing branch-specific local correlation changes. Avoid assuming the exact coupling of the crossing is known more precisely than the bracket.

At J'=J the probability of an AB two-spin singlet is `0.43674041`; the AC/BC singlet probabilities are identical. The 3-spin doublet probability 0.87348082 is *not* the sum of the three two-spin singlet probabilities, because the three pair-singlet bases are nonorthogonal. The two physical doublets together span a four-dimensional state space (spin-1/2 pair times two flavor/chirality states).

The reference classical product state with three spin-1/2 vectors 120 degrees apart has `<B_ij> = -1/8` and `P_1/2 = 3/4`; the ED value 0.87348 is higher, consistent with strong quantum antiferromagnetic local correlations. This does not by itself demonstrate a resonating-valence-bond phase.

## Interpretation and next discriminating tests

1. **Established:** local triangle doublet weight increases on the sampled ground branches; original AB-versus-C distinction in the two-point bond correlators disappears at isotropy.
2. **Established:** both total-spin crossings have a local triangle-doublet step. They cannot be identified as thermodynamic transitions from this finite cluster.
3. **Not yet established:** time-dependent or spatially coherent intertriangle bond resonance, spontaneous chiral order, a new phase, or suppression of 120-degree long-range order.
4. **Needed:** triangle reduced density matrices including chirality and doublet flavor, *connected intertriangle* dimer-dimer 4-spin correlators, scalar chirality-chirality 6-spin correlators, symmetry-preserving Fourier structure factors, and an extension to 0.6–0.9 with full-sector ground-state verification. The attached Codex execution specification gives precise tasks and numerical validation anchors.

## Source-level validation

- Numerical code and complete 11 ground + two competitor CSV included.
- 8x8 identity checks pass to norms of order 1e-16.
- Maximum dispersion of any given AB or AC/BC bond scalar among all 54 elementary triangles in point-group-averaged data is <= 2e-16.
- Reconstructed energies: E(J'/J=0.2) = -10.185949128477..., E(0.5) = -11.083587020439..., E(1.0) = -15.125973385092... consistent with verified ED references.

For reproducibility use `python check_local_identities.py` and `python analyze_triangle_existing.py` (uses supplied small matrix snapshots without large archives, or direct original archives when present).
