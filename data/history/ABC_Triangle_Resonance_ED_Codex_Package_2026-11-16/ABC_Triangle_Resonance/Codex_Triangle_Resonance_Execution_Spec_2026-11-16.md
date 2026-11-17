# CODEX EXECUTION TASK: 27-spin ABC-triangle quantum correlations, J' -> J

**Instruction:** Carry out the missing microscopic ED calculations below; do not stop after reviewing this specification. The supplied ZIP contains exact 27-spin preanalysis, numerical CSV, a figure and a reproducible independent 8x8 spin-1/2 test. Reuse the previously validated 27-spin ED codes/checkpoints. Deliver numerical outputs plus validation logs. Treat incomplete results as incomplete rather than claiming a phase or resonance.

## Scientific question

On the 27-spin, 3x3 periodic stuffed-honeycomb Heisenberg model, the hub and rim labels cease to define distinct bond types at J'=J, where the lattice is isotropic triangular. Is the local doublet (AB,C) <-> (BC,A) <-> (CA,B) a **physically useful trimer-resonance correlation** beyond trivial restoration of lattice symmetry, and what is its relationship to 120-degree magnetic correlations?

Hamiltonian:

`H = J sum_<AB> S_A · S_B + J' sum_<AC,BC> S_i · S_j`; 27 spins: 9 A, 9 B, 9 C/hub. Each A/B rim site has 3 rim and 3 hub neighbors; each C/hub has 6 rim neighbors. 27 AB bonds and 54 C-rim bonds, 54 elementary ABC triangular faces.

**Geometry and normalization MUST be inherited from `ed.geometry(L=3, jp=...)` and `observables.py` supplied in the verified archives.** Unit-cell momentum is not primitive-triangular momentum: primitive physical K folds to the three-site-cell Gamma at J'=J.

## Already completed; please reproduce, not recompute as if missing

The package includes `triangle_27spin_existing.csv`, exact 27x27 SU(2)-scalar two-point matrix snapshots, and `analyze_triangle_existing.py`. Reproduced from the supplied static/crossing ZIPs (11 ground coupling values plus 2 competing branches):

| J'/J | Local doublet P_1/2 | <AB> | <AC>=<BC> | abs(<D>) |
|---:|---:|---:|---:|---:|
| 0.2 | 0.76996325 | -0.37033550 | -0.01730469 | 0.35303081 |
| 0.2609375 | 0.79620800 | -0.35713360 | -0.04358920 | 0.31354440 |
| 0.3921875 | 0.82991759 | -0.32546230 | -0.08470704 | 0.24075526 |
| 0.5 | 0.83199015 | -0.32302122 | -0.08748200 | 0.23553922 |
| 1.0 | 0.87348082 | -0.18674041 | -0.18674041 | < 1e-12 |

For each ABC triangle, `B_ij = S_i · S_j`, `S_tri = S_A+S_B+S_C`, `omega = exp(2 pi i/3)`, and 

`P_1/2 = (15/4 - S_tri^2)/3 = 1/2 - (2/3)(B_AB + B_BC + B_CA)`,

`D_tri = B_AB + omega B_BC + omega^2 B_CA`,

`chi_tri = S_A · (S_B x S_C)` in that *oriented* ordering.

**Critical exact spin-1/2 identities**, independently checked at 8x8 Hilbert-space level in `check_local_identities.py`:

- `chi_tri^2 = (3/16) P_1/2`.
- `(D_tri^dag D_tri + D_tri D_tri^dag)/2 = (9/8) P_1/2`.
- `D_tri^dag D_tri - D_tri D_tri^dag = 3 sqrt(3) chi_tri` (provided D and chi use the same oriented cyclic order).
- Pair singlet probability `P_s(ij) = 1/4 - <B_ij>`.

**Do not treat local <chi^2> or symmetrized <D^dag D> as independent indicators of quantum resonance: they contain no more information than P_1/2.** Equal local bond correlations at J'=J also do not by themselves establish coherent dynamic resonance. The chirality-odd part requires a 3-spin expectation, not deducible from real, two-body static matrices alone.

## REQUIRED tasks, in priority order

**Task A. New J'/J values and ground-state sectors.** Extend verified 27-spin calculations to `J'/J = 0.6, 0.7, 0.8, 0.9` (optionally 0.95). Do not assume that the S_tot=5/2 cell-Gamma branch persists up to J'=1 (which has S_tot=1/2). Screen all inequivalent momentum and symmetry families and report S_tot, energy, residuals and relative branch energy gaps. Track possible additional finite-cluster level crossings in 0.5 < J'/J < 1. At J'=1 use verified primitive K ground-state benchmark and distinguish its cell-Gamma folding from physical primitive Gamma.

**Task B. Local triangle reduced state, doublet flavor and coherence.** At J'/J=0.5, 0.7, 0.9, 1.0 (and if feasible 0.2), construct the oriented ABC-triangle reduced density matrix rho_ABC for the verified ground state (8x8), or compute a complete equivalent set of local spin-rotational scalars. Use a clearly specified point-group and spin-multiplet average and also retain per-eigenstate results if symmetry-breaking/chiral sectors are mixed. Check trace=1, positive semidefiniteness, Hermiticity and projectors. Report doublet/quartet weights; local chirality expectation and variance; and a 2x2 flavor pseudospin matrix in an **orthonormalized** doublet basis. Do not treat the three nonorthogonal singlet-pairing vectors as orthogonal. Quantify chirality-pseudospin polarization, purity, and the limitations of a local mixed state as a 'resonance' diagnostic. Recover exactly the known bond values and P_1/2.

**Task C. Genuine intertriangle connected 4- and 6-spin correlations.** Use oriented triangles t=(x,y,mu) on the 3x3 periodic cell (six faces adjacent to each C site, 54 in total). Orient all faces consistently with the physical triangular plane, e.g. CCW. Compute the complex matrix

`C_D(t,t') = <D_t^dag D_t'> - <D_t^dag><D_t'>`,

and the chirality correlation

`C_chi(t,t') = <chi_t chi_t'> - <chi_t><chi_t'>`.

Report separate statistics for *identical, overlapping,* and *disjoint* triangle pairs so overlapping-site operator algebra cannot masquerade as long-range coherence. Include at least nearest-neighbor and several maximal-distance disjoint triangle pairs. In the case of overlapping triangles, handle operator order exactly (do not silently assume chi_t and chi_t' commute). Provide Hermitian, positive-semidefinite connected correlation matrices. Fourier transform the 6 triangle-orientation channels to 3x3 unit-cell q (9 cell momenta), and at J'=J distinguish physical primitive-momentum components after basis-phase unfolding. Resolve the leading eigenvalues/eigenvectors of dimer and chirality structure factors and compare to conventional physical-K 120-degree spin structure factors computed with established normalization.

**Task D. Optional dynamical confirmation** (after static tasks): Evaluate low-energy resolvents `A_DD(q,omega)` and `A_chichi(q,omega)`, with pole energies/residues and exact zeroth and first moments. Do not infer lifetimes/damping from artificially broadened ED delta peaks.

## Required validation

1. Existing 11 ground and two competitor rows reproduced from archived 27x27 scalar matrices, max absolute error ideally <1e-8. At 0.2, 0.5 and 1, reconstructed energies `27 J <AB> + 54 J' <AC>` must match -10.185949128477, -11.083587020439 and -15.125973385092 respectively.
2. All 54 elementary ABC faces identified, no double counting. Every face has one AB and two A/B-C edges. Confirm 27+54 bonds.
3. Exact local spin identities above evaluated by explicit 8x8 operator algebra and by 27-spin expectations, precision <=1e-8 or document limitation.
4. At isotropy J'=J: AB=BC=CA to available precision; `P_1/2=0.87348082` from existing two-point ED; `|<D>|~0` for symmetry-averaged state. Symmetry-averaged `<chi>` should be checked independently and distinguished from a single chiral momentum-partner state.
5. Compare sublattice and physical primitive momentum conventions. Include spin multiplet/point-group averaging and time-reversal partner treatment; do not lose chirality by averaging before inspecting individual states.
6. Four-/six-spin connected covariance matrices must be Hermitian; covariance eigenvalues nonnegative within numerical tolerance. Cauchy-Schwarz bounds hold; trace and projector tests pass for rho_ABC.
7. Report H eigenvector residuals and S_tot^2 eigen-residuals for every new ground checkpoint. Save full repro scripts and raw per-face/per-pair CSV/NPZ, all momentum-resolved matrices, plots, README, validation log, and a concise physics interpretation with 'established / suggested / not yet established' tiers.

## Available reference archives (from earlier Drive shares)

- Stage 1–2 static ED, 27-spin, contains `hub_static_stages12/static_j*_realspace.npz` and checkpoint code: https://drive.google.com/file/d/1_lx7QD0fpxuq0SO-nrv8lLSZJCW0twek/view?usp=drivesdk
- Level crossings + verified 27x27 real-space matrices for 0.2375–0.425: https://drive.google.com/file/d/1h_GiE26iTnfj2FzQ-y5uBct53eR9OqTA/view?usp=drivesdk
- Isotropic triangular ground verification (physical K): https://drive.google.com/file/d/15V5uOv_KRn-b5fuiEBLrV8OS4i_PE5c3/view?usp=drivesdk

This package includes the small realspace correlation matrices required to independently reproduce all preanalysis **without downloading these large archives**. To compute NEW higher-point ground-state correlations, the original verified eigenvectors and kernels are still needed.

## Outputs

1. `triangle_extended_ground.csv` (0.2–1.0 scanned), `triangle_rdm.csv` (+ rho_ABC files), `intertriangle_dimer_correlators.csv`, `intertriangle_chirality_correlators.csv` and momentum-resolved complex matrices.
2. Convergence/validation log, script(s) rerunnable against raw ED checkpoints, figures suitable for article drafting.
3. A grounded physics discussion: distinguish trivial equivalence of three local bonds, elevated local doublet occupancy, *genuine connected resonance/valence-bond correlations*, and classical/quantum 120-degree magnetic order. Do not call it a distinct phase based on one finite torus.

**Start calculations immediately. If memory/time constraints prohibit part C, complete A/B and the subset of C at J'=J, and report exactly which quantities remain uncomputed.**
