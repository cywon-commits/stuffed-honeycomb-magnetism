# Codex execution request: finite-size validation of hub polarization and rim-spin compensation

## Goal
Carry out calculations, not only a feasibility review. Test whether the mechanisms observed in the 27-spin spin-1/2 stuffed-honeycomb antiferromagnetic Heisenberg model survive changes in cluster size. Focus on `r=J'/J < 0.7` and separate well-supported finite-size results from thermodynamic extrapolations. **Run the affordable 18-spin phase first; treat 36 spins as a conditional HPC-scale phase.**

## 1. Model and conventions

Use the exact Hamiltonian and geometry from the existing *verified* Codex 27-spin repository, with `J=1`:

`H = J sum_<AB rim bonds> S_i·S_j + J' sum_<hub-rim bonds> T_R·S_i`.

There are 3 spin-1/2 sites per unit cell: rim A, rim B, hub H; each hub has six neighboring rim sites, each rim site has three neighboring hub sites. Avoid double-counting bonds. Do not replace the model with a different dice-lattice convention. Maintain periodic boundary conditions and the same reciprocal-space conventions.

The earlier 27-spin calculation (`Ncell=9`, `N=27`) found ground-state total-spin sequence, in sampled regimes, `1/2 -> 7/2 -> 5/2 -> 3/2 -> 1/2` as r rises from weak coupling to about 0.7. The first two 27-spin crossings were bracketed around ~0.25772 and ~0.39032. These are **finite-size level crossings**, not confirmed bulk transitions. Validate reference data from actual checkpoints/CSV rather than using approximate numbers from prose.

## 2. Geometry and momentum prerequisite: mandatory

For each candidate supercell defined by integer lattice vectors `(a,b)`, `(c,d)` in the 3-site-cell triangular Bravais basis, verify:

- `Ncell = |ad-bc|` and `Nspin = 3*Ncell`;
- periodic bonds, hub/rim coordination, unique bond count and translational automorphisms;
- admissibility of the 27-spin weak-coupling momentum `k0=(1/3,2/3)` in 3-site-cell reciprocal coordinates, i.e. `k0·T1` and `k0·T2` are integers;
- shortest nonzero periodic lattice distance, aspect ratio, and geometry-induced finite-size artefacts;
- all numerically relevant momenta for rim Néel, hub Γ/finite-q and eventual triangular K reference. Distinguish physical primitive triangular momentum from 3-site-cell momentum folding.

Candidates, **subject to automated verification**:

- `N=18`, `Ncell=6`: `T1=(3,0)`, `T2=(2,2)`, determinant 6, K-compatible. This is geometrically skewed and can have unusually short winding vectors; search alternative integer supercells and select the most isotropic momentum-compatible shape(s), or explicitly document if all are strongly skewed.
- `N=27`, `Ncell=9`: reproduce original verified 3x3 geometry.
- `N=36`, `Ncell=12`: `T1=(2,2)`, `T2=(-2,4)`, determinant 12, K-compatible, equal lengths and 60-degree angle in the triangular Bravais metric; verify exact geometry and momentum phase factors before computation.

Do not assume K-compatible clusters at N=18 and 27 have comparable shape effects. If practical, run a second 18-spin supercell shape as a finite-geometry control. N=24/30 may be used as auxiliary non-K-compatible geometry controls, **not** silently combined into a K-sector extrapolation.

## 3. PHASE A (mandatory, feasible): 18-spin full ED sweep

- Implement or adapt the exact same matrix-free Heisenberg action, momentum sectors and total-spin measurement for the chosen 18-spin supercell(s).
- Run a coupling sweep `r=0.10, 0.15, 0.20, 0.225, 0.25, 0.275, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.675, 0.70`; refine around all level crossings, using energy intersections and spin/momentum sectors. Resolve crossing brackets to width <=0.005 wherever resources allow.
- **Verify the actual ground state** by sufficient sector coverage: at least all `Sz >= 0` sectors needed to certify candidate spin sectors or a mathematically justified SU(2) reduction plus momentum-space coverage. Do not label the lowest energy of one fixed `Sz`/momentum sector as the global ground state.
- For near-degenerate spin multiplets, report lowest energies for competing `S_tot` and cell momentum sectors, energy gaps and Lanczos residuals. Verify spin labels by `<S_tot^2>` and preferably `||[S_tot^2-S(S+1)]psi||` or ladder operators.
- Export ground-state checkpoints and static spin–spin correlation matrices at least at `r=0.2, 0.3, 0.4, 0.5, 0.6, 0.68`, plus every crossing endpoint.

## 4. Observables and normalization (identical for all sizes)

Define `S_H=sum_{hub i} T_i`, `S_R=sum_{rim i} S_i`, `N_H=Ncell`, `N_R=2*Ncell`, `N=3*Ncell`. Report:

1. `S_tot` and `m_tot = S_tot/N` (also give `sqrt(S_tot(S_tot+1))/N` optionally; do NOT mix the two conventions).
2. `H2=<S_H^2>`, `R2=<S_R^2>`, `HR=<S_H·S_R>`; verify `H2+R2+2*HR = S_tot(S_tot+1)`.
3. Hub pair polarization `p_H = (H2 - (3/4)*N_H) / (N_H*(N_H-1)/4)`; this is 1 for maximal ferromagnetic hub total spin and 0 for uncorrelated hub spins. Values may be negative.
4. `m_R,uniform^2 = R2/N_R^2`, and rim Néel correlation `m_N^2 = <(sum_R (S_A,R-S_B,R))^2>/N_R^2`, with clearly documented Fourier/local normalization. If the legacy ED quantity has a different normalization, provide both and conversion formula.
5. `S_h(q)` at Γ and relevant finite q, all momentum sectors and peak momenta (specify whether vector/spin-summed or z-only).
6. Optional: normalized intra-rim and hub–rim bond correlations, order-parameter distribution/structure factor in full set of allowed q, and low-energy gap for competing branches.

**Key caution:** the exact `S_tot` quantization changes between even N (18, 36: integer) and odd N (27: half-integer). Compare intensive quantities and scaling trends, NOT one-to-one spin plateaus or crossing values. S_H and S_R are generally not conserved separately. `H2` and `R2` do not prove pure subsystem spin eigenstates.

## 5. PHASE B: cross-size diagnostics

- Recompute or ingest verified 27-spin results under precisely the same definitions, local geometry and Fourier conventions. Reproduce their known energy and correlation benchmarks before cross-size comparison.
- Plot vs r, 18 and 27 separately with distinct markers: `p_H`, `m_N^2`, `m_R,uniform^2`, `m_tot`; annotate each size's level crossings and momentum change.
- Quantify whether the strong rise in hub polarization, subsequent increase of rim uniform total-spin correlations, and decline of `m_tot` recur in both sizes; check alternative 18-site shapes.
- If the 18-spin result disagrees strongly, analyze cluster shape, allowed q, and even/odd spin quantization before attempting the 36-spin calculation.
- Do not assign thermodynamic critical couplings from 18/27 data alone. Clearly separate exact observations from rotor/effective-field interpretations.

## 6. PHASE C (conditional): 36-spin HPC feasibility + selected points

First perform a **real** feasibility pilot at `r=0.30` or `r=0.45`: generate symmetry-reduced basis, estimate sector dimensions, RAM for vectors and matvec metadata, wall time per sparse matrix-free `H|v>`, and projected Lanczos iterations. Avoid full Lanczos basis storage; use restart/targeted symmetry sectors when valid. Check available memory before allocating. Note: `C(36,18)=9,075,135,300`; even after dividing approximately by 12 translations, individual momentum sectors are enormous, especially with complex vectors. Do not claim 36-spin ED has been completed without converged results and residuals.

If resources permit, calculate only `r=0.20, 0.30, 0.45, 0.60, 0.68`, with more refined values only when clearly useful. Coverage should certify the ground-state total spin and momentum. Report any restricted-sector calculation explicitly as such. If it is not feasible in the current environment, stop **only Phase C** after documented quantitative performance estimates; still complete and deliver Phases A/B, and propose the smallest credible HPC or tensor-network alternative.

## 7. Independent checks

- Pair-bond/coordination tests, full and reduced-sector Hamiltonian comparisons on small clusters, SU(2) rotational invariance and exactly `S(S+1)` total-spin eigenvalues, Hermiticity, and checkpoint reconstruction.
- Validate `H2+R2+2HR` sum rule to numerical tolerance.
- Cross-check selected static correlators from operator application and from the real-space correlation matrix.
- Run `r=0` and `r=1` sanity checks where meaningful and affordable (ensure `r=1` reproduces the intended triangular geometry).
- Record absolute eigenpair residual norms and tolerances; level-crossing brackets must exceed reported energy uncertainty.

## 8. Outputs

Produce `README.md`, `PHYSICS_REPORT.md`, `RUN_LOG.md`, `geometry_manifest.csv`, `gs_18_scan.csv`, `crossing_brackets_18.csv`, `size_comparison.csv`, `results_36.csv` (if available), eigenvector residual and completeness tables, source code, scripts, tests and checkpoints. Generate publication-quality PDF/SVG/PNG plots for ground-state spin, hub polarization, rim Néel and uniform components, and finite-size geometry (include 18/27/36 if available).

Provide an executive summary clearly distinguishing:
- Verified 18-spin results;
- Comparison with verified 27-spin results;
- Confirmed 36-spin results **or** quantified feasibility limitations;
- Inferences and hypotheses about thermodynamic physics.

Save reproducible outputs into a ZIP package and link it. Do not stop after reviewing these instructions. Start with the mandatory 18-spin calculation and proceed through the milestones in order.
