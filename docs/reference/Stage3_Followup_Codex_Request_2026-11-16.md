# Codex task: Complete the weak-coupling 27-spin dynamical spectral matrix

**Objective.** Continue, rather than repeat, the verified Stage 3 calculation in `stage3_verified_dynamics.zip`. Resolve the remaining low-energy hub-accessible spectral weight and its rim–hub cross residues near `J'/J = 0.2375`, then perform a controlled comparison with `J'/J = 0.2`. Avoid further phase-boundary scans unless needed to validate a state or sector.

## Model, source and target sectors

- Spin-1/2 stuffed-honeycomb Heisenberg model, 27 spins on nine three-site cells, `J=1`, periodic boundary conditions. Reuse the geometry, Hamiltonian kernel, symmetry routines, Fourier phases, and existing verified checkpoints. Follow `stage3_verified/dynamics/CONVENTIONS.md` exactly.
- At `J'/J=0.2375`: ground state `E0=-10.224890752094577`, cell momentum `k0=(1/3,2/3)`, `S_tot=1/2`, `Sz_tot=1/2`. Probe `k=(0,1/3)` takes the forward longitudinal target to cell `k_target=(1/3,0)`, still `Sz=1/2`. The reverse probe targets `(1/3,1/3)`.
- At `J'/J=0.2`: ground state `E0=-10.1859491284774`, same cell source momentum and `Sz` sector, same probe and target conventions. Reuse `hub_static_stages12_verified.zip` and the direct cross-moment reference `stage3_cross_moments_direct_27spin.csv` for validation.
- Define `O_n^z(k) = sum_R exp[-2 pi i(k1 x+k2 y)](S_A,R^z-S_B,R^z)/3`, `O_m^z` with `+` instead of `-`, and `O_h^z` with `T_R^z`; `sqrt(Nc)=3`, `Nc=9`. For target states `|f>`, report `v_a(f)=<f|O_a^z(k)|0>` and the phase-invariant residue matrix `R_ab(f)=conj(v_a(f))*v_b(f)`. Keep **complex** off-diagonal elements.

## Task A: Complete accessible low-energy spectral weight at 0.2375 (highest priority)

1. Start from the **two independently verified eigenpairs**, not the unverified 72-dimensional block-Krylov quadrature nodes. Verified excitation energies `omega/J=0.0002045122243590214` and `0.000235370036847371`; both have `S_f=1/2`. Their `h` weights are `0.009547913747585264` and `0.08181721183444425`. These are verified states, but are **not proven to be the two lowest levels of the entire target sector**.
2. In the cell target `(1/3,0)`, find **additional distinct, hub-accessible eigenstates** below `0.1J` relative to `E0`, prioritizing high hub spectral weight. Use symmetry-resolved, seeded, deflated/restarted Lanczos or a better justified matrix-free method. Re-orthogonalize explicitly against previously verified eigenvectors. Include the spin sectors allowed by a single `z` spin operator, `S_f=1/2` and `3/2`, and check all symmetry irreps/characters that can carry probe weight. Probe-seeded Krylov spaces may be used, but check for missed low-energy weight. Avoid dense factorization of the 2,228,700-dimensional target sector.
3. For each **converged** pole, export `omega`, absolute energy, cell momentum, `S_f`, rotation/reflection character when applicable, direct original-H residual, `S^2` residual, `n/m/h` diagonal weights, and complex `3x3` residue matrix. Prefer `||H psi-E psi|| <= 1e-8 J` and a comparably strong `S^2` check; report any looser acceptance transparently. Check that low-lying states are mutually orthogonal and distinguish genuine degeneracy from a repeated Ritz root.
4. Track cumulative hub weight `W_h(<Omega)` for thresholds `Omega/J = 0.01, 0.04, 0.05, 0.1`, with the **verified portion** separate from rigorously inferred spectral-weight bounds. Continue as resources permit, prioritizing the `0.1J` window. It is acceptable to report partial completion with certified residual spectral bounds; never label quadrature nodes or unconverged Ritz values as eigenpairs.

### Exact 0.2375 benchmarks (longitudinal `zz`, forward `k=(0,1/3)`)

- `M0_hh=0.33103928776800506`, `M1_hh/J=0.010068702554030905`, `M2_hh/J^2=0.01099926194038919`.
- Already-verified two-pole hub weight `W_known=0.09136512558202951` (27.60% of `M0_hh`).
- Excluding those poles, `M0_hh_remaining=0.23967416218597554`, `M1_hh_remaining/J=0.010047492568788176`, and remaining spectral centroid `0.04192146736698221J`.
- By positivity, **at least `0.1391992364980938` further hub weight must lie below `0.1J`**, i.e. at least 42.05% of the total hub `M0`. This bound does not determine how many poles carry that weight. It is a useful independent stopping/consistency diagnostic, not a presumed numerical result.
- Additional `M0/M1/M2` matrix references: `stage3_verified/dynamics/weak_moments.json` and `weak_exact_moments.csv`; all full `3x3` moments must be reconstructed by verified poles plus an explicitly labeled remainder, not by assigning fictitious energies to remaining weight.

## Task B: Matched coupling comparison at 0.2

1. Using **identical** `n/m/h` operators, forward/reverse cell momenta, `Sz`, spectral normalization, and residual tolerances, obtain the lowest resolvable hub-accessible poles at `J'/J=0.2` and their complex cross residues. Use the verified `0.2` checkpoint as source; revalidate its original-H residual and total spin.
2. Check the exact hub benchmarks `M0_hh=0.33177055287103`, `M1_hh/J=0.0069218759870593`, `M2_hh/J^2=0.0075353314103994`; reproduce the full `4x4` `n,m,h,Phi` moments in `stage3_cross_moments_direct_27spin.csv` where feasible.
3. Compare at fixed cell probe `k=(0,1/3)` the excitation energies (in both `J` and `J'^2/J` units), residues `R_hh`, `R_nn`, `R_mm`, and complex `R_nh`, normalized relative `R_nn/R_hh`, and any continuous spectral-weight transfer. Do not compare poles simply by sorted index: where possible match their spin/symmetry quantum numbers and overlap/adiabatic continuity. Make explicit which states remain unresolved.
4. Note as a reference that the **hub `M0` changes by only about -0.22%** between 0.2 and 0.2375, whereas the exact hub `M1` increases about 45.46%; these are **integrated moments**, not measured pole shifts.

## Cross-checks and deliverables

- Independently compute `M0=B^dagger B`, `M1=B^dagger(H-E0)B`, `M2=[(H-E0)B]^dagger[(H-E0)B]`; compare with full-precision static/cross sum rules. Enforce Hermiticity and positive-semidefiniteness of each physical pole's residue matrix and of the residual moment matrices when applicable. Check `k` versus `-k` using the correct partner source state, and explicitly exclude elastic ground-manifold transitions (the stated generic target has none).
- Retain numerical precision in CSV/JSON (no rounding of complex residues prior to sum checks), a detailed residual log and reproducible command line/code, and one concise physics report. Include a clear table of **verified poles vs approximate spectral quadrature**, and cumulative weight plots without confusing visualization broadening with intrinsic damping.
- Scientific interpretation: determine whether the previously unresolved low-energy hub spectral weight has sizeable **rim cross residues**, or is almost entirely hub dominated. Compare 0.2 with 0.2375. Do not infer a thermodynamic Goldstone velocity, magnon linewidth, Néel phase transition, or self-energy imaginary part from a finite 3x3 torus.
- If compute limits prevent full pole enumeration, stop with the strongest verified subset, certified remainder moments/bounds, and a targeted proposal for the next run. **Do not claim a complete spectrum unless completeness is independently demonstrated.**

Please execute the calculation (not merely review this task), and return the results package and physics report.
