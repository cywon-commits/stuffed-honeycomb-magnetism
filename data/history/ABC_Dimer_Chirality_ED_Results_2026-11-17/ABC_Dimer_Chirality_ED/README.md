# 27-spin scalar dimer/chirality dynamical ED

This package contains actual matrix-free ED results, exact direct-H moments, ground-manifold projections, and fully reorthogonalized block-Krylov spectral estimates. It does not identify approximate quadrature nodes as converged excitation energies. Read `report.md` for coverage and unresolved weight.

## Run

Python 3.12, NumPy, SciPy, Matplotlib and g++ with OpenMP are required. The measured environment had two CPU cores and 8 GiB RAM. Install `requirements.txt` if needed. Run with BLAS threads restricted:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python run_all.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python run_all.py --tests
```

The default reconstructs tables and figures from saved numerical results and performs validation. `--tests` also runs independent full-basis small-cluster comparisons and supplied local/reference scripts. The physical 27-site/two-up-spin test includes direct operators, an independently assembled 351-state Hamiltonian, and spectral moment sums. The 12-site test checks full operator seeds and reduced density matrices.

To rerun the expensive calculations:

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python run_all.py --compute
```

This reruns the priority queue, optional coupling extension, direct local tests and pole grouping. It writes large seed/eigenvector checkpoints locally; allow several GiB disk space and sustained CPU time. Existing `source_j1.npz` is reused; remove it to repeat the documented refinement of the original supplied checkpoint. The original checkpoint remains in `inputs/`. C++ kernels compile automatically from source. To force compilation for a different CPU, remove the `.so` files before running. `-march=native` binaries are machine-specific; source is included.

## Geometry and Fourier conventions

There are 27 spins, 81 bonds (27 rim AB, 54 hub bonds), 54 CCW elementary faces, and nine ABC cells. Work is in N_up=14, Sz=+1/2; the unreduced dimension is 20,058,300. A primitive triangular basis with 60-degree vectors has cell vectors (1,1) and (2,-1). Thus physical Q maps to cell momentum k=(Q1+Q2,2Q1-Q2) modulo integers. Physical K=(1/3,2/3), K′=(2/3,1/3) both fold to cell Γ. Physical Γ is distinct. This torus has no exact triangular M points.

Cell-face operators use exp[-2πi(k1*x+k2*y)]/3. At isotropy three primitive folds and two up/down face fields are obtained using the unitary geometric coefficients in `geometry_weights`. Each primitive face field has normalization 1/sqrt(27). Its normalized uniform up/down chirality combination is sum over 54 CCW faces divided by sqrt(54). Away from isotropy all six cell-face channels are retained for each D, Ddag and χ. A physical-Q *probe form factor* can still be applied to these fields, but primitive momentum is not an eigenstate quantum number there. All target eigenmomenta are cell momenta for J′≠J.

The first vertex of each CCW triangle fixes the geometric circular-dimer gauge. Odd face index has D_CCW=ω² D_ABC† and χ_CCW=-χ_ABC. Ddag fields use the conjugate geometric gauge. No Fourier normalization is changed to force agreement with the reference eigenmode weights. Complex mode coefficients are saved with the residues.

## Operators, residues and moments

Spins are Pauli/2. D=B_AB+ω B_BC+ω² B_CA, χ=S_A·(S_B×S_C), P=1/2-(2/3)sum B. All are SU(2) scalars. From the pure S=1/2 source only S_final=1/2 contributes, with ΔSz=0. Other members of a spin multiplet do not add transitions in this fixed-Sz scalar calculation. Spatial symmetry is verified independently of spin.

For b_a=O_a|g>, the stored convention is R_ab=<b_a|f><f|b_b>=a_a* a_b, where a_a=<f|b_a>. Every residue matrix is complex Hermitian positive semidefinite. Indices are ordered D_up,D_down,Ddag_up,Ddag_down,χ_up,χ_down at isotropy, and D_mu0..5,Ddag_mu0..5,χ_mu0..5 in cell calculations. Operator gauge changes transform the matrix and are not a physical change in response.

M0=B†B, M1=B†(H-Eg)B, M2=[(H-Eg)B]†[(H-Eg)B] are calculated by actual H action. M0 carries operator-squared units, M1 an additional J and M2 an additional J². The supplied 27-spin numerical reference gives M0, including zero-energy subtraction; it does not provide independent numerical M1/M2 tables for the 27-spin system. New M1/M2 are checked by independent tiny-cluster spectral sums and block-Krylov moment reproduction.

The K/K′ orbital ground pair is verified at isotropy. In each target sector, all matching verified ground vectors are projected out of seeds. The zero-energy matrix is saved separately. Source spectra are averaged with weights 1/2 only after individual deflation. At Γ the conjugate source response follows exact conjugation; finite-Q sources are solved separately. The subleading Γ chirality zero weight is the diagonal selected-state expectation (orbital-ensemble variance); the physical-K dimer zero weight can be an off-diagonal K↔K′ transition. Away from isotropy the previously verified unique fixed-Sz cell ground vector is removed at cell Γ.

Local circular-dimer sum rules include zero frequency. The difference of positive-energy local integrated weights equals 3sqrt(3)<χ> minus the circular difference in ground-projected weights. This identity applies to a single oriented face, or a correctly reconstructed all-Q sum, rather than each Q separately. `fresh_local_sum_rules.json` records both total and deflated versions with new source RDMs.

## Solvers and convergence

The Hamiltonian-vector product is on-the-fly C++/OpenMP. No full 27-spin Hamiltonian is constructed. `eigsh` uses implicitly restarted ARPACK (the complex interface uses Arnoldi). Source-derived scalar seeds cover accessible channels; all returned states receive fresh H and S² residual checks. Higher-spin states are recorded as forbidden checks, not spectral poles. A rank-one +25J projector shifts already known ground states out of the eigensolver search without changing their orthogonal excited states. Other compatible little-group reductions preserve the Hamiltonian.

Converged physical poles require direct H residual <1e-8J and S² residual <1e-7. The source refinement and actual per-vector errors are saved. A low Ritz value alone never qualifies. Arnoldi iteration counts are recorded as iterative H actions; restart-cycle counts are not exposed by SciPy. Verification H/S² calls are counted separately where possible.

Exact degeneracies require an entire subspace residue. Complex ARPACK vectors can be nonorthogonal within a degenerate eigenspace: use V(V†V)^-1V†, not a naive sum of their individual outer products. `pole_subspace_data.npz` stores their Gram matrices and seed amplitudes. `degenerate_pole_groups.json` and canonical tables use basis-independent grouped residues. Conjugate spatial partner sectors are explicitly marked and receive fresh residual checks when vectors are available. Physical-Q star exports are exact symmetry-derived results, not additional eigensolves.

Block-Krylov calculations use two-pass full reorthogonalization and bounded bases, saving numerical quadrature nodes in separate files. Their M0/M1/M2 reproduction checks quadrature consistency but does not independently converge each node. All plots distinguish these estimates from verified eigenstate lines. No artificial broadening is used, and no physical linewidth or lifetime is inferred.

## Files and cost interpretation

- `poles.csv`, `ensemble_poles.csv`: verified scalar eigenspaces and orbital averages.
- `residues.npz`, `complex_residues.csv`, `responses.json`: full complex matrices and fixed reference-mode coefficients.
- `static_moments.csv`, `ground_projection.csv`: exact direct-H moments and separate zero-energy contributions.
- `mode_summary.csv`, `mode_poles.csv`, `remaining_weight_bounds.csv`: resolved/unresolved weights and positive-moment cutoff bounds.
- `symmetry_partner_poles.csv` and companion NPZ: exact sixfold-star partners and K/K′ comparison.
- `mode_convergence.csv`, per-solve progress/JSON/logs, `validation.json`: direct residuals, spin checks and validation coverage.
- `computational_costs.csv`, `queue_records.json`, `extension_records.json`, `environment.json`: measured costs. Early seed metadata measures setup/operators/moments before checkpoint compression; queue task walls include complete process/file-writing costs. Later seed records explicitly include writing. CPU-seconds are blank where they were not instrumented; CPU quota/model and process peak RSS are always distinguished from estimates. Summed task walls can overlap and are not total elapsed turnaround.

The compact archive retains original/refined ground checkpoints, input references, source, logs, exact moments, amplitudes/Gram matrices and all result tables. Large transient scalar seeds and excited-state eigenvector arrays stay in the computation workspace; `--compute` regenerates them. Cached reproduction uses the saved subspace data, so the compact archive does not need transient vectors. Per-file SHA256 records give source and input fingerprints.

The finite 27-site spectrum establishes finite-system scalar dynamics only. It does not establish thermodynamic coherent resonance, chiral or dimer order, a phase boundary, or magnon damping. The physical K magnetic static reference is a different tensor channel and cannot be equated to scalar spectral weights.

The optional J′/J=0.7 search targets the one-dimensional cell-Γ symmetry sectors needed for the two physical-Γ chirality face form factors. E-character sectors and finite cell momentum are not independently enumerated there; their exact moments may still be tabulated, with unresolved weight explicit.

Scalar projected weights that are negative solely by roundoff below 1e-12 are clipped to zero in the leading-mode summary and plots. Original complex matrices, their tiny negative eigenvalues and validation errors are preserved for audit. Operator-call counts refer to the named target model; auxiliary source-validation calls are not included in all early target counters.
