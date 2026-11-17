# Verified triangular-lattice ED: corrected physical K and Gamma

All production calculations use N=27, N_up=14, Sztot=1/2, J=J'=1 and H=sum_b J_b S_i dot S_j. The full fixed-spin basis has 20,058,300 states. No full or sparse reduced 27-site Hamiltonian was constructed. The original 3x3 torus has three sites per cell: it is NOT a nine-site single-sublattice 3x3 triangular torus.

## Correct momentum labels and conclusion

In standard triangular reciprocal primitive coordinates, physical K=(1/3,2/3), with partner K'=(2/3,1/3). Physical Gamma=(0,0). They were independently verified from primitive nearest-neighbor vectors and the structure factor, then freshly benchmarked with the same optimized Lanczos implementation. The label K is NOT applied to primitive (1/3,1/3), or to the old requested code-cell sector (1/3,1/3).

Physical K has E=-15.125973385092 J and direct residual 8.296e-09 J. Physical Gamma has E=-14.529841875985 J and residual 8.047e-10 J. K is lower by 0.596131509107 J. The global ground-state primitive momenta are K and K'; their energies are equal by verified exact symmetry. The ground energy per site is -0.560221236485 J. Only one partner was independently timed.

Both physical K points fold to (0,0) in the original code's three-site cell coordinates. Therefore the original code-cell Gamma block contains both the physical Gamma and physical K states. This distinction is essential when interpreting a cell-symmetry-only solve.

There are no exact M points on this torus. No allowed finite-size momentum is labeled M, and no M eigenvalues are reported. The unavailable targets (1/2,0), (0,1/2), (1/2,1/2) would require different periods or boundary conditions.

## Independent primitive-vector and structure-factor verification

Choose e1=(1,0) and e2=(1/2,sqrt(3)/2) in Cartesian coordinates. Write integer primitive coordinates (u,v), with squared length u*u+v*v+u*v. The graph embeds with cell vectors a1=e1+e2, a2=2*e1-e2 and sublattice offsets A=0, B=e1, H=e2. The torus periods are 3*a1=(3,3), 3*a2=(6,-3) in primitive coordinates. Their determinant has magnitude 27. Original site numbering 3*(3*x+y)+s is preserved.

Periodic distance enumeration independently reproduces all 81 distinct undirected edges: 27 rim bonds and 54 hub bonds. Every site has degree six, without duplicate edges or self bonds. At J'=J the six nearest-neighbor displacements are +/-e1, +/-e2, +/-(e1-e2). This gives the isotropic triangular antiferromagnet; when J' differs from J only the cell translations are retained.

The nearest-neighbor structure factor is

f(q)=2[cos(2*pi*q1)+cos(2*pi*q2)+cos(2*pi*(q1-q2))]
    =|1+exp(2*pi*i*q1)+exp(2*pi*i*q2)|^2-3.

Thus its global minimum is -3 at (1/3,2/3) and (2/3,1/3), modulo reciprocal lattice vectors. f(Gamma)=6, whereas f(1/3,1/3)=0. verify_momenta.py extracts the six nearest-neighbor displacements from the graph's periodic edges, evaluates all 27 allowed momenta, verifies this analytic formula, and finds exactly those two minima. momentum_verification.json records the checks. This geometric argument identifies the Brillouin-zone K points independently of the ED energies.

Primitive momentum q folds to code-cell momentum k=(q1+q2,2*q1-q2) modulo integers. The code-cell Gamma block contains primitive (0,0), (1/3,2/3), (2/3,1/3). The old code-cell target (1/3,1/3) contains primitive (2/9,1/9), (5/9,7/9), (8/9,4/9); these are now named cell11_fold1/2/3, without a high-symmetry K label.

The original cell group is Z3 x Z3: allowed code-cell sectors are k=(n1/3,n2/3), n1,n2=0,1,2. All cell orbits have length nine at N_up=14, so each cell block has 2,228,700 states. At J'=J the full primitive translation group is Z3 x Z9, with 27 elements. Since 14 is not divisible by three, all its state orbits have length 27, giving 742,900 states per primitive block. The physical K/Gamma CSV explicitly gives n1,n2 for primitive k=(n1/3,n2/3), their folded code-cell momentum, and the primitive block dimension. General allowed primitive momenta on this 27-site tilted torus can also have ninths; it is not a primitive 3x3 nine-site torus.

The original cell PBC require 3*k1 and 3*k2 to be integers. Every requested half-integer M target fails this condition. In primitive coordinates the torus translation group also has odd order, so an order-two M character is impossible.

## Completed measured results and complete momentum coverage

| Primitive momentum q | Folded code momentum k | E/J | direct residual/J | iterations | wall s | mean H@v s | peak RSS MiB |
|---|---|---:|---:|---:|---:|---:|---:|
| (0, 0) | (0, 0) | -14.529841875985 | 8.05e-10 | 100 | 80.07 | 0.532 | 1423.3 |
| (1/3, 2/3) | (0, 0) | -15.125973385092 | 8.30e-09 | 70 | 39.72 | 0.330 | 1085.9 |
| (2/9, 1/9) | (1/3, 1/3) | -14.442647528132 | 6.93e-09 | 100 | 63.80 | 0.363 | 1426.8 |
| (5/9, 7/9) | (1/3, 1/3) | -14.611700724316 | 2.60e-09 | 90 | 53.48 | 0.322 | 1312.4 |
| (8/9, 4/9) | (1/3, 1/3) | -14.177170425055 | 6.64e-09 | 130 | 96.60 | 0.334 | 1767.0 |
| (0, 1/3) | (1/3, 2/3) | -14.266683821872 | 5.09e-10 | 120 | 102.95 | 0.474 | 1654.0 |


The first two rows are fresh, consistently optimized physical Gamma and K benchmarks. Their primitive sectors are respectively k=(0/3,0/3) and k=(1/3,2/3). Each has dimension 742,900. The hardware was an Intel Xeon Platinum 8573C with a two-CPU quota, two OpenMP threads and one BLAS thread, under an 8 GiB RAM cap. Peak RSS above is measured process memory, not a nominal array estimate. All rows are CONVERGED with a normalized direct residual below 1e-8 J. The displayed digits are numerical estimates at the residual verification scale.

The finite torus's exact 60-degree rotation and reflection preserve the weighted graph at isotropy. symmetry.py verifies the site permutations and bonds. Their six momentum stars (sizes 1,6,6,6,6,2) exhaust all 27 primitive momenta, as recorded in momentum_classes.json using integer numerators over nine. The six table rows represent these stars, so the ground-state comparison covers every allowed momentum. Exact symmetry partners were not independently timed. A single-vector Lanczos run does not measure eigenvalue multiplicities inside a block.

physical_K_Gamma.csv contains the requested physical-sector comparison, primitive_results.csv contains all six classes, and JSON files contain full precision and convergence histories. convergence.png/svg show recurrence residuals; black points mark the final directly recomputed residuals. Checkpoint readback independently verified all saved vectors; checkpoint_verification.json records the results. The bundle includes K_primitive.npz for direct ground-vector verification without repeating Lanczos. Other vectors remain in the workspace and can be regenerated.

A residual proves proximity to the spectrum, not a rigorous interval certificate of the lowest root. The lowest-root calculation uses deterministic random starts, complete reorthogonalization, stable convergence histories, independent small-cluster validation and complete symmetry coverage. No unconverged Ritz values are called eigenvalues.

## Independent Hamiltonian validation

An AST-isolated extraction of the uploaded geometry exactly matches the new bond list, up to ordering. Small-cluster references use independently enumerated bit states and direct spin-exchange rules. At 12 sites all explicit U-dagger H U momentum blocks agree with the matrix-free operator, including complex sectors and nontrivial stabilizers. The isotropic reference bond list is independently constructed by periodic nearest-neighbor distances; an anisotropic J'=0.73 calculation additionally validates the cell-only symmetry path. Maximum block Frobenius error: 2.168e-14. Maximum combined complete-spectrum error: 4.885e-14.

A three-site periodic quotient is checked against the analytic spectrum [-2.25,-2.25,2.25] J, with its periodic bond multiplicities retained explicitly. The 27-site fully polarized state gives 81/4=20.25 J. A dense 12-site block independently checks the Lanczos ground eigenpair. Tests also check unconverged status and momentum incompatibility rejection. All tests pass; see tests.json and test.log. No full 27-site matrix is built in validation.

The projected basis uses coefficients exp(-2*pi*i*q.g)/sqrt(L_r) on distinct orbit states T_g r, hence translation eigenvalue exp(+2*pi*i*q.h) under T_h. A bond flip of representative i to T_g r_j contributes (J_b/2)*sqrt(L_i/L_j)*exp(-2*pi*i*q.g) to gather row i. The independent projectors verify the sign and normalization. Stabilizer character checks are implemented even though the 27-site N_up=14 orbits are all free.

## Implementation, convergence monitoring and benchmarks

The optimized matrix-free kernel is C++ with OpenMP through ctypes. It uses bitwise bond flips, cached diagonal entries, split combinadic rank tables and state-to-orbit/translation maps. Row-gather parallelism has no concurrent output writes. It never stores the exchange destination list or the Hamiltonian matrix. Numba was unavailable, so the new implementation has no Numba dependency. No speedup versus the uploaded Numba program is claimed or measured.

Hermitian Lanczos performs two full reorthogonalization passes per step, using BLAS conjugate-transpose operations without materializing a conjugated Krylov basis. Every ten iterations it records the lowest tridiagonal Ritz value and recurrence residual. A small estimated residual triggers reconstruction, normalization and a fresh H@v residual check; only a direct residual below 1e-8 J marks convergence. At the bounded maxiter limit, the best pair is saved with converged=false. The current production cap is 180 vectors. This is a bounded full-basis solver, not an implicitly restarted solver; larger caps require proportionally more memory.

All active table measurements used the same optimized solver. solve_seconds includes Lanczos, reorthogonalization, direct residual checks and final compressed checkpoint, excluding model setup and two preliminary benchmark/Hermiticity products. Mean H@v time averages solve products. setup_seconds and a separate one-shot cold-cache product time are retained in JSON. RSS is getrusage process high-water memory. Production solves ran sequentially; timing variability is visible and no controlled cross-sector performance ratio is inferred.

One complex vector uses 11,886,400 bytes (11.34 MiB). The int32 state-orbit and uint8 translation maps use 100,291,500 bytes (95.65 MiB). A 180-vector Krylov array has 2,139,552,000 bytes of nominal payload (1.99 GiB); RSS depends on how many rows are touched before convergence. None of the active calculations hit the memory or iteration limit. Complete environment versions are in environment.json. Compiler: g++ 12.2.0, flags -O3 -march=native -fopenmp -shared -fPIC. The machine-specific binary is regenerated from source.

The original uploaded source is retained unchanged for provenance. Its K=(1/3,1/3) labels refer only to its chosen cell character, not the physical primitive K point. Its real embedding is mathematically valid for a lowest energy but duplicates eigenvalues; those duplicates are not physical degeneracies. Its first product timing includes Numba compilation, and its ARPACK maxiter parameter is not a matrix-vector product count. The new implementation replaces these limitations with explicit momentum labels, separate timings and per-iteration/direct-residual monitoring.

## Frequency-dependent correlations: feasibility assessment only

Broadened dynamical structure factors appear feasible using continued-fraction Lanczos seeded with S_q^alpha|psi0>. Longitudinal correlations stay at N_up=14; transverse operators require N_up=13 or 15 and correctly projected seed vectors. Those spin-operator applications and spectral sum-rule checks remain additional work. No frequency-dependent data were calculated.

A three-vector recurrence needs about 34 MiB of vector payload plus maps and workspace. Hundreds of products have an estimated kernel cost of a few minutes at the measured throughput, before seed construction, reorthogonalization and I/O. Saving 500 complex Krylov vectors alone would need about 5.54 GiB; larger bases can exceed the 8 GiB cap. Selective reorthogonalization, two-pass methods or disk-backed vectors can support longer spectra. Required length and resolution must be converged against broadening eta and checked using spectral weights and zeroth/first moment sum rules. These costs are estimates, not completed correlation benchmarks.

## Reproduction

Install requirements.txt and a C++/OpenMP compiler. From the extracted directory:

```sh
OPENBLAS_NUM_THREADS=1 python test_ed.py
OPENBLAS_NUM_THREADS=1 python verify_momenta.py
OPENBLAS_NUM_THREADS=1 python symmetry.py
OPENBLAS_NUM_THREADS=1 python run_all.py
OPENBLAS_NUM_THREADS=1 python verify_checkpoint.py K_primitive.npz
python make_report.py
MPLCONFIGDIR=./.mpl python plot_convergence.py
```

run_all.py skips previously converged result files; move the six named JSON files to force new solves. benchmark_physical.py always reruns the physical Gamma and K benchmarks. Checkpoints contain normalized Ritz vectors but are not a full Lanczos restart facility. verify_all.py recomputes residuals for all available checkpoint vectors using a freshly constructed operator. make_report.py refuses to report completion unless every active result passes the direct residual threshold. The SHA256SUMS manifest records bundle integrity. Both original and corrected labels are explained here; no current table assigns M to an allowed sector.
