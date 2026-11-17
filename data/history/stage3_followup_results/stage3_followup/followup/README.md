# Stage 3 follow-up: 27-spin weak-coupling dynamical ED

Start with report.md, verified_poles.csv, complex_residues.csv, cumulative_weight_bounds.csv and final_validation.json. Every physical pole is independently checked with the original Hamiltonian and total-spin Casimir. The prior 72-node block quadrature is retained only as explicitly approximate reference data.

All poles use source M=1/2, cell k0=(1/3,2/3), probe (0,1/3), target (1/3,0). n=(A-B)/3, m=(A+B)/3, h=H/3 with negative cell Fourier exponent. Raw complex amplitudes depend on eigenvector gauge; residue matrices conj(a_i)*a_j do not. Both little-group reflection characters are searched. Generic cell momenta have no C6 rotation label. There are no triangular M points on this 3x3 torus.

## Cached reproduction

From the full archive root:

    OPENBLAS_NUM_THREADS=1 python followup/run_followup.py

Add --tests for the small-complex-sector deflation/bound tests. Cached reporting does not rerun expensive eigensolvers. The compact results archive contains all numerical tables, scripts and validation records; the full archive also contains the source and new eigenvector checkpoints. Cached reporting from the compact archive reuses the delivered verified overlap cache. The checkpoint archive is required for fresh vector-overlap recomputation, revalidation and expensive ED.

## Recompute

    OPENBLAS_NUM_THREADS=1 python followup/run_followup.py --tests --compute

Verified files and vectors are reused. Add --fresh to move the delivered new-eigenpair records/checkpoints into a backup directory and rerun the new solves. The original six input checkpoints remain intact. Requires Python, NumPy, SciPy, Matplotlib, g++/OpenMP; environment.json and requirements.txt record the session versions. Model construction compiles kernel.cpp on demand. No full 27-spin Hamiltonian is assembled.

The eigensolver is SciPy's implicitly restarted Arnoldi/Lanczos interface for Hermitian complex operators, with explicit physical-vector projector deflation, ncv=64, tol=2e-12, and direct residual acceptance. Each iterative callback has a 1200-second wall limit; verification/checkpoint I/O can extend final wall time. Some comparison and final continuation solves use a one-time exact Casimir polynomial to start in S=1/2 or S=3/2. Original-H, spin, reflection and orthogonality checks determine acceptance, not the seed label. All candidate states are rechecked using two-pass Gram-Schmidt and independently formed bit-diagonal seeds.

Source checkpoints live in sibling crossing/ and static/ directories. New eigenvectors are in momentum-plus-reflection reduced bases; recreate CyclicModel(cell_n=(1,0), angular=r, up=14, jp=j) before acting on them. Direct_inputs.py computes source-to-target Fourier seeds at representative spin configurations, projects the compatible reflection, and verifies the known residues. This reproduces the older full-cell/C++ operator seeds more economically; final_validation.json records the independent amplitude and full-moment checks.

## Bounds and limitations

Moment sums of listed poles are incomplete. spectral_remainders.json contains exact operator moments, verified sums, unresolved matrices, individual PSD and block-Hankel tests. Cumulative bounds use positivity, a global spectral upper bound, and one-sided variance inequalities. They are conditional on the verified numerical input and include a 1e-8 safety margin; they are not formal interval-arithmetic certificates. No unresolved mass is assigned a fictitious energy.

The direct four-channel reference CSV named in the specification was not attached. The supplied diagonal benchmarks are reproduced, and all n,m,h,Phi moments are independently checked using real-space Gzz and bond f-sums. Phi is linearly dependent on n,m, so its moment matrix is rank deficient. Correct conjugate-source validation is separate from reversing the probe around the same source.

One duplicate comparison worker was stopped before export. Its partial compute is disclosed in duplicate_worker_note.json; the corresponding early raw progress log can contain interleaved lines/null padding. Accepted-root JSON and final validation are authoritative for H calls, timing and eigenpairs. Setup, validation and overlapping jobs are not a dedicated-machine benchmark. No broadening, damping, Goldstone velocity or phase transition is inferred.
