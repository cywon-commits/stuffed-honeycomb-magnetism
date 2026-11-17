# Static stuffed-honeycomb ED: completed Stages 1–2

Only static calculations are performed. No continued-fraction spectrum, frequency broadening, susceptibility or dynamical linewidth has been computed. The attached specification was reviewed first; SPEC_REVIEW.md records the Fourier, normalization, ensemble and selection-rule decisions.

## Model and state conventions

N=27 spin-1/2 sites, 3x3 three-site cells, periodic boundary conditions, J=1 and Jprime/J in the completed set [0.2, 0.5, 1.0]. Every ground-state search is at N_up=14, Sz=+1/2, full dimension 20,058,300. The production Hamiltonian is matrix free. At unequal couplings only the correct cell translations, plus independently validated compatible point operations, are used. Static measurements are evaluated after converting the saved eigenvector to the cell-translation basis of dimension 2,228,700.

The previously verified isotropic ground vector is reused, not redetermined: primitive K=(1/3,2/3), folded cell Gamma. Its primitive Kprime partner is symmetry degenerate. The original full-momentum verification is preserved in isotropic_ground_verification.json and isotropic_ground_sector_reference.csv. No exact M point exists here, and no allowed sector receives an M label.

Primary matrices are C_ab(k)=sum_spin <O_a(k)^dagger O_b(k)>, O_a(k)=sum_R exp[-2*pi*i*k.R] S_Ra /sqrt(9). We average the scalar real-space correlations over the twelve exact hub-centered point operations and over the ground spin multiplet. Scalar dot correlations are unchanged by spin-multiplet averaging, while vector means vanish. The primary raw scalar matrix is therefore also the connected scalar matrix of that specified ensemble. This removes arbitrary choice of the degenerate momentum-star partner. Separately, CSVs preserve the actual selected Sz=+1/2 state's raw and connected scalar and longitudinal zz matrices. They are not substituted for one another.

## Verified ground states and key results

| Jprime/J | E/J | cell momentum representative | S | rim m_N squared | hub S(Gamma) | 120-degree intensity | energy residual |
|---:|---:|---|---:|---:|---:|---:|---:|
| 0.2 | -10.185949128477 | (1/3,2/3) | 0.5 | 0.24532220 | 0.09517246 | 2.24468141 | 3.28e-14 |
| 0.5 | -11.083587020439 | (0,0) | 2.5 | 0.22257082 | 2.53725956 | 3.25999050 | 3.16e-14 |
| 1 | -15.125973385092 | (0,0) | 0.5 | 0.16383234 | 1.97524731 | 2.94898207 | 8.30e-09 |


From isotropy to Jprime/J=0.2, the rim Neel squared amplitude changes from 0.16383234 to 0.24532220; hub S(Gamma) changes from 1.97524731 to 0.09517246. At Jprime/J=0.5 the corresponding values are 0.22257082 and 2.53725956. These are equal-time finite-size correlations, with all definitions and state averages held fixed.

- Jprime/J=0.2: hub peak S=0.99531166 at cell momentum numerator(s) [[0, 1], [0, 2], [1, 0], [1, 1], [2, 0], [2, 2]] over 3; rim NN=-0.37033550, hub-rim bond=-0.01730469, all-distinct hub-pair mean=-0.08185344; uniform rim/hub Gamma cross=-0.01208820.
- Jprime/J=0.5: hub peak S=2.53725956 at cell momentum numerator(s) [[0, 0]] over 3; rim NN=-0.32302122, hub-rim bond=-0.08748200, all-distinct hub-pair mean=0.22340745; uniform rim/hub Gamma cross=-1.08303935.
- Jprime/J=1: hub peak S=1.97524731 at cell momentum numerator(s) [[0, 0]] over 3; rim NN=-0.18674041, hub-rim bond=-0.18674041, all-distinct hub-pair mean=0.15315591; uniform rim/hub Gamma cross=-1.94746953.


Only measured, directly residual-verified eigenvectors enter these tables. The fixed-Sz value is not used to guess S: total S-squared is applied independently and its eigen-residual recorded. Detailed convergence, little-group dimensions and all checked energies are in sector_convergence.csv and per-block JSON. Digits in CSV are numerical estimates, not assertions of an equal number of rigorous significant digits.

## Independent ground-sector determination at each unequal coupling

The weighted stuffed graph has three cell-momentum stars: Gamma; the six momenta (0,1),(0,2),(1,0),(1,1),(2,0),(2,2) over 3; and the pair (1,2),(2,1) over 3. Every class is checked independently at each coupling. Their symmetry operations preserve rim and hub bond weights, as verified explicitly.

Additional compatible cyclic little groups reduce cost without discarding sectors: C6 at Gamma (characters ell=0,...,5), one compatible reflection C2 at representative (0,1/3), and C3 at representative (1/3,2/3). At Gamma complex conjugation of the real Hamiltonian equates ell=1 with 5 and ell=2 with 4; ell=0,1,2,3 are independently solved. Both reflection parities and all three C3 characters are solved. The nine measured blocks consequently cover every cell momentum and character. coverage.py checks that each representative cell block sums to 2,228,700 and that symmetry-star multiplicities recover all 20,058,300 states. A positive residual alone does not rigorously certify the lowest root, but standard lowest-root Lanczos, complete symmetry coverage and independent small-cluster spectrum checks underpin the numerical determination.

Jprime/J=0.2: nine representative symmetry blocks converged; summed measured solve wall time 1258.34 s; largest per-process RSS 997.7 MiB. Ground file j0.2_n12_r0, cyclic order 3, character index 0.

Jprime/J=0.5: nine representative symmetry blocks converged; summed measured solve wall time 1086.58 s; largest per-process RSS 1415.0 MiB. Ground file j0.5_n00_r3, cyclic order 6, character index 3.


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

| Jprime/J | analysis wall s | pair kernel s | analysis peak RSS MiB | S-squared eigen-residual |
|---:|---:|---:|---:|---:|
| 0.2 | 61.66 | 18.99 | 466.9 | 5.80e-13 |
| 0.5 | 23.90 | 7.30 | 441.4 | 4.56e-13 |
| 1 | 49.06 | 11.93 | 452.9 | 1.72e-08 |


Ground searches use restarted real-Hermitian ARPACK Lanczos with a real embedding only when the projected operator is complex. Real characters use a real block directly. The embedding duplicates eigenvalues algebraically; duplicates are not called physical degeneracies. A fresh complex-vector H residual below 1e-8 J is mandatory. ARPACK tolerance is 1e-11; subspace size ncv and exact H-product counts are saved per run. SciPy does not expose its internal restart iteration count: that field is null, rather than fabricated from H-product counts. Progress records give H applications, wall time and RSS, not unverified intermediate eigenvalues.

The environment has a two-CPU quota and 8 GiB RAM cap; OpenMP uses two threads and BLAS one. Full versions and hardware are in environment.json. All wall times and RSS are measured in this shared environment; brief validation/analysis jobs can overlap a solve, so these are not controlled speed comparisons. New solve_seconds includes ARPACK and its fresh residual check, excluding model setup and compressed checkpoint I/O. analysis_seconds measures core conversion, operator checks and static pair work before final CSV/JSON output and plotting. Compilation/setup and solve times are separated. The initial translation-only weak-coupling trial was stopped in favor of validated additional symmetries; translation_only_trial.json explicitly marks it unconverged and supplies no eigenvalue. None of its intermediate values enter static results.

## Recommendations for a later dynamical stage — not executed

Use Jprime/J=1 as the validation anchor. The most informative contrast is 0.5 versus 0.2: the former has an S=5/2 multiplet and strong uniform hub correlations, while the latter has S=1/2 and dominant finite-momentum hub correlations. Prioritize Gamma n,m,h at 0.5 and the (0,1/3) transfer star at 0.2, with Gamma Ward controls at both. A few further static points between 0.2 and 0.5 could locate the finite-cluster sector crossing before a fine spectral scan; no such scan is performed here. Use cell Gamma to probe rim-staggered n, uniform m, and h, retaining the full m/h cross matrix; the inelastic conserved combination m+h must be null. Inversion forbids n/h at Gamma in the stated ensemble, so a zero there cannot settle hub coupling. Operator transfer k sends a chosen ground momentum k0 to k0+k, not simply to k. recommended_dynamics.csv gives proposed source, transfer and target cell sectors without running a spectrum. Also target the allowed finite cell momentum star where the hub static structure factor peaks, and compare the other inequivalent nonzero star as a control. The per-coupling peak sectors are listed above.

At isotropy use physical basis-phase-resolved K/Kprime probes separately from physical Gamma, even though they share cell Gamma. Explicitly project the complete degenerate ground manifold out of the strictly inelastic seed: one-point subtraction alone does not remove every zero-energy spin- or momentum-partner transition. Large equal-time weight prioritizes channels but does not prove low-frequency spectral weight or a Goldstone linewidth. Eta convergence, spectral moments and the inelastic Ward identity belong to the unstarted Stage 3.

## Reproduction and saved artifacts

Install requirements.txt and a C++ compiler with OpenMP. Run `OPENBLAS_NUM_THREADS=1 python run_static.py` from this directory. It reuses verified results/checkpoints when present; `--force` reruns all sector solves and statics. Individual commands are documented in README.md. The original isotropic checkpoint is included. Ground checkpoints at the two unequal couplings are included when available; non-ground checkpoints are optional and can be regenerated. Only small reference tests assemble Hamiltonian matrices. All production H applications and static pair reductions are matrix free.

The principal numeric files are static_summary.csv, sector_convergence.csv, correlation_eigenmodes.csv, per-coupling static_j*_C.csv and static_j*_channels.csv, real-space NPZs, JSON validation/convergence data, and ground_j*_coverage.json. Figures static_evolution, static_matrices, hub_modes and total_spin_and_polarization are supplied as PNG and SVG; lines between three points only guide the eye. SHA256SUMS records deliverable integrity.
