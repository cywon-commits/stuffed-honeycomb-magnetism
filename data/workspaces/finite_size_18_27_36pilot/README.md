# Finite-size hub polarization and rim compensation

This is a completed 18-spin ED sweep and comparison with verified 27-spin data, plus a **restricted high-Sz 36-spin feasibility pilot**. No low-spin or global 36-spin ground-state energy is claimed. See `PHYSICS_REPORT.md` and `feasibility36.json`.

## Reproduction

Requires Python, NumPy, SciPy, Matplotlib, and g++ with OpenMP. The measured allocation was two CPU cores and 8 GiB RAM.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python run_all.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python run_all.py --tests
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=2 python run_all.py --compute --tests
```

The default regenerates the CSV/JSON, validates every 18-spin ground checkpoint by fresh H/S² actions and full-basis reconstruction, and makes PNG/SVG/PDF figures. `--tests` additionally runs independent full 12-spin Hamiltonian comparisons and the independent 36-site/two-down-spin operator test. C++ kernels compile automatically. If moving to another CPU, remove the `.so` files to force compilation with that CPU's `-march=native` flags.

`--compute` resumes the full 18-spin sweeps, crossing refinement, independent full-space checks, 27-spin checkpoint audits, and the 36-spin high-magnetization pilot. Existing converged checkpoints are reused. For a genuinely cold rerun, use an additional copy of the package with its `checkpoints` directory removed. Do not delete the verified `reference27` inputs. Runtime and roundoff may differ; residual criteria and grouped quantum numbers are the comparison standard.

## Geometry

Site IDs are `3*cell_index+s`, s=0(A),1(B),2(H). Coordinates are the original triangular primitive coordinates `(x+2y+s_is_B, x-y+s_is_H)`. Cell vectors are (1,1) and (2,-1); their metric is three times the 60-degree triangular Bravais metric. Integer supercell vectors in `geometry_manifest.csv` are in this **three-site-cell basis**.

The 18-spin K18 torus uses T1=(3,0), T2=(2,2). It is equivalent to an orthogonal rectangle with site-spacing lengths 3 and 3sqrt(3). Exhaustive index-six HNF enumeration finds three valid weak-k0-compatible lattices; all are point-group equivalent. Thus no independent, more isotropic K-compatible 18-spin torus exists. C18 is a separate non-k0-compatible geometry control, T1=(1,2), T2=(0,6), reduced lengths sqrt(12),sqrt(21), angle 79.1066 degrees. It must not be combined silently into a K-sector extrapolation.

Weak-coupling cell k0=(1/3,2/3) and physical triangular K=(1/3,2/3) use different reciprocal bases. Physical Q maps to cell k=(Q1+Q2,2Q1-Q2) modulo integers; physical K folds to cell Gamma on every ABC-cell torus. `momenta` lists every allowed cell q satisfying T_i·q integer. Γ is (0,0). Labels M are not assigned to arbitrary sectors.

Each valid graph has 3*Ncell unique AB edges and 6*Ncell unique hub–rim edges; hub coordination is six and rim coordination is three hub plus three rim neighbors. All translation permutations are explicitly checked as weighted bond automorphisms. At r=1 an independent primitive-coordinate nearest-neighbor construction checks all 18/27/36 graphs against the intended triangular torus.

## ED and completeness

The new 18-spin C++ action is matrix-free and uses bit-mask orbit representatives and actual periodic translations. It stores lookup metadata, not a Hamiltonian matrix. N_up=9, Sz=0 gives 48,620 unreduced states and all six cell-momentum blocks (dimensions sum exactly to 48,620). Every integer-S SU(2) multiplet has M=0, so this decomposition includes every total-spin ground-state candidate. No inference from one fixed momentum alone is used. Each sector has an independent cold random start and five lowest eigenpairs from restarted ARPACK, followed by direct H and S² residuals.

Near crossings, numerical H eigenvectors can have small admixtures of neighboring spin sectors. Joint H/S² resolution in near-degenerate subspaces, followed by H diagonalization within each S, removes this ambiguity. The three repaired diagnostic points and negligible energy changes are in `near_crossing_spin_resolution.json`. Other returned states were checked directly and already satisfied the final spin tolerance. Global checks at r=.4 and 1 also solve the full 48,620-dimensional M=0 space without translation reduction.

The specified grid is augmented by .0125 spacing and refinement of all detected quantum-number changes. Crossing endpoints, both competing energies, residual envelopes, and measured secant diagnostic points are recorded. The authoritative objects are the measured sign-changing intervals, not interpolated critical couplings. No claim excludes arbitrarily narrow re-entrant states between two unchanged sampled endpoints.

## Observables

All primary correlators are vector/spin summed: G_ij=<S_i·S_j>, G_ii=3/4. The selected fixed-M state suffices for these SU(2) scalars. Primary definitions are identical at N=18 and 27:

- H2=sum_hh G, R2=sum_rr G, HR=sum_hr G; H2+R2+2HR=S(S+1).
- m_tot=S/N; sqrt[S(S+1)]/N is saved separately. Neither is the chosen state's <Sz>/N.
- p_H=[H2-(3/4)N_H]/[N_H(N_H-1)/4]. Maximal hub spin gives 1; uncorrelated hubs give 0; negative values are valid.
- m_R,uniform²=R2/N_R².
- m_N²=sum_ij n_i n_j G_ij/N_R², n_A=+1,n_B=-1,n_H=0.
- S_h(q)=sum_RR' exp[-2πi q·(R'-R)] G_HR,HR'/N_H. All allowed q and maxima are exported. `S_h_zz_SU2_average=S_h/3` is a spin-multiplet average, not the selected M=0 state's zz value for S>0.

For legacy C_ab(q)=sum_RR' exp[-2πi q·(R'-R)]G_aR,bR'/Ncell, F_N(Γ)=(1,-1,0) C(Γ) (1,-1,0)^T=4Ncell*m_N². The normalized (A-B)/sqrt(2) mode is 2Ncell*m_N². Both conversions are included.

H2 and R2 do not assign conserved subsystem spins. Nonzero subsystem-S² variances are provided for representative 18-spin states and directly audited 27-spin checkpoints. No rotor approximation is used in the ED.

## Data and scope

- `gs_18_scan.csv`: all 18-spin global sampled minima, gaps and normalized observables (r=1 sanity rows are separate from the r<=.7 study).
- `sector_completeness_18.csv`, `all_low_levels_18.csv`: exhaustive momentum minima and the returned low eigenpairs, integer-S labels, direct H/S² residuals.
- `crossing_brackets_18.csv`: measured crossing intervals, branch energy differences, endpoint residuals, original brackets, and secant estimates.
- `checkpoints/K18`, `checkpoints/C18`: all sampled reduced eigenvectors, r/q/T/energy/spin and complete real-space G matrices. The full vector is reconstructed by `Basis.expand`; no expanded vector needs to be stored.
- `size_comparison.csv`, `hub_structure_factors.csv`: common normalization and all allowed q for 18/27.
- `reference27`: immutable checkpoint-derived input matrices, metadata and five actual audited ground checkpoints, with source SHA256. Original 27-spin indexing is explicitly mapped to the new cell order.
- `reference27_checkpoint_audit.json`: fresh H/S² and subsystem operator checks, not prose-only reference values.
- `results_36.csv`: actual **high-Sz pilots only**; `feasibility36.json` contains exact low-spin dimensions, measured pilot performance and separately labeled extrapolations.
- `tests.json`, `sanity_checks.json`, `cached_reproduction_validation.json`, checkpoint reconstruction CSVs and `validation.json`: independent checks and numerical errors.
- `computational_costs.csv`, phase records and `RUN_LOG.md`: measured costs and their scopes. Point timings exclude final checkpoint compression. H-call counts are target-H actions and do not count every auxiliary S² call. Summed walls overlap because the two shapes ran concurrently; they are not single elapsed turnaround. No low-spin 36-spin benchmark is represented as measured.

`python pilot36.py` regenerates the real N=36 high-magnetization pilot safely. It checks memory before allocation. Exact Sz=0 momentum dimensions use Burnside character traces and do not require generating the enormous low-spin basis. One complex translation-sector vector exceeds the 8-GiB allocation; Γ can be real but one real vector is already 5.63 GiB, and the present combinatorial metadata is 42.26 GiB. Point-group/HPC alternatives are design estimates, not completed ground-state calculations.

The archived measurements support finite-size hub polarization followed by rim compensation, with substantial shape/parity dependence. They do not establish thermodynamic long-range order, a bulk critical coupling, or a universal spin-plateau sequence.

The manifest separately tests physical triangular M admissibility at r=1: both chosen 18-site tori contain only (1/2,1/2), the original 27-site torus contains none, and the 36-site candidate contains all three. Away from isotropy a folded cell sector is not automatically a primitive-M eigenstate. Sparse 27-spin curves are drawn as markers only; connecting lines for dense 18-spin data are visual guides, not interpolated results.
