# Targeted Stage 3 dynamical ED

All residues refer to longitudinal Sz=1/2 sources, with 1/sqrt(9) Fourier normalization. No broadening is applied. These finite-cluster peaks do not determine thermodynamic damping. Read CONVENTIONS.md for the reviewed translation, spin, Fourier, and sum-rule definitions.

## Second crossing: resolved visible rim Neel poles

| Jprime/J | omega/J | n-z residue | fraction of n M0 | largest state residual/J |
|---:|---:|---:|---:|---:|
| 0.3875000 | 0.005756116931 | 1.315149452287 | 41.156489% | 6.98e-11 |
| 0.3921875 | 0.003820478184 | 1.316313489511 | 42.203797% | 7.40e-11 |

The competing S=7/2, chi=0, F=+ and S=5/2, chi=3, F=- eigenvectors pass fresh original-H, total-spin, rotation and reflection checks. This establishes substantial low-energy dipolar weight, rather than symmetry allowance alone. The low poles carry only about half a percent of the n first moment; appreciable high-energy weight coexists with the soft low-energy peak. The forbidden m/h amplitudes are below 1.3e-16. The first crossing has DeltaS=3 and no direct single-spin transition.

Gamma n elastic weight is below 1.1e-30. Fixed-M m/h elastic weights are retained separately in the JSON. Projecting the same-M ground state suffices for longitudinal response: other members of its spin multiplet have different M and cannot contribute. The m+h seed equals v/6 within 6e-17; its inelastic part vanishes. Scalar full-multiplet elastic weights are not subtracted from zz spectra.

## Weak coupling j=0.2375

Source cell (1/3,2/3), probe (0,1/3), target (1/3,0), all at M=1/2. Neither ground corner partner lies in the target, so the elastic ground-manifold contribution is zero. Exact forward/reverse three-channel moment matrices and complex cross terms are in weak_exact_moments.csv.

The target-sector calculation produced 2 independently verified eigenpairs by spatial-star transfer of saved symmetry-character minima. It is a partial low-energy spectral matrix, not a complete spectral decomposition. Individual complex 3x3 residues are in weak_complex_residues.csv; original-H and S-squared residuals are in weak_verified_poles.csv. Target verification: 2 H calls, 80.03 measured solve seconds, 658.9 MiB process peak RSS.

| omega/J | S | n residue | m residue | h residue | H residual/J |
|---:|---:|---:|---:|---:|---:|
| 0.000204512224 | 0.5 | 9.65615e-05 | 3.25612e-05 | 0.00954791 | 3.76e-12 |
| 0.000235370037 | 0.5 | 0.00078656 | 0.000269871 | 0.0818172 | 6.43e-14 |

Exact moments minus the verified-pole sums are saved as moment_remainders. They quantify unsaturated spectral weight and are not treated as validation failure. Higher-energy poles and dipole-forbidden levels have not been exhaustively enumerated.

A separate 72-dimensional block-Krylov quadrature approximates the full energy-resolved matrix. Its nodes are explicitly **not verified eigenvalues** and must not be combined with the verified poles as extra states. It reproduces exact zeroth/first/second moments with matrix-norm errors [2.8445917799073156e-15, 7.716229546889964e-15, 2.444529058792141e-14]. Quadrature time 191.55 s, peak RSS 3297.5 MiB. No convergence of its individual nodes is claimed.

## Moment validation and computational limits

All 30 available full-precision CSV comparisons of M0zz and symmetrized M1zz agree within 6.193e-12. Forward M1(k) and transposed reverse M1(-k) differ by at most 1.977e-14. There is no supplied M2 reference; exact M2 is newly computed from DB norms. A separate 12-spin seed test checks Fourier phases and complete spectral M0/M1/M2 sums (errors below 1.1e-14).

The two-CPU quota and matrix-free H applications limit high-energy eigenpair enumeration. Exact moments, verified low poles, and moment-complete approximate quadrature are deliberately distinguished. Timings are per operation and include contention; they are not a dedicated-machine performance estimate. The discarded cold and probe-started weak searches have no exported Ritz values; solver_restart_note.json documents its excluded partial cost. No linewidth or thermodynamic transition is inferred. Optional first-crossing spectra are deferred.
