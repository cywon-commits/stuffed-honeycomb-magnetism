# Stage 3 follow-up review — 27-spin weak-coupling hub–rim dynamics

Date: 2026-11-16. Source: user-provided `stage3_followup_results.zip`, `stage3_followup_inputs.zip`, `stage3_followup_new_j02375.zip`, `stage3_followup_new_j020.zip`. Original result report `stage3_followup/followup/report.md`. Review independently sums CSV residues and compares against exact moment records; it did not repeat the expensive 27-spin eigenvector solves.

## Model and normalization

27 spin-1/2 sites, 3-site cell, 3×3 PBC; J=1. Source Sz_tot=1/2 and cell momentum k0=(1/3,2/3); longitudinal transfer k=(0,1/3) to target cell (1/3,0). Spin-operator Fourier factors normalized by 1/3 for 9 unit cells. Channels n=A-B, m=A+B, h=hub.

## Verified 27-spin poles

At J'/J=0.2: six converged poles (two previously known, four new) with omega/J approx 0.0000981875, 0.0001085245, 0.0010747878, 0.0146971670, 0.0156755074, 0.0157083647. Distinct total spins S_f=1/2 or 3/2. Hub spectral weight summed over these = 0.3292012341.

At J'/J=0.2375: seven converged poles (two previously known, five new) at omega/J approx 0.0002045122, 0.0002353700, 0.0018313718, 0.0210402297, 0.0226815266, 0.0227302241 and 0.0410993462. The final pole has tiny hub intensity 1.3552e-6. Hub weight below 0.04J (first six poles) = 0.3272222631.

## Exact moments and resolved fractions (from direct CSV sums)

| Quantity | J'/J=0.2 | J'/J=0.2375 |
|---|---:|---:|
| exact hub M0 | 0.331770552871 | 0.331039287768 |
| exact hub M1/J | 0.006921875987 | 0.010068702554 |
| exact hub M2/J² | 0.007535331410 | 0.010999261940 |
| verified hub weight under 0.04J | 0.329201234121 (99.2256% M0) | 0.327222263066 (98.8470% M0) |
| all verified hub M0 fraction | 99.2256% | 98.8474% |
| all verified hub M1 fraction | 37.8242% | 37.4677% |
| all verified hub M2 fraction | 0.5206% | 0.7371% |
| all verified rim Neel M0 fraction | 0.3688% | 0.5414% |
| magnitude of summed low-pole cross M0(n,h) versus exact | 98.3732% | 98.0264% |
| unresolved hub M0 | 0.002569318750 | 0.003815669493 |
| unresolved hub centroid/J | 1.67505 | 1.65009 |

The remaining small high-energy spectral fraction accounts for ~62% of hub M1 and over 99% of M2. The high-energy centroid is an intensity-weighted mean, not a single pole. The complete spectrum is not enumerated.

At 0.2375, the original two poles had combined hub M0=0.0913651256 (27.60% of hub M0). Additional verified poles contribute 0.2358584927, closing the earlier bound requiring at least 0.1391992365 below 0.1J. The reported moment/positivity bound constrains the missing hub weight below 0.1J to at most approximately 0.00020812; the 98.8474% directly verified under the cutoff is at least 99.9364% of the actual weight within the cutoff.

## Coupling dependence

Six bright branches match between the couplings by spin quantum numbers, little-group reflection, and physical-vector overlaps >0.9986. For bright low-energy poles the intensity ratio n/h is ~0.0065 at J'/J=0.2 versus ~0.0094–0.0098 at 0.2375, a relative increase of about 45–50%, even though n remains at the ~1% intensity level. The very lowest excitation energy increases by more than simple J'^2 scaling; several modes near 0.015–0.023J approximately follow a J'^2 trend. No single universal scaling is established from two points.

The resolved cross-spectral M0(n,h) sums to:

- J'/J=0.2: -0.0230106264158 - 0.0132851913555i, vs exact -0.0233911435604 - 0.0135048830313i.
- J'/J=0.2375: -0.0278314252320 - 0.0160684808496i, vs exact -0.0283917697126 - 0.0163919958863i.

Thus the hub-weighted cross overlap is mostly supplied by resolved low-energy poles, although these poles carry <1% of the total rim-Neel spectral weight. These are *coupled response* cross residues, not an independently extracted irreducible hub susceptibility or a Goldstone self-energy.

## Numerical verification

All four supplied ZIP files passed archive-integrity tests. The source `final_validation.json` records status PASS with 13 states, full M0/M1/M2 independent checks near 1e-15 (matrix norms), max new original-H residual 1.31e-11, max S² residual 9.78e-10, independent seed amplitude error ~1.53e-15, reflection and orthogonality checks. Those residuals were checked by Codex and read from its records, *not independently recomputed here* from the large eigenvector archives. CSV sums, fractions, and channel comparisons above were independently recomputed in this review.

## Conclusions and next steps

1. Weak-coupling hub spectral weight at finite cell momentum is very low energy and concentrated in six bright eigenstates; both spin-1/2 and spin-3/2 target sectors contribute.
2. Rim-Neel weight at the same momentum remains mostly outside this narrow hub-dominated window. Nonzero hub–rim residues grow moderately with coupling; mode intensity ratio is small, about 0.6–1.0%.
3. At this one finite momentum and N=27, no intrinsic damping, Goldstone velocity, or bulk ordered phase can be inferred.
4. Next prioritize the reduced coupled response and symmetry-constrained extraction near the true rim Goldstone momentum, with finite-size and twist/momentum sampling. The separate near-J'=J triangle-resonance hypothesis requires triangle/dimer/chirality observables in the strong-coupling region and is not tested by these weak-coupling spectra.

See `Stage3_WeakCoupling_Dynamics_2026-11-16.png` for a journal-style summary plot.
