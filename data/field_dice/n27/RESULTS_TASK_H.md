# 27-spin Task H results: wedge scaling and the ferrimagnet threshold

Specification: `docs/guides/27SPIN_3차_의뢰서_H_KO.md`
Couplings J = (1, j + d/2, j − d/2). Raw data: `task_H.csv` (933 rows) · Derived values: rows H1 and H2 of `summary_n27.csv` · Unit test: U10 in `unit_tests.json`

## Verdict

### H1 — scaling of the first-jump threshold d_c: **neither linear nor quadratic; d_c ∝ j^1.38**

| j | **d_c interval (N = 27)** | d_c / j | d_c / j² | Linear prediction 0.056–0.058 j | Quadratic prediction 0.28–0.39 j² | Overlap (linear / quadratic) |
|---|---|---|---|---|---|---|
| 0.100 | **(0.004375, 0.0046875]** | 0.044–0.047 | 0.44–0.47 | 0.0056–0.0058 | 0.0028–0.0039 | no / no |
| 0.125 | **(0.0059375, 0.00625]** | 0.048–0.050 | 0.38–0.40 | 0.0070–0.0073 | 0.0044–0.0061 | no / yes (upper edge only) |
| 0.150 (Task C) | (0.0075, 0.00875] | 0.050–0.058 | 0.33–0.39 | 0.0084–0.0087 | 0.0063–0.0088 | – |
| 0.175 | **(0.009375, 0.0096875]** | 0.054–0.055 | 0.31–0.32 | 0.0098–0.0102 | 0.0086–0.0119 | no / yes |
| 0.200 (Task C) | (0.01125, 0.0125] | 0.056–0.063 | 0.28–0.31 | 0.0112–0.0116 | 0.0112–0.0156 | – |

- **Decision rule (section 2 of the specification).** "Both j = 0.1 and 0.125 overlap linear only" does not hold, and neither does the quadratic version. At j = 0.1 the interval overlaps neither prediction. So, as the rule requires, d_c = a·j^p was fitted to all j values.
- **Fit.** Least squares on the midpoints of the five intervals gives **p = 1.38, a = 0.108**. The range of p for which one power law passes through all five intervals at once is **1.26 ≤ p ≤ 1.42**. The linear (p = 1) and quadratic (p = 2) laws both fall outside this range.
- **Interpretation.** d_c/j rises steadily with j (0.045 → 0.059) and d_c/j² falls steadily (0.45 → 0.30). The 27-spin first-jump threshold therefore lies between the two pure scalings in this range of j. The earlier "linear" reading came from the two points j = 0.15 and 0.2 alone, which overestimated the coefficient (0.056–0.058). The "quadratic" reading of the manuscript draft does not fit the small-j points either.
- **S at the jump:** 1/2 → 3/2 for every j. The S = 1/2 ground state is at K, K′, and the S = 3/2 ground state is at Γ. S was confirmed with S² at both ends of every interval.

### H2 — ferrimagnet threshold d_full: **S jumps directly from 3/2 to 9/2, and d_full ∝ j (p = 1.01)**

| j | First jump (S: 1/2→3/2) | **d_full (S: 3/2 → 9/2)** | d_full / j | d_full / j² | Theory 0.34–0.37 j² | Classical j² |
|---|---|---|---|---|---|---|
| 0.10 | (0.004375, 0.0046875] | **(0.0159375, 0.01625]** | 0.159–0.163 | 1.59–1.63 | 0.0034–0.0037 | 0.0100 |
| 0.15 | (0.0075, 0.00875] | **(0.0240625, 0.024375]** | 0.160–0.163 | 1.07–1.08 | 0.0077–0.0083 | 0.0225 |
| 0.20 | (0.01125, 0.0125] | **(0.0321875, 0.0325]** | 0.161–0.163 | 0.80–0.81 | 0.0136–0.0148 | 0.0400 |

- **Staircase.** For all three j the ground spin is S = 1/2 → 3/2 → **9/2**. No S = 5/2 or 7/2 ground state appears; at the second jump S goes straight to the ferrimagnetic value N/6 = 9/2. The thresholds d(→5/2), d(→7/2) and d(→9/2) are therefore the same interval. S was confirmed with S² at both ends of every interval (3/2 / 9/2).
- **Fit.** **d_full = 0.164 · j^1.01**. The range of p consistent with all three intervals is 0.99–1.03, so d_full is linear in j. d_full/j is constant to within 0.159–0.163.
- **Comparison with theory.** The hub-saturation prediction 0.34–0.37 j² is 2.2–4.5 times smaller than the measured d_full and has a different scaling (j² against the measured j¹). The classical non-collinear wedge j² does not match either (d_full/j² = 1.6 → 0.8).
- **Comparison with 18 spins** (`pd_wedge_18.csv`, K18: S jumps 0 → 3 = N/6, directly to the ferrimagnet). δ_c(K18) = 0.0212, 0.0315, 0.0418 at j = 0.10, 0.15, 0.20, which is 0.209 j. The 27-spin d_full is 0.161 j. **At both sizes the threshold for reaching the ferrimagnet is linear in j**, and the coefficient falls from 0.21 to 0.16 (ratio 0.77) on going from 18 to 27. C18 (0 → 1) does not reach the ferrimagnet, so it cannot be compared here.

## Grid of H2 ground-state levels (lowest per S^z; S = 3/2 except where marked)

The distance from the S^z = 5/2 level to the S = 3/2 ground state shrinks steadily with d. Just before the jump it is below 1e-3, and at the jump the S^z = 3/2 … 9/2 levels all become degenerate together. The full table is in `task_H.csv`.

| j | d | E0 | E(5/2) − E0 |
|---|---|---|---|
| 0.10 | 0.0047 / 0.0097 / 0.0147 / 0.0159 | −10.11835 / −10.11902 / −10.12015 / −10.12050 | 3.5e-3 / 2.7e-3 / 1.4e-3 / 2.4e-4 |
| 0.10 | 0.01625 (S = 9/2) | −10.12079 | 0 |
| 0.15 | 0.0088 / 0.0163 / 0.0213 / 0.0238 | −10.14618 / −10.14803 / −10.14987 / −10.15097 | 7.9e-3 / 5.7e-3 / 3.6e-3 / 8.6e-4 |
| 0.15 | 0.0244 (S = 9/2) | −10.15171 | 0 |
| 0.20 | 0.0125 / 0.0200 / 0.0275 / 0.0322 | −10.18663 / −10.18918 / −10.19287 / −10.19572 | 1.4e-2 / 1.1e-2 / 6.9e-3 / 7.5e-4 |
| 0.20 | 0.0325 (S = 9/2) | −10.19603 | 0 |

## Unit test U10 (passed)

Both ends of the Task C interval at j = 0.2 were re-solved from scratch, without checkpoints.

| d | S | Re-solved E | Task C E | Difference |
|---|---|---|---|---|
| 0.01125 | 1/2 (levels, S²) | −10.186327697145138 | −10.186327697145 | 1.4e-13 |
| 0.0125 | 3/2 (levels, S²) | −10.186628949714953 | −10.186628949715 | 4.8e-14 |

U1–U9 passed in the first and second requests; the results are in `unit_tests.json`.

## Method and caveats

- **H1:** at every d, S^z = 1/2 and 3/2 were computed over all 9 momenta. S > 1/2 when E(3/2) ≤ E(1/2) + 1e-8. The starting interval was [0, 0.02] and the bisection resolution 0.0003125 (requirement: ≤ 0.0004). At both interval ends, S^z = 1/2, 3/2, 5/2 were rebuilt with eigenvectors over all momenta and S was confirmed with S².
- **H2:** the grid starts at the upper end of the first jump with Δd = 0.0025 and runs until S = 9/2 (stopping at 2j² at most). Every change of S was bisected down to 0.0003125 (requirement: ≤ 0.0005). At each d the lowest levels were found for S^z = 3/2 upward (all momenta) until the energy rises.
  - The S^z = 1/2 check was dropped at grid points partway through (j = 0.15, 0.2); the specification allows computing only the sectors that set S. This assumes S does not fall back to 1/2 after the first jump. It was checked with S^z = 1/2 (Γ, K, K′) at both ends of every threshold interval, and no warning appeared.
- **Quality:** all 933 rows converged, with ‖Hψ − Eψ‖ ≤ 1.0e-9.
  - The S² residual is at most 8.9e-7, at the lower end of the j = 0.10 first-jump interval (d = 0.004375). There the S = 1/2 and S = 3/2 levels in the same sector are only 6.8e-7 apart, so a small admixture is inevitable. The S value itself (0.5000) is unambiguous. Every other row has an S² residual ≤ 1.2e-7.
- **Sizes:** only 18 and 27. The exponent p is fitted over a narrow window, j = 0.10–0.20, and is not a thermodynamic-limit exponent.
- **Wall time:** 2026-11-18 23:35 to 2026-11-20 00:06 UTC (about 24.5 hours, including two restarts that resumed from cache with no loss), on 4 cores.
