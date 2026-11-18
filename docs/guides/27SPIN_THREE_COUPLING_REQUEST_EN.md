# Execution request: 27-spin ED of the three-coupling triangular antiferromagnet

Date: 2026-11-18. Target paper: PRB, "honeycomb–dice–triangle" S3 phase diagram.
Research plan (Claude Doc): https://claude.ai/artifact/EoRrJbnJmc8erKYMbU6eGj
Stage-3 code and 18-spin data: `src/field_dice/`, `data/field_dice/` (this branch).

## 1. Model and conventions

Spin-1/2 Heisenberg model on the triangular lattice split into three sqrt3 x sqrt3
sublattices A, B, C (C = "hub"; site index 3*cell + s, s = 0 A, 1 B, 2 C):

```
H = J_AB sum_{A-B} S.S + J_BC sum_{B-C} S.S + J_CA sum_{C-A} S.S - h sum S^z ,   all J >= 0
```

Each site has 3 bonds of each type touching its sublattice (coordination 6 when all J > 0).
Bond lists are exactly those of `geometry()` in
`src/history/hub_ground_state_crossings_verified/hub_ground_state_crossings/ed.py`:
rim bonds A_R-B_R, A_R-B_(R-a1), A_R-B_(R-a2) carry J_AB; hub bonds from C_R go to
A_R, B_R, A_(R+a1), B_(R+a1-a2), A_(R+a1-a2), B_(R-a2); **hub bonds to s = 0 carry J_CA,
hub bonds to s = 1 carry J_BC.** This is the only change to the model. The verified kernel
already takes one coupling per bond (`js` array), so `kernel.cpp` needs no change.

Parametrization used in the paper (S3 E representation, Jbar = 1):

```
J_k = 1 + 2 eps cos(phi - 2 pi k / 3),   k = 0 (AB), 1 (BC), 2 (CA)
phi = 0   : J_AB > J_BC = J_CA  (honeycomb-direction mirror; old r = J'/J = (1-eps)/(1+2eps) < 1)
phi = 60  : J_AB = J_BC > J_CA  (dice-direction mirror;      old r = (1+eps)/(1-2eps) > 1)
```

`src/field_dice/simplex.py: couplings(eps, phi)` returns (J_AB, J_BC, J_CA).
Cluster: N = 27, 3x3 cell torus (`ed.py`, L = 3), cell translations only
(`primitive=False`). Point-group reduction is optional and valid only on the mirrors.

S^z sectors: n_up = 14 (S^z = 1/2), 15 (3/2), 16 (5/2), 17 (7/2), 18 (9/2), 19 (11/2), 20 (13/2).
Dimensions per cell-momentum sector: about 2.23 M, 1.93 M, 1.45 M, 0.94 M, 0.52 M, 0.25 M, 0.10 M.
The 1/3 plateau is S^z = N/6 = 9/2.

## 2. Code change (small)

1. Add `geometry3(L, jab, jbc, jca)` that reproduces `geometry(L, jp)` bond by bond but
   assigns J_CA to hub bonds ending on s = 0 and J_BC to those ending on s = 1.
2. Let `Model` accept `couplings=(jab, jbc, jca)` and use `geometry3`. Keep `primitive=False`
   unless jab = jbc = jca.
3. Reuse the existing Lanczos with explicit residual check and the existing S^2 action
   (as in `cross_scan.py` / `observables.py`) for the total-spin label of reported states.

## 3. Unit tests (must pass before production)

| Test | Expected |
| --- | --- |
| U1. `geometry3(3, 1, r, r)` equals `geometry(3, r)` bond for bond | identical lists |
| U2. Verified 27-spin values, J = (1, r, r), S^z = 1/2, cell Gamma | r = 0.275: E0 = -10.290051743981 (S = 7/2); r = 0.35: -10.472017603235 (S = 7/2); r = 0.425: -10.732527324629 (S = 5/2, C6 character 3); source `data/history/hub_ground_state_crossings_verified/.../initial_points_table.md` |
| U3. High-S^z sectors, minimum over all 9 cell momenta | match `data/field_dice/reference_N27_highSz.csv` (independent unreduced sparse ED, 5 coupling points incl. off-mirror and J_CA = 0) to 1e-9 |
| U4. Saturation field E(27/2) - E(25/2) | equals 3S[sum J + sqrt(sum J^2 - sum_{i<j} J_i J_j)] = 9S Jbar (1 + eps) exactly |
| U5. Lieb–Mattis edge J = (1, 1, 0) and (1, 0.5, 0) | ground spin S = 9/2 (= N/6) |
| U6. Hermiticity and translation covariance of the reduced operator on random vectors | < 1e-12 relative |

## 4. Production tasks

All energies in units of the coupling scale stated in each row; report every eigenpair with
residual ||H psi - E psi|| < 1e-8 and S(S+1) from the S^2 action.

### Task A — bulk plateau on the honeycomb mirror (prediction P7')
- Couplings J = (1, r, r), r = 0.4, 0.6, 0.8.
- Sectors S^z = 1/2 ... 13/2, lowest level in every cell momentum (S^z = 1/2 and 3/2 may be
  limited to Gamma and (1/3,2/3), (2/3,1/3) after checking one r with all 9).
- Output plateau edges by Maxwell construction:
  h_lo = max_{m<9/2} (E_9/2 - E_m)/(9/2 - m), h_hi = min_{m>9/2} (E_m - E_9/2)/(m - 9/2).
- Decision: the rotor–hub finite-size formula predicts W_27 = (1/chi - 3r)/9 with chi = 0.14–0.16
  (W_27 ≈ 0.55–0.67 at r = 0.4, 0.48–0.60 at 0.6, 0.42–0.53 at 0.8). W_27 clearly above this band
  indicates a bulk quantum plateau for r < 1; inside it, finite-size tower only.

### Task B — upper plateau edge ratio (prediction Q7)
- Points: (eps, phi) = (0.1, 30°), (0.2, 30°); dice J = (1, 1, 0); mirror point J = (1.1, 1.1, 0.8).
- Sectors S^z = 7/2, 9/2, 11/2 (all 9 momenta).
- Classical window (selected UUD = sublattice opposite the weakest bond pointing down):
  h_pm = 3S [J_min ± sqrt((J_1 - J_min)(J_2 - J_min))].
- Report h_hi / h_+ and (h_lo - h_-). At N = 18 the ratio is 1.51–1.52 at every off-mirror point
  and h_lo - h_- = -0.08 ... -0.13. Report whether the ratio stays direction-independent at N = 27.

### Task C — partial-disorder wedge (prediction Q4); the most expensive task
- Couplings J = (1, j + d/2, j - d/2), j = 0.15, 0.20, 0.30.
- For each j, bisection in d on [0, 0.08] to tolerance 0.002 (about 6 evaluations per j).
- At each evaluation: lowest level in S^z = 1/2, 3/2, 5/2, 7/2, 9/2. For S^z = 1/2 use cell
  momenta Gamma, (1/3,2/3), (2/3,1/3) (the weak-coupling ground state is at (1/3,2/3)); run all 9
  momenta at d = 0 and at the final bracket. Ground spin S = largest S^z whose minimum equals the
  global minimum within 1e-8.
- Report d_c(j) = first d at which S > 1/2, and S after the jump.
- Decision: N = 18 gives d_c = 0.21 j (K18), 0.19 j (C18). Finite-size scaling d_c * N_H = const
  predicts d_c ≈ 0.14 j at N = 27 (0.021, 0.028, 0.042 for j = 0.15, 0.2, 0.3). A bulk wedge from the
  Seifert–Vojta effective hub model predicts d_c → 0.34–0.37 j^2 (0.0076, 0.0135, 0.030; see
  `src/field_dice/wedge_theory.py`). Classical non-collinear wedge: d = j^2. j = 0.15 and 0.2 separate
  the two predictions best; j = 0.3 is a consistency point.
- Correction (2026-11-18, later the same day): an earlier version of this request quoted 0.13 j^2. That
  number counted each hub pair once against a table that counts it twice and dropped the transverse term.

### Task D — resonance efficiency (prediction Q5)
- E_hex: J = (1, 0, 0), lowest level in S^z = 9/2 (honeycomb singlet + 9 free hub spins up).
- E0 at dice J = (1, 1, 0) and at (eps, phi) = (0.560, 30°) (0.97 of the domain edge eps_max = 0.577), S^z = 9/2, Gamma.
- f = E0 / ((J_AB + J_BC + J_CA) E_hex). N = 18: f = 0.82 (dice), 0.83 (phi = 30° edge).
  LSWT: 0.89. Classical: 1.0 and 0.98. Report f_27 and the direction of change with N.

## 5. Deliverables

- `data/field_dice/n27/` : one CSV per task with columns
  `task, J_AB, J_BC, J_CA, eps, phi_deg, n_up, Sz, k1, k2, E, S, H_residual, S2_residual`,
  plus `summary_n27.csv` with the derived quantities (h_lo, h_hi, d_c, f).
- `src/field_dice/n27/` : the modified solver, the unit tests (U1–U6) and the run scripts.
- `RUN_LOG.md` : commit hash, compiler flags, CPU/RAM, wall time and peak RSS per task.
- Unconverged levels must be labeled as such and never used in a derived quantity.

## 6. Resource notes

Peak memory per Lanczos is dominated by the S^z = 1/2 sector (2.23 M complex amplitudes,
about 36 MB per vector; the existing look-up tables need a few hundred MB). Previous runs with
2 cores / 8 GiB took ~2 s per H action at this size. Task C needs roughly
3 j × 6 d × (3 + 4 smaller sectors) Lanczos runs; Tasks A, B, D are lighter. Run A, B, D first,
then C. If time is short, run C at j = 0.15 and 0.2 (0.021 vs 0.0076, 0.028 vs 0.0135).

36-spin results are not assumed. If they appear later they are used only as a reference.
