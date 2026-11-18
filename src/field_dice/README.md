# field_dice — magnetic field and the honeycomb–triangle–dice interpolation

Stage-3 calculations for the PRB paper on the stuffed-honeycomb Heisenberg
antiferromagnet in a field (r = J'/J from the honeycomb side through the
triangular point r = 1 to the dice limit J/J' -> 0).

Geometry (sites, bonds, K18/C18 tori) is imported unchanged from
`src/workspaces/finite_size_18_27_36pilot/geometry.py`.

| File | Purpose |
| --- | --- |
| `lswt.py` | One-magnon saturation field h_sat(r) and LSWT stability of collinear UUD states |
| `windows.py` | Classical (LSWT) field window of the A-down and H-down UUD states vs r |
| `ed_sz.py` | ED in every S^z sector (H = J H_rim + J' H_hub), C2 inversion and primitive translation t |
| `validate.py` | Reproduces the verified S^z = 0 ground energies of `gs_18_scan.csv` (agreement ~1e-14) |
| `scan.py` | Full scan: lowest energy per sector; plateau-sector levels with C2 parity, hub magnetization, t eigenvalues |
| `analyze.py` | Plateau edges h_lo, h_hi, width and h_sat from the scan |
| `simplex.py` | S3 fundamental domain (eps, phi): h_sat = 9 S Jbar (1+eps), classical ground state, f_cl, 1/J triangle-inequality boundary |
| `three_j.py` | Three independent couplings: Lieb-Mattis edge check, partial-disorder wedge near the honeycomb vertex, resonance index f |
| `ed3.py` | Shared helpers: S^z sectors with H_AB, H_CB, H_CA split; lowest levels |
| `lswt3.py` | LSWT for collinear states with three couplings: stability window, zero-point energy |
| `windows3.py` | Classical plateau window h_pm = 3S[J_min +- sqrt((J1-J_min)(J2-J_min))] over the domain |
| `f_lswt.py` | Resonance efficiency f from LSWT vs ED in the Lieb-Mattis region |
| `wedge.py` | 18-spin partial-disorder wedge: jump asymmetry delta_c(j) -> data/field_dice/pd_wedge_18.csv |
| `ref27.py` | Independent unreduced N=27 reference energies for S^z >= 13/2 (unit tests of the 27-spin solver) |
| `scan_simplex.py` | 18-spin ED over 55 points of the fundamental domain -> data/field_dice/simplex_scan_18.csv |

Reproduce (about 2 min on 2 cores):

```bash
cd src/field_dice
python3 lswt.py        # P6: h_sat = S max(6J+3J', 9J')
python3 windows.py     # P7: classical windows
python3 validate.py    # unit test against the repository's 18-spin data
python3 scan.py        # writes data/field_dice/sz_scan_18.csv, plateau_levels_18.csv
python3 analyze.py     # writes data/field_dice/plateau_summary_18.csv
python3 simplex.py     # Q1, Q2, Q5 classical checks
python3 scan_simplex.py  # fundamental-domain scan (about 2 min)
```

Conventions: J = 1 (or J' = 1 for the dice-side rows with J < 1), S = 1/2,
h includes g mu_B, H = ... - h sum S^z. Plateau sector S^z = N/6 (= 3 for N = 18).
Results and their interpretation are tracked in the research-plan Claude Doc.

Next stage: `docs/guides/27SPIN_THREE_COUPLING_REQUEST_EN.md` (27-spin request; Korean summary alongside).
