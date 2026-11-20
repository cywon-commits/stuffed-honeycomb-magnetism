# RUN_LOG — 27-spin three-coupling ED (requests 1: Tasks A–D, 2: Tasks E–G)

## Environment

- Machine: claude.ai cloud container, 4 vCPU Intel Xeon @ 2.10 GHz (AVX-512), 15 GiB RAM, no swap
- Compiler: g++ 13.3.0 (Ubuntu 24.04), flags `-O3 -march=native -fopenmp -std=c++17`
- Python 3.13 + numpy (scipy not available; the verified `ed.py` imports scipy, hence the C++ solver below)
- Solver: `src/ed36/ed_sym.cpp` (symmetry-adapted matrix-free Lanczos, OpenMP), driven by `src/field_dice/n27/{core,tasks,tests}.py`
- Run scripts: `src/field_dice/n27/run_all.sh` (tests + tasks A–D), `run_second.sh` (E alongside C with 2 threads; G, F after C)
- Threads: 4 (OMP), except Task E: 2 threads while Task C used 4

## Commits

| stage | commit (branch `claude/stoic-ptolemy-l3eac5`) |
|---|---|
| request 1 runs started (tests U3–U5, tasks A–D) | `eb17f8a` |
| task C resumed after container restart | `878d949` |
| second request code (E–G, U7–U9) | `4f6dacf` (E), `abcf466` (G, F) |
| third request (Task H, U10) | `e06d9fa` (start), `325d5e8` (H2 S^z = 1/2 skip) |

Every solver result is cached (JSON lines, outside git), so reruns of the same scripts reproduce the CSVs without recomputation.

## Wall time and memory

| item | wall time | notes |
|---|---|---|
| U1, U6 | 3 s | |
| U3 + U4 | 193 s | high S^z, 5 coupling points |
| U2 | 494 s | 3 points, S^z = 1/2, Gamma |
| U5 | 2322 s | 2 Lieb–Mattis points, S^z = 1/2 … 13/2, all momenta |
| U7, U8, U9 | < 5 s | |
| Task A | 5766 s (1.6 h) | r = 0.4, 0.6, 0.8 |
| Task B | 1548 s | 4 points |
| Task D | 127 s | |
| Task C | ≈ 10.2 h | 05:27–15:41 UTC incl. one container restart (resumed from cache) |
| Task E | 932 s | 2 threads, concurrent with Task C |
| Task G | 2611 s | |
| Task F (r = 0.2, 0.3, 0.5) | 3387 s | |
| Task F (r = 0.7, 0.9) | 4383 s | |
| U10 | 94 s | third request |
| Task H (H1 j = 0.10, 0.125, 0.175; H2 j = 0.10, 0.15, 0.20) | ≈ 24.5 h | 2026-11-18 23:35 – 11-20 00:06 UTC, two container restarts resumed from cache; H2 grid points after 12:20 UTC skip S^z = 1/2 |

Peak resident memory per solver process: ≤ 106 MiB (largest sectors: 1.12 M states). Solver runs: 1524 pass-1 (eigenvalue), 450 pass-2 (vector, residual, S^2), 17 observable runs.

## Convergence policy

Pass-1 Lanczos to Ritz residual 1e-9 (max 2000 steps); every reported eigenpair rebuilt and checked: ||H psi − E psi|| ≤ 1.0e-9 for all rows of task_A … task_H; S from the explicit S^2 action (S2_residual ≤ 1.04e-8). No unconverged level occurs in any derived quantity (all rows `converged = True`).
