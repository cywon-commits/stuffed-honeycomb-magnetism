# RUN_LOG — 36-spin Task I

## Environment

- claude.ai cloud container, 4 vCPU Intel Xeon @ 2.10 GHz (AVX-512), 15 GiB RAM, no swap
- g++ 13.3.0, `-O3 -march=native -fopenmp -std=c++17`; OpenMP 4 threads
- Python 3.13, numpy, scipy 1.18.1 (installed for the rotor fit, `src/field_dice/rotor_aniso.py`)
- Solver `src/ed36/ed_sym.cpp`, driver `src/field_dice/n36/task_i.py`, runner `src/field_dice/n36/run_i.sh`
  (resumable; every solver result is cached outside git)

## Commits

| stage | commit |
|---|---|
| driver, U13, I4 dimensions | `f385dd1` |
| run started (U11, U12, I1) | `f385dd1`; resumed after a container restart at 06:16 UTC (same commit) |

## Wall time and memory (I1, r = 0.8; pass 1 = eigenvalue, pass 2 = eigenvector + observables)

| S^z | sectors solved | largest sector | pass-1 wall time (sum) | pass 2 + observables | peak RSS |
|---|---|---|---|---|---|
| 8 | 9 (all, equivalence check) | 7,061,460 (K) | 1.5 h | 5 min | 0.53 GiB |
| 7 | 6 representatives | 16,689,036 (K) | 2.4 h | 14 min | 1.2 GiB |
| 6 | 6 representatives | 34,770,357 (K) | ≈ 4.5 h (one sector resumed after restart) | 27 min | 2.6 GiB |
| 5 | 6 representatives | 64,188,600 (K) | 9.5 h | 53 min | 4.9 GiB |

Total I1: 19.1 h (2026-11-20 06:16 – 11-21 01:21 UTC).
Typical cost: about 1.2 µs per state per H action (real Gamma sectors), and about 2× that for complex K sectors. Lanczos needed 64–160 steps.

Unit tests: U11 26 s, U12 61 min, U13 < 1 s, U14 < 1 s.

## I2 (r = 0.6)

| S^z | sectors solved | pass-1 wall time (sum) | pass 2 + observables | peak RSS |
|---|---|---|---|---|
| 8 | 9 (all, equivalence check) | 1.2 h | 5 min | 0.53 GiB |
| 7 | 6 representatives | 1.8 h | 11 min | 1.3 GiB |
| 6 | 6 representatives | ≈ 3.4 h (one K sector resumed after a container restart at 04:28 UTC) | 24 min | 2.7 GiB |
| 5 | 6 representatives | 8.3 h | 48 min | 5.0 GiB |

Total I2: 13.4 h of solver time (2026-11-21 01:21 – 17:51 UTC including the restart).
U14 rerun after I2: passes (8 states, max |zz_AB + t_AB − ⟨H_AB⟩/36| = 1.3e-14).
