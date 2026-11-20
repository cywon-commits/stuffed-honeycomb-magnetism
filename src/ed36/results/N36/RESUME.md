# Paused N = 36 run — how to resume

**Status (paused 2026-11-18T02:20Z by request):**

| | |
|---|---|
| run | r = 0.3, Γ, C6v A1, spin-flip z = +1 (singlet sector), 31,541,930 states |
| spec | `specs/N36_nup18_fulllittle_k0/s000.spec` (regenerate with `python3 symmetry.py N36 --full-little --k 0 --outdir specs`) |
| pass 1 (Lanczos) | **finished**: 226 iterations, E0 = −13.629646604664240 (E0/N = −0.37860129), Ritz residual 9.7e-7 |
| pass 2 (ground-state vector) | interrupted at step 110/226; its 757 MB checkpoint was only in container scratch and is not kept here |
| pair-orbit observables | not yet measured |
| log | `r0.3_Gamma_A1_zp.snapshot.txt` |

Note: this is the lowest state **in this sector only**. It is not yet a proven
36-spin ground state, because the other momentum and irrep sectors have not been computed.

## Saved checkpoint

`checkpoints/r0.3_Gamma_A1_zp.lz`: 3,640 bytes, the pass-1 Lanczos coefficients
(α, β for 226 steps, flag finished = 1).
SHA256 `e591fa7ec649bfac9af98408a2071567e2274a7d81467b5a1186936b84aa404a`.

With this file the solver skips pass 1 (≈ 2.5 h saved) and only rebuilds the
Krylov vectors from the same seeded start vector (pass 2, ≈ 226 H-actions ≈ 2.5 h
on 4 cores). Then it measures the pair-orbit correlations (≈ 5 min). The start
vector uses 4 fixed RNG streams, so the thread count of the new machine does not
matter. This was verified on K18: pass 1 with 4 threads, pass 2 with 1 thread
gave an identical E0 and residual, with norm 1.000000000000.

## Resume

```bash
cd src/ed36
g++ -O3 -march=native -fopenmp -std=c++17 -o ed_sym ed_sym.cpp
python3 symmetry.py N36 --full-little --k 0 --outdir specs
WORK=/path/with/2GB/free            # scratch for the pass-2 checkpoint and the 250 MB ground-state vector
mkdir -p "$WORK" && cp results/N36/checkpoints/r0.3_Gamma_A1_zp.lz "$WORK/"
setsid nohup ./run_n36.sh specs/N36_nup18_fulllittle_k0/s000.spec 0.3 r0.3_Gamma_A1_zp "$WORK" >/dev/null 2>&1 </dev/null &
tail -f results/N36/r0.3_Gamma_A1_zp.log        # expect "resumed pass-1 checkpoint at iteration 226 (finished=1)"
```

When it finishes:

```bash
python3 analyze.py N36 specs/N36_nup18_fulllittle_k0/pairs.txt results/N36/r0.3_Gamma_A1_zp.log
```

This gives S (from S²), ⟨S_H²⟩, ⟨S_R²⟩, ⟨S_H·S_R⟩, p_H, m_R,uniform², m_N²,
hub S(q) and the energy recomputed from G as a check. If the container is
interrupted again during pass 2, rerun the same `run_n36.sh` command. It resumes
from `$WORK/r0.3_Gamma_A1_zp.p2`, which is written every 10 steps.

Requirements: about 4 GiB of RAM, g++ with OpenMP, python3 with numpy.
