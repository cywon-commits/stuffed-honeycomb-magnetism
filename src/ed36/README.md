# ed36 — symmetry-adapted matrix-free ED for N = 36 (and smaller) clusters

A new low-|Sz| solver for the spin-1/2 stuffed-honeycomb Heisenberg model

```text
H = J Σ_rim S_i·S_j + J' Σ_hub–rim S_i·T_R ,   J = 1, r = J'/J
```

built to answer whether the 36-spin Sz = 0 problem, previously judged
infeasible (8 GiB / 2-core environment, translation symmetry only), can be
solved on a modest machine. It uses the **same geometry module** as the earlier
work (`src/workspaces/finite_size_18_27_36pilot/geometry.py`, 36-site torus
`T = [[2,2],[-2,4]]`).

## Method

| ingredient | detail |
|---|---|
| symmetry group | full bond-type- and hub-preserving automorphism group of the cluster, found by graph search: 144 for N = 36 (12 translations × 12 point operations), 108 for N = 27, 48 for N = 12. Sz = 0 adds the global spin flip (×2). |
| sectors | for every momentum k: little group of k; if its point part is non-abelian (C6v at Γ, C3v at K) the cyclic subgroup (C6 / C3) is used so that **all** irreps are 1D and the sector set is complete (2D irreps appear as pairs of complex 1D sectors). Burnside's lemma gives every sector dimension; Σ dims = full C(N, N↑) is checked for N = 12, 18 (K18, C18), 27, 36. |
| representative search | hub sites occupy the top bits. Every symmetry maps hubs to hubs, so a table over the 2^12 hub patterns gives the minimal hub image and the few group elements realizing it (mean 1.86 of 288 for N = 36). Only those are applied to the 24 rim bits (3 byte-tables). |
| index lookup | sorted representatives, buckets keyed by (hub block, top rim bits) with ~2 states per bucket, then a short binary search. |
| matvec | pull form (no write races), OpenMP, real arithmetic when all characters are real. |
| eigensolver | Lanczos, eigenvalues every step, Ritz vector/residual when E0 settles; pass 2 rebuilds the Krylov vectors to form ψ and reports ‖Hψ − Eψ‖. Pass-1 checkpoint (`--ckpt PREFIX`, every `--ckevery` steps) allows restart after container loss. |

## Files

- `symmetry.py` — cluster, automorphisms, sectors, Burnside dimensions, spec export (`python3 symmetry.py N36 --outdir specs`).
- `ed_sym.cpp` — solver. Build: `g++ -O3 -march=native -fopenmp -std=c++17 -o ed_sym ed_sym.cpp`.
- `validate_small.py` — N = 12: union of all 36 dense sector spectra vs. the full unsymmetrized Sz = 0 spectrum.
- `validate_reference.py` — Lanczos vs. archived K18/C18 per-momentum energies and the verified 27-spin ground states.
- `results/` — validation JSON and benchmark logs.

Usage:

```bash
./ed_sym specs/N36_nup18/s000.spec --mode dim            # basis only
./ed_sym specs/N36_nup18/s000.spec --mode bench --r 0.3  # timed H-actions + Hermiticity
./ed_sym specs/N36_nup18/s000.spec --r 0.3 --tol 1e-9 --ckpt ck/s000_r0.3 --out gs.vec
```

## Validation

See `results/validation_N12.json` and `results/validation_reference_18_27.json`.

- **N = 12, r = 0.3 / 0.77 / 1.0:** the union of the 36 sector spectra equals the full 924-state Sz = 0 spectrum, max |ΔE| = 3e-14; Hermiticity error ≤ 3e-15; every sector dimension equals its Burnside count.
- **N = 18 (K18, C18), r = 0.1 / 0.3 / 0.5 / 1.0:** lowest energy in every one of the 6 momenta agrees with the archived NEW_18_SPIN_ED data to ≤ 1.8e-14; global ground energies agree to ≤ 1.1e-14; eigenvector residuals ‖Hψ − Eψ‖ < 1e-9.
- **N = 27 (Sz = 1/2, no spin flip), r = 0.2 / 0.2375:** global ground energies −10.185949128477 and −10.224890752095 agree with the verified 27-spin checkpoints to 1.4e-14 and 5e-15; the ground-state momentum is the reference K point (found at −k, i.e. the opposite Fourier sign convention).

## Measured N = 36 benchmark (this environment)

Environment: 4 vCPU Intel Xeon @ 2.1 GHz (AVX-512), 15 GiB RAM, no swap, g++ 13.3, OpenMP 4 threads.

| quantity | Γ, C6-trivial, flip-even sector (`s000`, contains A1) |
|---|---|
| dimension | 63,049,929 (= Burnside prediction) |
| table + basis build | 2.7 s |
| basis + lookup memory | 1.35 GiB |
| H-action | 68.3 s, 66.0 s |
| Hermiticity \|⟨x₂,Hx⟩ − ⟨Hx₂,x⟩\| | 8e-17 |
| peak RSS with 4 real vectors | 3.9 GiB |

Full log: `results/bench_N36_s000.log`. Every one of the 72 Sz = 0 sectors has
63.0 M states (Γ / K, real or complex) or more where the little group is
smaller; dimensions are listed by `symmetry.py N36`.

## Feasibility assessment (estimates, not measurements)

- **Memory is no longer the obstacle.** A real Γ sector Lanczos (3–4 vectors) needs ≈ 3.5–4 GiB; complex sectors roughly double the vector part (≈ 6–8 GiB). Both fit in 15 GiB. The translation-only approach needed ≥ 42 GiB.
- **Time is the constraint.** With ≈ 67 s per H-action, a Lanczos run of 150–300 steps plus the eigenvector pass costs ≈ 5–11 h per (sector, r) point on 4 cores. Covering all 72 Sz = 0 sectors at one r would take weeks here. A targeted study (lowest few sectors at a handful of r values, guided by the 18/27-spin results) is feasible over days; a dense r scan is still an HPC job, but now needs only ~16 GiB per process instead of 128–256 GiB.
- Larger total spin is cheaper: the Sz = S sectors are several times smaller than Sz = 0.
- The cloud container is reclaimed after inactivity; long runs must use `--ckpt` and be restartable.

## Not yet implemented

- S² / spin-multiplet identification and the observables required by `docs/guides/HPC_36SPIN_EXECUTION_REQUEST_EN.md` (⟨S_H²⟩, ⟨S_R²⟩, structure factors, …).
- Splitting Γ A1/A2 (and B1/B2) with the full C6v 1D characters, which would halve the most important sector to ~31.5 M states.
- Thick-restart / block Lanczos for several low-lying states per sector.
