# Stuffed-honeycomb static ED

Start with **physics_report.md** (concise results), **SPEC_REVIEW.md** (definition review), and **technical_report.md** (validation, symmetry and resources).

The Hamiltonian and production static calculations are matrix free. No spectral/dynamical calculation is invoked. Dense matrices occur only in independent small-cluster tests.

## Reproduce

Requirements: Python 3.12, requirements.txt, g++ with OpenMP. From the extracted folder:

```sh
OPENBLAS_NUM_THREADS=1 python run_static.py
```

Converged files are reused. `--force` reruns the sector solves and observables; runtime is measured in minutes and may vary substantially by hardware. The measured configuration uses ncv=32 for the 0.2 Gamma/reflection blocks and ncv=64 for the remaining blocks; run_static.py preserves these choices. The 600-second default per-block safety limit can be raised with `solve_sector.py --seconds ...`; timed-out calculations are never labeled converged.

Independent validation:

```sh
OPENBLAS_NUM_THREADS=1 python test_observables.py
OPENBLAS_NUM_THREADS=1 python test_cyclic.py
OPENBLAS_NUM_THREADS=1 python test_solver.py
```

Reuse the original verified isotropic vector:

```sh
OPENBLAS_NUM_THREADS=1 python analyze_state.py --jp 1 --primitive --file isotropic_primitive_ground.npz
```

Run one unequal-coupling block, with cell momentum k=(n1/3,n2/3), including its compatible cyclic character:

```sh
OPENBLAS_NUM_THREADS=1 python solve_sector.py --jp 0.2 --n 0 1 --cyclic --angular 0 --ncv 64
```

`cyclic_model.py` uses C6 rotations at Gamma (character ell=0,...,5), C3 rotations at cell (1/3,2/3), and a compatible reflection C2 at cell (0,1/3). `angular` labels the character index even for the reflection, where 0/1 mean +/- parity. Gamma ell=4,5 follow from complex conjugates of ell=2,1. `coverage.py` requires all nine representative blocks before declaring a global ground sector. `ground_j*_coverage.json` names the winning checkpoint.

For a winning cell/cyclic checkpoint:

```sh
OPENBLAS_NUM_THREADS=1 python analyze_state.py --jp 0.2 --file j0.2_n12_r0.npz
```

Use the file actually identified in ground_j0.2_coverage.json; this example is not a substitute for the coverage check. The program recognizes cyclic metadata and converts the vector to the correct cell basis before sublattice measurements.

Rebuild derived observables, tables and figures:

```sh
OPENBLAS_NUM_THREADS=1 python enrich_results.py
python build_report.py
OPENBLAS_NUM_THREADS=1 python plot_static.py
```

The primary matrices are SU(2)-multiplet and point-group averaged scalar dot products. CSV `convention` also distinguishes the selected Sz=+1/2 raw/connected scalar and longitudinal zz matrices. Cell momentum is always labeled by n1,n2 over three. Physical primitive K/K' are separate basis-phase combinations at cell Gamma; no sector is labeled M.

All ground searches and scalar measurements have N_up=14. The +/- operators are tested on small fully enumerated sectors; an additional full-size static global-spin ladder check at Jprime/J=0.5 uses N_up=15,16,17 solely to confirm the S=5/2 multiplet. It can be repeated with `OPENBLAS_NUM_THREADS=1 python verify_spin_multiplet.py`. No transverse spectral or continued-fraction computation is performed.
