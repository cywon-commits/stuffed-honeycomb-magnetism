# Stuffed Honeycomb Quantum Magnetism

Exact diagonalization (ED), quantum spin correlations, and dynamical spectra of the spin-1/2 stuffed-honeycomb Heisenberg antiferromagnet, interpolating between the honeycomb-host/hub regime and the isotropic triangular-lattice limit.

## Hamiltonian

```text
H = J sum_<i,j in rim> S_i . S_j + J' sum_<i in rim,R in hub> S_i . T_R
r = J'/J
```

The 3-site cell has two honeycomb rim sites (A,B) and one hub site (H). Reference numerical work includes 18-spin cluster scans, 27-spin ground-state and spectral calculations, and a preliminary 36-spin high-magnetization benchmark.

## Research questions

- How does a weakly polarized/disordered hub subsystem develop collective ferromagnetic correlations?
- How do increasing rim uniform-spin correlations compensate hub polarization?
- How do low-energy hub excitations couple to the rim's spin dynamics?
- How do local ABC-triangle doublets, dimer correlations and chirality responses connect to the triangular-lattice limit?

**Scientific caution:** Finite-cluster level crossings are not automatically bulk phase boundaries; finite spectral broadening is not physical magnon damping.

## Repository layout

- `src/` — preserved source code from the 18-, 27- and triangle-ED stages (see per-stage directories).
- `data/` — numerical CSV/JSON benchmarks and selected compact checkpoint material.
- `docs/` — HPC handoff, physics reports, and execution specifications.
- `scripts/` — environment and numerical verification helpers.
- `external_data/` — manifest and links for large binary eigensystem checkpoints and full historical archives.

The historical source packages contain some duplicated versions; they are retained by stage rather than silently mixed into one implementation. **The 36-spin low-spin solver is not implemented; the existing 36-spin code is a feasibility pilot.**

## Large artifacts

A complete archival package (including large eigenvector checkpoints) is maintained separately from git because storing hundreds of MB of binary checkpoints in the repository would be difficult to version and review. See `external_data/README.md`.

## Start here

Read `docs/README_START_HERE.md` and `docs/HPC_36SPIN_EXECUTION_REQUEST_EN.md` to reproduce the existing work or continue a 36-spin HPC calculation.

## Data provenance

Numerical files, validation results and scripts are those provided during the Quantum Magnetism research project through 2026-11-17. Please inspect each stage's README and convergence logs before drawing conclusions. This repository is research code, not a validated production software release.
