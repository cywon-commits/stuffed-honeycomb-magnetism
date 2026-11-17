# 27-spin ABC triangle: completed ED preanalysis + Codex-ready follow-up

Date: 2026-11-16

## What's completed

- Reused the numerical 27x27 SU(2)-scalar spin-correlation matrices from the verified Stage 1/2 and level-crossing ZIPs.
- Identified 54 ABC elementary triangular faces from the original geometry, 27 AB bonds and 54 hub-rim bonds.
- Calculated per-face bond correlations, pair singlet probabilities, local triangle doublet occupancy, local dimer order parameter `<D>`, `<chi^2>`, and symmetrized local dimer intensity for **11 verified coupling values**; included 2 excited competing branches near the two crossings as branch contrasts.
- Reconstructed E0 from bond correlations; for 0.2 / 0.5 / 1 matched original verified ED to available precision.
- Independently confirmed all reported spin-1/2 local operator identities on an explicit 8x8 Hilbert space.
- Prepared a three-panel, source-data-driven figure and an ED execution specification for missing higher-point correlations.

## What's NOT completed

No new 27-spin diagonalization has been run here; no 4-/6-spin **intertriangle** dimer/chirality correlations or local 3-site RDM from raw eigenvectors have been computed. These require the source eigenvectors and additional ED operators; request those from Codex using `Codex_Triangle_Resonance_Execution_Spec_2026-11-16.md`. An individual 3-spin local doublet occupancy is **not** proof of a coherent trimer resonance phase.

## Reproduce existing results offline

```
python check_local_identities.py
python analyze_triangle_existing.py --input-dir /directory/containing/original/zips
python plot_triangle_results.py
```

If original large ZIPs are absent, `analyze_triangle_existing.py` automatically reads the embedded snapshots in `reference_realspace/`. All printed values and plots can be reproduced with these small inputs. Archive manifest with per-snapshot SHA256 hashes is in `reference_manifest.csv`.

## Files

- `triangle_27spin_existing.csv` — complete grounded table (11 ground+2 competitor rows)
- `triangle_local_ED_figure.png`, `.svg` — publication style figure
- `analyze_triangle_existing.py` — full local analysis from verified real-space correlation matrices
- `check_local_identities.py` — spin-1/2 local algebra test
- `plot_triangle_results.py` — deterministic source-data plotting
- `reference_realspace/` — original 27x27 correlation matrices only (not ED eigenvectors)
- `Codex_Triangle_Resonance_Execution_Spec_2026-11-16.md` — missing ED tasks and numerical tests
- `Triangle_Resonance_Preanalysis_Report_2026-11-16.md` — physics interpretation and caveats
