# ABC trimer dimer/chirality dynamical ED handoff — 2026-11-16

**Codex: run the calculations**, not just review. Start with `Codex_Execution_Spec.md`.

Contents:
- `Theory_and_Preanalysis_Report.md`: completed analytic/group-theoretic findings, independently calculated 8×8 triangle dynamics and recomputed 27-spin static spectral targets.
- `Codex_Execution_Spec.md`: exact Hamiltonian, operator definitions, SU(2)/point-group selection rules, physical/three-site-cell momentum conventions, prioritized 27-spin spectral runs, strict validation criteria and deliverables.
- `code/triangle_exact_dynamics.py`: no imported repository dependencies; deterministic exact numerical 8×8 triad test generating spectral moments and cross residues.
- `code/static_mode_summary.py`: derives 27 physical-Q 2×2 face-channel inelastic M0 matrices from the supplied static data, requiring numpy and scipy and repository original `ensemble` module.
- `results/`: generated numerical CSVs and local-algebra validation JSON.
- `figures/`: preview plots of isolated triangle spectral gap and 27-spin inelastic M0 target strengths.
- `reference/ABC_triangle_ED/`: verified original ED kernels, geometry and original j=1 K ground checkpoint. The checkpoint is **not** newly generated, and its published residual should be reverified.
- `reference/Codex_stage_triangle_report.md`: previous complete 27-spin static analysis.

Reproduce quick tests at root:

    OPENBLAS_NUM_THREADS=1 python code/triangle_exact_dynamics.py
    OPENBLAS_NUM_THREADS=1 python code/static_mode_summary.py

**Important:** The j=.9 ground eigenvector is not included; use your previous Codex workspace or re-solve with `CyclicModel` and check E/J=-14.205410071922, S=1/2, Gamma/C6 character 3. More detailed instructions in execution spec.

The local isolated-trimer gap Δ=J-J' is **not** evidence of a 27-spin soft mode or coherent resonance phase. A finite-spin ED scalar line is not a physical lifetime/damping rate.
