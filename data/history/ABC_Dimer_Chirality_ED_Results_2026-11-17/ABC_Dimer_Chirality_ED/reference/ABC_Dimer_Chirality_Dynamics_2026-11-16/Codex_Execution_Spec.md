# Codex execution specification — scalar dimer/chirality dynamics on 27-spin stuffed honeycomb

**Execution request, not a review request:** Perform new **frequency-resolved** ED of the scalar dimer and scalar chirality channels, prioritized by exact local symmetry and supplied 27-spin static spectral matrices. Reuse provided code, validated j=1 checkpoint, and static correlation data. Provide reproducible scripts and converged discrete poles with complex residues. Do not stop after analyzing this specification.

## A. Physics definitions and exact, mandatory tests

Original 27-spin Hamiltonian, **J=1, r=J'/J**, ABC unit triangle:

    H = sum_<A,B> S_A.S_B + r sum_<A,C>,<B,C> S_i.S_j.

There are 18 rim (A,B) and 9 hub (C) spins, 9 three-site cells, 54 overlapping elementary triangular faces, and 27 primitive triangular sites at r=1. All spins are 1/2. Use the original site indexing, periodic boundaries, and primitive-vs-cell momentum maps from `reference/ABC_triangle_ED/triangle/`.

For oriented face t=(A,B,C), define B_AB=S_A.S_B, etc., ω=e^{2πi/3},

    D_t = B_AB + ω B_BC + ω² B_CA,
    χ_t = S_A · (S_B × S_C),
    P_t = 1/2 - (2/3)(B_AB+B_BC+B_CA).

The **physical CCW** face orientation differs from ABC for the odd orientation index μ. Reuse `local_geometry.matrices` and `ensemble.unfold`, which handle these phases. Never silently replace χ_ABC by χ_CCW or vice versa.

Exact spin-1/2 local identities (verify with supplied standalone code):

    P²=P,  χ²=(3/16)P,
    (D†D+DD†)/2=(9/8)P,
    D†D-DD†=3 sqrt(3) χ,
    [χ,D]=-(sqrt(3)/2)D,   D²=0.

Each D/χ operator is a **global SU(2) scalar**. Hence from pure ground spin S_g, its matrix elements to *all* S_f≠S_g are exactly zero. This is crucial for spin-sector restrictions and for rejecting errors. At j=1, verified S_g=1/2; at j=.9 S_g=1/2 (C6 character 3). At j=.5 S_g=5/2, so target S_f=5/2 there. Scalars DO NOT have the same spin selection as S^z neutron scattering.

Under local C3: D -> ω² D, D† -> ω D†, χ -> χ. Under reflection exchanging two corners: D -> D† (up to fixed triangle-label phase), χ -> -χ; D_re and D_im form E, χ is the A2 pseudoscalar of C3v. Time reversal sends χ -> -χ, and leaves the two *Hermitian real-bond* components D_re and D_im invariant. The group-character selection rule is Γ_f ⊂ Γ_O ⊗ Γ_g at high-symmetry momenta. At generic Q there is usually no mirror little group and do not force D–χ cross residues to vanish.

The exact isolated-triangle doublet Hamiltonian in AB-singlet/triplet-flavor basis is

    Htri = (-1/4-r/2) I - (1-r)/2 τ_z,
    D=-3/4 τ_z + (3i/4) τ_x,   χ=(sqrt(3)/4) τ_y.

The isolated triangle **inelastic** poles at r<1 have ω/J=1-r; W(χ)=3/16, W(D)=9/16, plus W_el(D)=9/16. They are NOT the 27-spin poles. At r=1 the local states are exactly degenerate; no isolated-triangle nonzero tunneling frequency exists. Intertriangle interactions, which are not weak near r=1, determine the 27-spin spectrum.

## B. Input files and verified baseline

- `reference/ABC_triangle_ED/K_primitive.npz`: full 27-spin ground checkpoint at r=1, from previously verified triangular 27spin ED. The reused checkpoint has reported H residual ≈8.30e-9 J and S² residual ≈1.72e-8. Recheck these rather than claiming greater precision.
- `reference/ABC_triangle_ED/triangle/{ed.py,cyclic_model.py,sectors.py,observables.py,native.py,triangle_ops.cpp,kernel.cpp,local_geometry.py,ensemble.py,measure.py}`: verified source implementation, compile C++ kernels as needed.
- `reference/ABC_triangle_ED/triangle/raw_observables_j1.npz`: 27-spin static complex Gram matrices.
- `reference/ABC_triangle_ED/triangle/isotropic_ground_manifold.npz`: known ground K/K' zero-frequency spectral projection and excited-state covariance.
- `results/27spin_static_mode_targets.csv`: independently regenerated positive-energy M0 at each physical Q and two primitive up/down-face modes. Run `python code/static_mode_summary.py` to reproduce.
- `code/triangle_exact_dynamics.py` and `results/isolated_triangle_spectra.csv`: independent 8x8 theoretical calibration. Run it first.
- `reference/Codex_stage_triangle_report.md` and `reference/Codex_stage_triangle_README.md`: previous Codex report, ground manifold cautions.

**Important missing data:** The previous triangle-results ZIP included j=.9 ground-state energy/sector metadata but not the full j=.9 eigenvector. Do not assume a j=.9 vector is present. If a verified prior Codex checkpoint can be recovered, use it. Otherwise re-solve the j=.9 cell-Γ C6-character-3 state, target E/J=-14.205410071922 and S=1/2, and save H/S² residuals. Source routine `measure.source(.9)` expects `p0.900000000_n00_r3_e0.npz`; the current source bundle does not include that file.

## C. Prioritized dynamical experiments (execute sequentially)

**First, r=1, physical Q=(0,0), scalar χ channel.** The strongest connected physical Γ chirality mode has total M0=0.3681687647653339, and its known K/K' ground-pair weight is zero in this leading eigenvector. Use the **leading of the two up/down-face eigenchannels** and include its orthogonal partner. Compute χχ, along with D–χ cross response and relevant angular/face symmetry labels. A *subleading* χ Γ mode DOES carry 0.0184710714683 ground-pair M0; explicitly subtract it.

**Second, r=1, small physical finite Q** in the sixfold star of (8/9,7/9), equivalently (1/9,2/9) up to sign/point-group operation. The leading dimer eigenchannel has inelastic M0=2.75762542346748 and **no known K/K' ground-pair contribution**. Compute D†D and D D†, D–χ and χχ at the same physical Q. Include both up/down-face modes, do not renormalize the 2.757625 figure after selecting a different q convention.

**Third, r=1, physical K=(1/3,2/3) and K'=(2/3,1/3).** These fold to the three-site-cell Γ yet are *different from physical Γ*. The known K/K' ground pair produces a **0.120894263834** static zero-energy dimer contribution in one subleading face mode. Remove it explicitly before interpreting low-energy peaks. Compare two momentum partners.

**Fourth, r=.9, optionally r=.7** with the cell symmetry appropriate for r≠1, then r=.5 if feasible. At nonisotropic r the 27 primitive translations are NOT symmetry operations; only 9 three-site-cell translations apply. Reuse all six face/orientation fields at each cell q and **do not use the primitive-only 27-Q block diagonalization** except at r=1. In r=.9 evaluate cell Γ and first accessible finite cell-q sectors. Track whether spectral peaks move and whether their D/χ overlap changes.

## D. Spectral definition and algorithm

For each selected ground |g> in the fixed-Sz=1/2 sector, and symmetry-covariant scalar operators O_a(Q) (a=D,D†,χ with up/down-face labels), construct raw complex seeds

    b_a(Q) = O_a(Q)|g>,
    A_ab(Q,ω)=Σ_f <g|O_a†(Q)|f><f|O_b(Q)|g> δ(ω-E_f+E_g),
    M_p^{ab}(Q)=<b_a|(H-E_g)^p|b_b>.

Compute M0, M1, M2 from **actual H action on the seed vectors**, not numerical integration of a broadened plot. M0 must match supplied exact static Gram matrices within reported numerical precision (target absolute 1e-6 or tighter). All moment and residue matrices must be Hermitian and positive semidefinite; complex cross residues are essential and cannot be reconstructed from magnitudes alone.

At r=1, source K and K' are a fixed-Sz ground orbital pair. For each source and transferred Q, construct P_G from *all verified ground states* in the target momentum and spin sector, then compute (1-P_G)b_a. Average separate source spectra only AFTER correctly projecting their poles, with symmetry weights specified. The zero-frequency part is a separate deliverable; do not throw it away or call it an inelastic resonance. Verify the decomposition M0=M0_elastic+M0_positive_energy against `isotropic_ground_manifold.npz`, including at K and primitive Γ.

For each Q, use symmetry-reduced matrix-free target Hamiltonian, spin and point-group labels if available. Prefer stable **thick-restarted block Lanczos** seeded with several independent O_a|g> vectors, or a selective targeted eigensolver with independent converged Ritz residuals. With 2 CPU and 8 GiB resource constraints, limit block size and store basis on disk or restart as necessary; a naive all-Q 200-step full block basis can exhaust memory. Converged isolated poles must satisfy ||(H-E_f)|f>||<1e-8 J (or report looser if source residual prevents), and scalar selection checks |<f|(S²-Sg(Sg+1))|f>| and transition forbidden amplitudes <1e-7. For degenerate multiplets, **sum the residue matrices across the degenerate eigenspace**; individual eigenvector residues are basis dependent.

For each converged pole export omega/J, S_f, target momenta (physical and cell), irreducible character/reflectivity when defined, energy residual, M0 contribution, **complete complex residue matrix** and overlap with the seeds. Compute accumulated weight versus omega, missing moment remainders and positive spectral bounds. Clearly separate discrete finite-cluster lines, artificial visualization broadening η, and genuine linewidth (NOT inferred from this calculation).

Scalar Wigner-Eckart: S_f=S_g. If an apparent scalar excitation has S_f=3/2 from the S_g=1/2 ground at j=1, this is a geometry, tensor-index, operator-order, or numerical bug.

## E. Exact sum rules and noncommuting channels

At each *single oriented triangle* (not blindly at the same Fourier Q):

    ∫ dω [A_{D†D}(ω)-A_{DD†}(ω)] = 3sqrt(3)<χ>.

**Include ω=0!** The positive-frequency difference alone equals RHS minus the difference in elastic weights. For the selected original K state at j=1 the reported |<χ_CCW>|=0.01849477979; hence absolute circular-dimer integrated asymmetry is 0.09610169481. K' reverses the sign; the symmetry-restored ensemble has zero local mean asymmetry. With Fourier fields, sum/weight over Q and restore oriented geometry before using this local identity.

Moments: M0^ab=<b_a|b_b>, M1^ab=<b_a|(H-Eg)|b_b>, M2^ab=<(H-Eg)b_a|(H-Eg)b_b>. For a non-Hermitian D do NOT replace M1 by one-half of a double commutator without including the opposite-order contribution. Verify Gram matrices and all local identities before Lanczos.

If a cutoff Ω is used, do not infer remaining spectral weight from a single centroid. Use exact M0,M1 positivity and report rigorous residual bounds wherever available.

## F. Required outputs and phase-claim guardrails

- `run.log`, per-sector dimension/H-call/time/RAM record, exact input SHA256, compiled kernels, source revision.
- `static_moments.csv/json`, `poles.csv`, `residues.npz` (complex), `ground_projection.csv`, `mode_convergence.csv` including direct eigenvector residuals, `validation.json`.
- Publication-ready plots: positive-frequency spectra of χ at physical Γ, D at small finite Q, D and χ at physical K/K', exact ground zero-energy weight separately; r=.9 comparison if completed.
- Short physics report addressing whether weight concentrates into isolated low-energy flavor modes, broad high-energy excitations or both; spin and point-group selection rules; role of 120-degree magnetic correlations; avoid interpreting static doublet purity as coherent tunneling or finite-size discreteness as a lifetime.
- Complete reproducible ZIP package. If full resolution is impossible, give fully converged partial poles, remaining weight and exact sum-rule bounds; do not label unconverged Lanczos quadrature Ritz nodes as physical excitation energies.

**Start immediately with exact small-cluster tests and the r=1 χ(physical Γ) experiment.**
