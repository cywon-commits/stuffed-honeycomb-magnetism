# Stage 3: 27-spin hub–rim cross spectral moments — direct ED verification

**Scope:** One momentum and coupling, J=1, J'/J=0.2, three-site cell probe k=(0,1/3). The 27-spin checkpoint `j0.2_n12_r0.npz` supplied with `hub_static_stages12_verified.zip` was reused. This is **not** a resolved dynamical spectrum and does not measure a physical linewidth.

## State, normalization and sector

- Ground state: cell k0=(1/3,2/3), S_tot=1/2, S^z_tot=1/2, E0/J=-10.1859491284774.
- The cell-momentum transfer k=(0,1/3) takes us to cell k_target=(1/3,0), still in S^z=1/2.
- Operator definitions (Nc=9): `h^z(k)=Σ_R exp(-ik·R) T_R^z/√Nc`, `n^z(k)=Σ_R exp(-ik·R)(S_A,R^z-S_B,R^z)/√Nc`, `m^z(k)=Σ_R exp(-ik·R)(S_A,R^z+S_B,R^z)/√Nc`; Φ(k) is the Fourier transform of the sum of the six rim spins neighboring each hub, with exact lattice geometry.
- Define |b_a> = O_a^z(k)|0>. The exact forward moments are `M0_ab=<b_a|b_b>`, `M1_ab=<b_a|(H-E0)|b_b>`, and `M2_ab=<(H-E0)b_a|(H-E0)b_b>`. Complex off-diagonal moments depend on Fourier phase convention. Reverse entries are complex conjugates.

## Direct 27-spin numerical results

| Component a, b=h | M0_ab | M1_ab / J | M2_ab / J² |
|---|---|---|---|
| n,h | -0.0233911435604 - 0.0135048830313i | -0.00173046899677 - 0.000999086741111i | -0.00418583795502 - 0.00241669467011i |
| m,h | -0.00883607564262 + 0.0153045319525i | -0.00173046899678 + 0.00299726022334i | -0.00341795313924 + 0.00592006849505i |
| Φ,h | -0.0498993704883 | -0.00692187598711 | -0.0144396973727 |
| h,h | 0.331770552871 | 0.00692187598711 | 0.00753533141040 |

The raw static and first moments from independent bond/correlation sum rules match the direct two-vector calculation. For h,h:

`M1_hh/M0_hh = 0.0208634429040 J` and `sqrt(M2_hh/M0_hh - (M1_hh/M0_hh)^2) = 0.149255468350 J`.

This is the variance of the **discrete spectral energy distribution**, not a lifetime or intrinsic linewidth. Because the centroid is small while its variance is comparatively large, low-energy spectral weight can coexist with a higher-energy tail.

The normalized magnitude of the n–h integrated overlap is `|M0_nh|/sqrt(M0_nn M0_hh) = 0.0614952` and for m–h is `0.0681689`. The small integrated cross-overlap does **not** imply small cross-spectrum at each individual excitation energy, because off-diagonal spectral residues can have different phases.

## Verification and limits

1. Constructed the source ground state in the 3-site cell translation basis from the independently verified C3-reduced checkpoint. The conversion preserved vector norm (~1 to 5e-15).
2. Applied each microscopic Fourier spin operator using the C++ seed kernel, thus avoiding a static-correlation-only inference.
3. Applied the *actual 27-spin Hamiltonian* in the target sector to the four spectral seed vectors, checking first and second cross moments directly.
4. Independent static-correlation sum-rule matrix is Hermitian and agrees with directly seeded M0, M1.
5. An exploratory Krylov/Lanczos run reached only ten steps before the computation limit, with unacceptably large lowest Ritz residuals (~0.067J and ~0.108J). Its apparent individual energies **must not** be quoted as physical pole positions. No converged frequency-resolved cross-spectral weights were obtained.

## Next full dynamic ED target

At J'/J=0.2 and cell probe k=(0,1/3), use a **block Lanczos** seed matrix from `n^z(k)|0>`, `m^z(k)|0>`, and `h^z(k)|0>` (optionally Φ), target momentum cell (1/3,0), Sz=1/2. Report discrete energies, complex cross spectral residues and independent residuals. Enforce M0, M1, M2 checks above, account for ground-state momentum partners, and keep plot broadening separate from a physical lifetime.

## Files

- `stage3_cross_moments_direct_27spin.csv` — complete Hermitian 4×4 M0/M1/M2 matrices.
- `cross_seed_validate_M2.py` — direct model/seed/H evaluation; run alongside the necessary source files and the ground-state checkpoint from the existing Codex ZIP.
