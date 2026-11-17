# Stage 3 preanalysis after 27-spin level-crossing search

**Date:** 2026-11-15 (Asia/Seoul)  
**Model:** spin-1/2 stuffed honeycomb Heisenberg, 3×3 three-site cells, 27 spins, J=1, J'/J variable.  
**Source:** `hub_ground_state_crossings_verified.zip` (Codex ground states and 27×27 real-space spin correlations).  
**Status:** New static-data-derived **exact spectral zeroth and symmetrized first moments**, plus symmetry selection rules. No frequency-resolved dynamical poles or intrinsic linewidths computed in this step.

## 1. Ground state crossings and selection rules

| observed J'/J bracket | low coupling ground | high coupling ground | single-spin selection |
| --- | --- | --- | --- |
| [0.256250, 0.260938] | S=1/2, cell k=(1/3,2/3), C3 character 0 | S=7/2, cell Γ, C6 character 0 | **Forbidden**, ΔS=3 exceeds 1, at any momentum |
| [0.387500, 0.392188] | S=7/2, Γ, C6 char 0, even reflection | S=5/2, Γ, C6 char 3, odd reflection | **Allowed by spin and spatial symmetry** for rim-staggered n(Γ); forbidden for m(Γ), h(Γ) |

The energy differences between the two competing branches in the second bracket are 0.00575611693 J (at 0.3875) and 0.00382047818 J (at 0.3921875), respectively. These are *branch gaps*; no spectral matrix element or definitive excitation threshold is claimed. A vanishing branch gap is not automatically a visible spectral soft pole. Single spin probes have |ΔS|≤1.

At cell Γ, the C6 rotation swaps A/B and changes n=A−B by a minus sign, while m=A+B and h remain invariant. Therefore n(Γ) transforms as angular character 3, m and h as character 0. The total spin m(Γ)+h(Γ) is conserved. For a nondegenerate orbital Γ multiplet, n(Γ) has no ground-spin elastic weight by inversion; the m/h Γ channels generally do.

## 2. Operator and moment definitions

For cell momentum k=(k1,k2), sites i=(x,y,s), Nc=9, define c^a_{i}(k)=c_s^a exp[-2πi(k1 x+k2 y)]/√9, with c^n=(1,−1,0), c^m=(1,1,0), c^h=(0,0,1). Let O_a^μ(k)=∑_i c_i^a S_i^μ. The SU(2)-summed equal-time correlation and symmetrized first moment are

\[
M^{(0),\mathrm{dot}}_{ab}(k)=\sum_{ij}c_i^{a*}(k)c_j^b(k)\langle\mathbf S_i\cdot\mathbf S_j\rangle,
\]
\[
M^{(1),\mathrm{dot}}_{ab}(k)=-\sum_{\langle ij\rangle} J_{ij}\langle\mathbf S_i\cdot\mathbf S_j\rangle [c_i^{a*}-c_j^{a*}][c_i^b-c_j^b].
\]

These are positive-frequency spectral moments for a symmetry-restoring ensemble where forward and reverse sectors coincide; otherwise the double-commutator expression must be interpreted as a symmetrized first moment and the direct Lanczos forward moment checked. At fixed S^z, the longitudinal f-sum replaces the bond dot correlator by (Gdot−Gzz)/2. For Γ, the inelastic m/h moment matrix has entries proportional to [[1,−1], [−1,1]], and the n–h first moment and equal-time overlap are zero by inversion. A nonzero n–h response at generic cell momentum is allowed.

The raw Gamma hub equal-time moment contains an SU(2) ground multiplet elastic contribution. Its exact scalar projection is

\[
M^{(0),\mathrm{elastic}}_{hh}(\Gamma) = \frac{\alpha_h^2S(S+1)}{N_c},\quad
\alpha_h=\frac{N_c[M^{(0)}_{hh}(\Gamma)+M^{(0)}_{mh}(\Gamma)]}{S(S+1)}.
\]

This expression assumes a unique orbital ground sector at the chosen cell momentum. The remaining **inelastic static weight** is raw minus this elastic part. For fixed-M longitudinal spectra, subtract the **longitudinal connected/disconnected** part separately; do not apply the full three-component elastic fraction to a single zz component.

## 3. New numerical results from verified static ED

All numbers below use the 27×27 correlation matrices in the crossing ZIP, with energies independently checked against the reported ground energies. `M1` is SU(2)-summed unless `zz` is indicated.

### 3A. First crossing: hub finite-q spectral centroid

Hub probe q=(0,1/3) (cell reciprocal coordinates):

| J'/J | S_ground | M0_hh(q) | M1_hh(q)/J | M1/M0 (J) |
| ---: | ---: | ---: | ---: | ---: |
| 0.237500 | 1/2 | 0.993118 | 0.030206 | 0.030415 |
| 0.256250 | 1/2 | 0.991823 | 0.035731 | 0.036026 |
| 0.260938 | 7/2 | 0.506760 | 0.068244 | 0.134668 |
| 0.275000 | 7/2 | 0.507356 | 0.072549 | 0.142994 |

The centroid changes sharply (~3.74× between the crossing bracket endpoints), but this **is not a single-pole excitation energy**. The hub q-probe cannot directly connect the two crossing branches because ΔS=3; the change reflects distinct source ground states. For q=(0,1/3) at low-spin 0.2375, source k0=(1/3,2/3) and target k=(1/3,0). The inelastic hub centroid 0.030415J provides a variational upper bound on at least one hub-accessible excitation in that sector, provided zero-energy ground-partner weight is absent, as in the saved source orbital ground sector.

At 0.2375 and q=(0,1/3), the three-component n–h cross moments are

- M0_nh = −0.0851753091−0.0491759877 i (complex Fourier phase convention);
- M1_nh/J = −0.00755152692−0.00435987610 i.

They verify that finite-q dynamical mixing is **allowed and has nonzero integrated moment**, but the energy-resolved spectral residues are not yet known. The fixed-M longitudinal versions of these moments are in the CSV.

### 3B. Gamma hub elastic versus inelastic weight

| J'/J | hub Gamma raw M0 | exact ground multiplet elastic M0 | inelastic M0 |
| ---: | ---: | ---: | ---: |
| 0.256250 | 0.104274 | 0.082738 | 0.021536 |
| 0.260938 | 2.708405 | 2.563446 | 0.144959 |
| 0.387500 | 2.672621 | 2.519172 | 0.153449 |
| 0.392188 | 2.553266 | 2.198518 | 0.354749 |

The ~26-fold raw hub Γ increase at the first crossing is **largely ground-multiplet weight**, not a 26-fold increase in inelastic susceptibility. Nevertheless, the inelastic equal-time hub Γ weight also grows by ~6.7×. This is an integrated (equal-time) effect, not a lifetime.

### 3C. Second crossing: Γ rim Néel spectral moment and soft-branch target

| J'/J | S_ground | M0_nn(Γ) | M1_nn(Γ)/J | M1/M0 (J) | opposing branch gap (J) |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.350000 | 7/2 | 8.567346 | 4.365352 | 0.509534 | — |
| 0.387500 | 7/2 | 8.549285 | 4.370004 | 0.511154 | 0.005756 |
| 0.392188 | 5/2 | 8.081627 | 4.104874 | 0.507927 | 0.003820 |
| 0.425000 | 5/2 | 8.058194 | 4.113744 | 0.510504 | — |

Even while an n(Γ)-allowed opposing level approaches close to the ground state, the n-channel **spectral centroid stays ~0.51J**. Thus a low-energy pole may have a small or large residue relative to the abundant high-energy weight; the first two moments alone do not determine its visibility. A resolved transition matrix element \(|\langle f|n^z(\Gamma)|0\rangle|^2\) is the crucial next observable. Γ m/h cannot couple the specific char0↔char3 crossing states.

### 3D. Suggested exact longitudinal acceptance tests for Stage 3

| source J'/J | probe q | target sector | fixed-M M0_zz | symmetrized M1_zz/J |
| ---: | --- | --- | ---: | ---: |
| 0.2375 | h^z(0,1/3) | cell (1/3,0), Sz=1/2 | 0.33103929 | 0.01006870 |
| 0.35 | n^z(Γ) | cell Γ, Sz=1/2, C6 char3 | 3.21077537 | 1.41363196 |
| 0.425 | n^z(Γ) | cell Γ, Sz=1/2, C6 char0 | 3.11184094 | 1.32047247 |

The source pointgroup-averaged static data fix the symmetrized moment; direct *forward* Lanczos moments must be checked independently when source state is not an invariant ensemble. The 0.3875 and 0.3921875 near-crossing values for n^z(Γ) are (M0zz, M1zz)=(3.195485,1.416573J) and (3.118946,1.316844J), respectively.

## 4. Recommended Stage 3 execution order

1. **Near second crossing first (high scientific specificity):** at J'/J=0.3875 and 0.3921875, compute the *lowest opposite-character competing eigenstate* and the exact n^z(Γ) transition residue from the appropriate ground state; report direct residual, target S=5/2 or 7/2, C6 character, reflection parity, and complex matrix element. The transition is allowed but its intensity is unknown. If the supplied eigenvectors for both levels at the same coupling are absent, recompute and verify them; do not infer weights from eigenvalues alone.
2. **Weak coupling and feedback control:** at J'/J=0.2375 (or the already verified 0.2), probe q=(0,1/3), source cell k0=(1/3,2/3), target (1/3,0), using h^z, n^z, m^z. Record spectral poles and complex cross residues; test exact M0/M1 moments. To isolate the feedback channel, retain the rim-uniform (m) as well as rim-staggered (n) operator.
3. **Crossing 1 comparison:** if budget allows, compare 0.25625 and 0.2609375 finite-q hub spectra; do not expect a direct dipolar pole between the crossing branches (ΔS=3).
4. **Conservation/elastic control:** at Γ deflate the full ground multiplet as appropriate. For inelastic (m,h) spectral matrix check m+h is null for every positive-energy pole. A Lorentzian broadening parameter is not a physical linewidth. Do not infer thermodynamic Goldstone damping from 27-spin discrete energies.
5. **Only then** move toward weak-coupling rim Goldstone feedback by comparing different J' and, ideally, system sizes; a microscopic irreducible χhub is not simply the full coupled-system χhh.

## 5. Provenance and validation

- Source: `hub_ground_state_crossings_verified.zip`, ground states identified by the `obs_p*_realspace.npz` and corresponding `obs_p*.json` files.
- Reproduction: `stage3_moment_targets.py`, no new eigensolver. Its numerical checks include E=∑bond J·Gdot within 5×10⁻⁹J (observed considerably tighter), Hermiticity, positive-semidefiniteness of M0 and M1, inversion-protected Γ n/h zeros, and Γ m+h f-sum null.
- CSVs: `stage3_targets_summary.csv` and `stage3_targets_moments.csv` (full 3×3 n,m,h moments including complex entries at Γ, (0,1/3), (1/3,2/3), both scalar and fixed-M longitudinal first moments).
- Limitations: Input static correlations are Codex's already verified ground states; this is **independent algebraic postprocessing**, not a second independent ground-state ED solver. Not a resolved dynamical spectral calculation. Selection-rule-allowed ≠ finite residue, and finite-size crossing ≠ thermodynamic phase transition.
