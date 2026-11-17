# Stage 3 preanalysis: exact frequency moments from 27-spin static ED

**Model:** 27 spin-1/2 stuffed-honeycomb Heisenberg model on 3×3 unit cells, \(J=1\).  
**Inputs:** `hub_static_stages12_verified.zip`, specifically the `static_j*_realspace.npz` correlation matrices, `static_j*.json` and the geometry in `ed.py`.  
**Status:** New *analytical sum rules and numerical moment evaluations*; **not** a dynamical Lanczos spectrum, linewidth, or thermodynamic extrapolation.

## 1. Definitions and exact first-moment identity

Let \(O^a_c(\mathbf k)=\sum_i c_i(\mathbf k) S_i^a\), where \(a=x,y,z\), and use \(c_{R,H}=e^{-i\mathbf k\cdot R}/\sqrt{9}\) with \(c_{R,A}=c_{R,B}=0\) for the hub operator. The static zeroth moment of the SU(2)-averaged three-component response is

\[
M_0^{\rm dot}(\mathbf k)=\sum_{a=x,y,z}\langle O^a_c(\mathbf k)^\dagger O^a_c(\mathbf k)\rangle.
\]

In an ensemble restoring inversion / time reversal, its positive-frequency first moment is

\[
M_1^{\rm dot}(\mathbf k)=\frac12\sum_{a=x,y,z}\langle[O^{a\dagger}_c(\mathbf k),[H,O^a_c(\mathbf k)]]\rangle
=-\sum_{\langle ij\rangle}J_{ij}\langle\mathbf S_i\cdot\mathbf S_j\rangle|c_i-c_j|^2.
\]

For a fixed-\(S^z\) longitudinal response, **do not assume** \(G_{zz}=G_{\rm dot}/3\) if \(S>1/2\):

\[
M_{1,zz}^{\rm sym}(\mathbf k)=-\frac12\sum_{\langle ij\rangle}J_{ij}\big[\langle\mathbf S_i\cdot\mathbf S_j\rangle-\langle S_i^z S_j^z\rangle\big]|c_i-c_j|^2.
\]

The double commutator directly gives the average of the forward and conjugate first moments; for a symmetry-averaged state the corresponding spectral moments coincide. A directly seeded Lanczos calculation should check \(\langle b|(H-E_0)|b\rangle\) independently.

Because there are **no hub–hub exchange bonds**, \(c_i-c_j\) along every hub–rim bond has unit modulus \(1/\sqrt{9}\) regardless of \(\mathbf k\). With 54 such bonds,

\[
\boxed{M_{1,hh}^{\rm dot}(\mathbf k)=-6J'\,G_{HR,\rm bond}},\qquad
\boxed{M_{1,hh}^{zz,\rm sym}(\mathbf k)=-3J'\,[G_{HR,\rm bond}-G_{HR,\rm bond}^{zz}]}.
\]

Both are **exactly independent of cell momentum** for this model and ensemble. These formulas are a powerful code acceptance test.

For \(\Gamma\), \(\mathbf S_{\rm tot}=\sqrt{9}(\mathbf m+\mathbf h)\) is conserved. Hence for every \(\Omega>0\)

\[
\mathcal A_{mm}=\mathcal A_{hh}=-\mathcal A_{mh}=-\mathcal A_{hm},\qquad
M^{(1)}_{(m,h)}(\Gamma)=\mu\begin{pmatrix}1&-1\\-1&1\end{pmatrix}.
\]

This applies separately to longitudinal and transverse channels with the appropriate normalization, after separating zero-energy ground-manifold transitions.

## 2. Numerical first moments: SU(2)-averaged three-component response

| \(J'/J\) | Hub–rim \(\langle\mathbf T\cdot\mathbf S\rangle\) | \(M_1^{\rm dot}/J\) at **every** cell \(\mathbf k\) |
|---:|---:|---:|
| 0.2 | -0.017304690 | **0.020765628** |
| 0.5 | -0.087482000 | **0.262446000** |
| 1.0 | -0.186740412 | **1.120442473** |

For the hub's sixfold generic momentum star represented by \(\mathbf k=(0,1/3)\), there is no coupling between ground-sector momentum partners at \(J'/J=0.2\). Thus \(M_0=S_h(\mathbf k)\) has no zero-energy ground-manifold contribution. The exact spectral centroid and the variational bound on the *lowest excitation with nonzero hub spin-operator matrix element* are

\[
\bar\omega=\frac{M_1}{M_0},\qquad 0<\Delta_{\rm hub,accessible}(\mathbf k)\le\bar\omega.
\]

| \(J'/J\) | \(S_h(0,1/3)=M_0\) | \(M_1/J\) | \(\bar\omega/J\) |
|---:|---:|---:|---:|
| 0.2 | 0.995311659 | 0.020765628 | **0.020863443** |
| 0.5 | 0.532603439 | 0.262446000 | **0.492760619** |
| 1.0 | 0.624625244 | 1.120442473 | **1.793783527** |

At \(J'/J=0.2\), \(J'^2/J=0.04J\), so the spectral centroid at the hub peak is **\(0.5216J'^2/J\)**. By positivity of the spectral measure, Markov's inequality yields

\[
\frac{\text{hub spectral weight below }\Omega_c}{\text{total hub spectral weight}}
\ge 1-\frac{\bar\omega}{\Omega_c},
\]

whenever \(\Omega_c>\bar\omega\). Consequently at least **47.8%** of the \((0,1/3)\) hub spectral weight is at energies below \(0.04J\), and at least **79.1%** lies below \(0.10J\). These bounds follow from the exact two moments, **not** from assumed spectral broadening.

## 3. Longitudinal moments: directly relevant to the Stage 3 first run

Use \(O_h^z(\mathbf k)\) in the ground \(S^z=1/2\) sector. For \(\mathbf k\neq0\), its disconnected expectation value vanishes.

| \(J'/J\) | \(M_{0,zz}(0,1/3)\) | \(M_{1,zz}(0,1/3)/J\) | \(\bar\omega_{zz}/J\) |
|---:|---:|---:|---:|
| 0.2 | 0.331770553 | 0.006921876 | **0.020863443** |
| 0.5 | 0.257735696 | 0.097971933 | **0.380125589** |
| 1.0 | 0.208208415 | 0.373480824 | **1.793783528** |

For \(J'/J=0.2\), the ground-state cell momentum is \(\mathbf k_0=(1/3,2/3)\). Under \(O_h^z(0,1/3)\), the target is

\[
\boxed{(S^z_{\rm tot},\mathbf k_{\rm target})=(1/2,(1/3,0))}.
\]

This is accessible without changing \(S^z\) sectors. Since \(S_g=1/2\), rank-one SU(2) selection allows \(S_f=1/2\) or \(3/2\). The longitudinal hub moment establishes an excitation at or below \(0.020863443J\) **in the target sector**, without running dynamical ED.

At \(J'/J=0.5\), \(S_g=5/2\), and longitudinal spin transitions may reach \(S_f=3/2,5/2,7/2\). The \(\Gamma\)-point *fixed-\(S^z\) connected* hub channel has

\[
M_{0,zz}^{\rm conn}(\Gamma)=0.145914428,\quad M_{1,zz}(\Gamma)=0.097971933J,
\quad \bar\omega_{zz}=0.671434173J.
\]

For the same coupling, the generic \((0,1/3)\) hub channel's centroid is only \(0.380125589J\). Hence it is worth computing finite momentum even though the **raw spin-vector static** structure factor peaks at \(\Gamma\).

## 4. Zero-energy ground-multiplet projection at \(\Gamma\)

For a unique orbital ground sector with total spin \(S\), the Wigner–Eckart theorem implies

\[
P_0\,\mathbf H_{\rm tot}\,P_0=\alpha_h\mathbf S_{\rm tot},\qquad
\alpha_h=\frac{\langle\mathbf H_{\rm tot}\cdot\mathbf S_{\rm tot}\rangle}{S(S+1)}.
\]

The spin-vector SU(2)-averaged \(\Gamma\) elastic weight for \(\mathbf h=\mathbf H_{\rm tot}/\sqrt{9}\) is

\[
M_{0,hh}^{\rm elastic}(\Gamma)=\frac{\alpha_h^2 S(S+1)}9.
\]

At \(J'/J=0.5\), this gives

- raw \(S_h(\Gamma)=2.537259562\);
- ground-spin multiplet elastic weight **2.175178034** (85.7% of the raw spin-vector structure factor);
- nonzero-energy three-component weight **0.362081528**;
- its spectral centroid **\(0.724825710J\)**.

This large spin-vector elastic weight is predominantly in *transverse transitions between magnetic quantum numbers* when a fixed \(S^z=1/2\) representative is used. **Do not subtract 85.7% from the longitudinal fixed-\(S^z\) dynamic seed:** its appropriate connected zeroth moment is instead \(0.145914428\), as above. At \(J'/J=1\), additional degenerate primitive \(K/K'\) ground sectors require explicit \(P_0\) projection, not only the one-multiplet formula.

## 5. Stage 3 implementation priorities and acceptance

1. **First:** \(J'/J=0.2\), \(\mathbf k=(0,1/3)\), \(h^z\) initial seed. Starting \(\mathbf k_0=(1/3,2/3)\), target \((1/3,0)\); \(S^z=1/2\) unchanged. Reproduce \(M_{0,zz}=0.331770553\), \(M_{1,zz}=0.006921876J\), and a spectral centroid of \(0.020863443J\). This guarantees an accessible low-energy mode and is more informative than a blind scan.
2. **Second:** \(J'/J=0.5\), both \(\mathbf k=(0,1/3)\) and \(\Gamma\), \(n^z,m^z,h^z\) block. At \(\Gamma\), explicitly project out \(P_0\) and enforce the positive-frequency conserved null combination \(m+h\). Check \(M_{1,zz}^{hh}(\Gamma)=0.097971933J\) and \(M_{0,zz}^{hh,\rm conn}(\Gamma)=0.145914428\).
3. **Third:** \(J'/J=1\) as implementation benchmark. At cell \(\Gamma\), distinguish physical \(K/K'\) from primitive \(\Gamma\) by the basis phases, and resolve degenerate ground states.
4. For every seed \(b=O^z(\mathbf k)|\psi_0\rangle\), independently compute \(\langle b|b\rangle\), \(\langle b|(H-E_0)|b\rangle\), and the corresponding static/double-commutator moments **before** producing a broadened spectrum. For complex \(\mathbf k\), check \(\mathbf k\) and \(-\mathbf k\) and the partner-state symmetry average where needed.
5. Record the discrete poles and their residues; broadening \(\eta\) is display-only and not a physical lifetime. Do not infer a rim Goldstone velocity from a single 3×3 torus.

### Independent algebra check

A separate 12-site \(L=2\), fixed-\(S^z=0\) sparse-ED test verifies the longitudinal f-sum directly against \(\langle b|(H-E_0)|b\rangle\):

- \(J'/J=0.2\), \(S=0\): both methods give **0.005239235066** (difference \(9.3\times10^{-17}\)).
- \(J'/J=0.5\), \(S=1\): both methods give **0.115720827095** (difference \(2.5\times10^{-16}\)). The \(S=1\) example specifically tests the correct \(G_{\rm dot}-G_{zz}\) formula instead of incorrectly assuming \(G_{zz}=G_{\rm dot}/3\).

## Files

- `stage3_first_moments.csv`: three-component moments, elastic projections where justified, and centroids.
- `stage3_first_moments_longitudinal.csv`: fixed-\(S^z\) longitudinal moments and connected static weights.
- `stage3_moment_analysis.py`: full reproduction from the Codex ZIP.
- `stage3_small_cluster_crosscheck.py`: independent 12-spin direct-H check.

**Scientific limitation:** These are finite-cluster spectral energy bounds. They establish neither a thermodynamic spin-liquid excitation nor Goldstone damping; Stage 3 is needed for resolved poles and rim–hub spectral-weight transfer.
