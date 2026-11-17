# Stuffed honeycomb 27-spin ED — observables and selection rules (v1.0)
Date: 2026-11-15. Model: J rim-rim, J' hub-rim, PBC 3x3 cells, 27 spins, spin 1/2.

## Evidence and scope

Grounding: `tri27_verified_results.zip`, specifically `tri27_verified/ed.py`, `report.md`, and `momentum_verification.json` (Codex output). Here formulas for observables and selection rules are analytical proposals, **not results of completed dynamical ED**.

## Lattice convention

Cell R=(x,y), x,y in Z_3. Site i=3*(3*x+y)+s, s=0(A),1(B),2(H).
Triangular primitive e1=(1,0), e2=(1/2,sqrt(3)/2) in Cartesian units.
Cell vectors a1=e1+e2, a2=2e1-e2; offsets r_A=0, r_B=e1, r_H=e2.
Rim edges A_R—B_R, A_R—B_{R-a1}, A_R—B_{R-a2}.
Six rim vertices of hexagon around H_R:
A_R, B_R, A_{R+a1}, B_{R+a1-a2}, A_{R+a1-a2}, B_{R-a2}.
N_c=9, N_r=18, N_h=9. Hub–rim Φ_R is sum of these six rim spin operators.

At J'!=J, only cell translations are symmetries; k=(n1/3,n2/3), n1,n2=0,1,2 (reciprocal coordinates conjugate to a1,a2). At J'=J, primitive translations also exist; physical Q=(q1,q2) folds to k=(q1+q2,2q1-q2) mod 1. Physical triangular K=(1/3,2/3) folds to cell Γ=(0,0), but is **not** the same primitive momentum as physical Γ. No exact triangular M is allowed on this odd-order torus.

## Core operators

For α∈{A,B,H}, a∈{x,y,z}:

O_α^a(k) = (1/sqrt(N_c)) sum_R exp[-2π i (k1 x+k2 y)] S^a_{R,α}.

Here T_R = S_{R,H}. Convention: if translation U_δ acts as U_δ S_{R,α} U_δ^{-1}=S_{R+δ,α}, then U_δ O_α(k) U_δ^{-1}=exp(+2π i k·δ) O_α(k). Thus a state with cell momentum k0 is sent to k0+k (mod integer reciprocal vectors). Confirm this against the projected-basis phases used by the solver.

Define n=O_A-O_B (rim staggered), m=O_A+O_B (rim uniform), h=O_H.

Define Φ(k)=F_A(k)O_A(k)+F_B(k)O_B(k), with
 F_A=1+exp(2πi k1)+exp[2πi(k1-k2)],
 F_B=1+exp[2πi(k1-k2)]+exp(-2πi k2).
Hence Φ(Γ)=3m(Γ) and Φ(k)=((F_A+F_B)/2)m(k)+((F_A-F_B)/2)n(k).
For physical long-wave continuum fields, include physical basis offsets when expanding q→0; the staggered *physical* form factor starts at q³, not at the naive linear term in F_A-F_B due to basis phases.

At isotropy, define ω3=exp(2πi/3). Physical K probe at cell Γ:
 O_phys(K)= (1/sqrt(27)) sum_R [S_A,R + ω3² S_B,R + ω3 S_H,R]
(componentwise in spin), since exp[-iK·a1]=exp[-iK·a2]=1 and offsets are 0,e1,e2.
This probes the three-sublattice 120° pattern. In the J'!=J model this is a diagnostic Fourier combination, not a new conserved primitive momentum.

## Static observables

Use the Hermitian 3×3 matrix
 C_{αβ}(k)=sum_{a=x,y,z} <δO_α^a(k)† δO_β^a(k)>.
Check C=C† and eigenvalues≥0. At finite fixed Sz=1/2, subtract possible expectation values, especially at k=Γ.

Rim Néel finite-size amplitude:
 m_N² = < [sum_R(S_A,R-S_B,R)]² > / N_r²
      = < n(Γ)†·n(Γ) > * N_c / N_r² (with disconnected term treated separately).
Rim NN correlations, hub-hub structure factor S_h(k)=sum_a <h_a(k)† h_a(k)>, and rim–hub cross correlators C_{Ah}, C_{Bh} are required.
Ferrimagnetism diagnostic:
 m_F²=<S_tot²>/27²=S_tot(S_tot+1)/27².
Do NOT infer net magnetization from hub S_h(Γ) alone: at J'=J, a 120° triangular pattern has a *uniform hub sublattice* while total moment cancels.
120° diagnostic: <O_phys(K)†·O_phys(K)> and/or projection of C_{ABC}(Γ) onto v=(1,ω3²,ω3); evaluate K and K' components.
Finite-size indicators ≠ thermodynamic phase transitions; check multiple geometries/sizes if possible.

## Dynamical response

Define zero-temperature positive-frequency spectral matrix (all components as needed):
 A_{αβ}^a(k,Ω) = sum_{n:E_n>E_0} <0|δO_α^{a}(k)†|n><n|δO_β^{a}(k)|0> δ[Ω-(E_n-E_0)].
A is positive semidefinite at each Ω in the sense of a matrix-valued measure. Corresponding retarded response is obtained by the Lehmann representation. Use a consistent susceptibility convention in which χ''(ω>0)≥0 for diagonal Hermitian channels; include negative-frequency partners when building χ^R.
Compute an independent dynamic spectral matrix for (n,m,h) by linear combinations of A/B/H, retaining cross terms and phases. Φ is derived from A/B and needs no additional independent seed.

Important q=Γ exact SU(2) Ward check:
 S_tot^a = sqrt(N_c)[m^a(Γ)+h^a(Γ)] and [H,S_tot^a]=0.
For all E_n>E_0, <n|m^a(Γ)|0>=-<n|h^a(Γ)|0>.
Therefore at finite Ω>0, in the (m,h) subspace:
 A_{mm}=A_{hh}=-A_{mh}=-A_{hm}
(under matching component/normalization conventions, up to numerical error).
The conserved total-spin combination has no **inelastic** weight at Γ; possible Ω=0 multiplet transitions require separate handling.

Hub-centered inversion exchanges A↔B, fixes hub sites up to cell permutation, and preserves all bonds at any J,J'. Thus at Γ n is inversion-odd and h/m inversion-even. In an inversion-invariant eigenstate or symmetry-averaged ensemble, C_nh(Γ)=χ_nh^R(Γ,ω)=0 by selection rules. Do not interpret this zero as no hub–Goldstone coupling. At J'=J a single primitive-K eigenstate is mapped to K' by inversion; averaging/checking the degenerate symmetry pair is needed to apply the rule.

## Spin/momentum sector selection in 27-site ED

N_up=14 corresponds to Sz=+1/2 and dimension C(27,14)=20,058,300. For J'!=J, every cell-momentum sector has dimension 2,228,700 (all up14 translation orbits length 9). The ground cell momentum k0 must be computed independently at each J'/J, especially across crossings.

O^z_α(k): ΔSz=0, goes from (N_up=14,k0) to (N_up=14,k0+k). O^-_α(k): ΔSz=-1, goes to (N_up=13,k0+k), with full-sector dimension 20,058,300 and cell-sector dimension 2,228,700. O^+_α(k): ΔSz=+1, goes to (N_up=15,k0+k), with full-sector dimension C(27,15)=17,383,860. At N_up=15 some orbits have length 3, and cell sector dimensions are 1,931,652 at Γ and 1,931,526 at each non-Γ k.

Spin-vector operators obey rank-one SU(2) selection: |S_f-S_g|≤1 and S_f+S_g≥1; with odd 27 sites, if S_g=1/2 only S_f=1/2 or 3/2 can couple. **S_g=1/2 itself must be checked**, not assumed from Sz=1/2 alone. For J'=J the verified ground state has primitive K or K', both folded to cell Γ; a physical momentum-transfer Q adds to the initial **primitive** momentum, while a 3-cell operator generally mixes three unfolded momenta if it lacks proper basis phases.

## Spectral sum rules and numerical acceptance

(1) Directly compare H action + symmetry/projected seed against a small L=2 fully enumerated reference for each spin operator (z,+,-) and both complex momentum sectors.
(2) Zeroth moment (including all E_n=E0 components where relevant): ∫dΩ A_{OO}(k,Ω)=<0|δO†(k)δO(k)|0>.
(3) First moment for an operator with positive-frequency spectral measure: ∫dΩ Ω A_{OO}(k,Ω)=<0|δO†(k)(H-E0)δO(k)|0>.
(4) Γ conservation identity above for every inelastic pole/energy bin.
(5) SU(2), inversion, time reversal and K/K' partner checks, where valid for the chosen state or ensemble.
(6) Static and dynamical correlators should be evaluated at the *same* coupling, state, geometry, and quantum-number conventions. Recompute ground-state momentum after every J' change.
(7) Continued-fraction broadening η and Lanczos step count must be converged; discrete ED δ-peaks are NOT a thermodynamic linewidth. Record η, moments, residuals, sum-rule error and available sector labels.

## Minimum next-run specification for Codex

Stage 1: Extend the verified code with reusable operators A(k),B(k),H(k), n(k),m(k),Φ(k), physical-K 120° probe; independently validate momentum phases and projection on L=2.
Stage 2: From a verified ground vector at J'=J (and later selected J'/J where primitive translation is NOT used), compute C_{αβ}(k), m_N², S_h(k), m_F² or <S_tot²>, and 120° structure factors. For ground spin multiplets, average over symmetries or report explicitly which Sz/partner state is used.
Stage 3: Implement the longitudinal spectral matrix A_{αβ}^{zz}(k,Ω) at a few k, using block Lanczos seeded with A/B/H spin-operator vectors. It remains in N_up=14. Test moments and conserved-spin null channel at Γ. Only then extend to transverse ± and their N_up=13/15 sectors.
Stage 4: Compare the full rim/hub spectral matrix to an isolated 18-rim-spin reference at J'=0. Refrain from calling an inverse-susceptibility difference a uniquely defined irreducible hub self-energy; specify projection and reference convention first.

Priority coupling points for *separate* tests (not full scanning): J'/J=1.0 (existing benchmark), an intermediate point such as 0.5, and weak coupling 0.2, subject to runtime. Especially near J'→0, expect near-degenerate hub multiplets and stringent convergence requirements.

## Scope warnings

- K at triangular primitive coordinates is not code-cell k=(1/3,1/3).
- The 27-site torus has no exact M point.
- 3-cell k=Γ at J'=J contains physical primitive Γ and K/K'; use basis-phase resolved physical probes to distinguish them.
- Finite-size static squared moment is not spontaneous order; finite-size level poles are not finite lifetime.
- A 3×3 torus alone cannot establish order destruction or true damping.
