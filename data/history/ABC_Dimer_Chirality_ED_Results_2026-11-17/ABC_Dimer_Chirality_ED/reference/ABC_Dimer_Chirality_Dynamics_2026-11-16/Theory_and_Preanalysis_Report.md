# ABC local-flavor theory and 27-spin dynamical target preanalysis

**Date:** 2026-11-16. **Completed here:** new exact 8×8 local analytic/numeric dimer/chirality dynamics; independent reprocessing of archived 27-spin static Gram matrices including known K/K' ground-pair subtraction; exact symmetry and commutator tests. **Not completed here:** converged 27-spin dynamical scalar poles/line shapes and thermodynamic mode dispersion.

## 1. Symmetry and spin selection

Each elementary ABC triangle has three spin-1/2 objects. Under global SU(2), its Hilbert space decomposes as (1/2⊕1/2)⊕3/2. The two copies of S=1/2 are the **two-dimensional flavor subspace**; the three singlet-pairing configurations are linearly dependent. Under cyclic C3, circular flavor states have phases ω and ω², combining as E under C3v. The chirality χ is an A2 pseudoscalar, D=B_AB+ωB_BC+ω² B_CA has C3 phase ω² and forms E with D†. D_re and D_im are both time-reversal even (they are real combinations of bond operators), whereas χ is time-reversal odd.

D and χ commute with the total spin SU(2), hence *exactly preserve the total S of a finite-size eigenstate*. In the isotropic 27-spin ground S=1/2, all scalar modes have S_f=1/2. This is markedly different from magnetic rank-one spin spectroscopy. A point-group selection Γ_f⊂Γ_O⊗Γ_g additionally filters states when the little group permits it. For j<1 the local AB anisotropy reduces the local mirror channel: D_im and χ are odd under swapping A and B, and connect the AB-singlet/triplet doublet flavors; D_re is even, diagonal in those flavors.

## 2. Exact doublet pseudospin algebra

In orthonormal doublet flavor basis |0>=AB singlet × C spin and |1>=AB triplet combined with C to total spin 1/2, at fixed M=+1/2,

```
H_doublet=-(J/4+J'/2)I -(J-J')τ_z/2
D= -3τ_z/4 + (3i/4)τ_x
χ= +(sqrt(3)/4)τ_y
```

The singlet-flavor |0> has E0=-3J/4; the triplet-flavor |1> has E1=J/4-J'. For J'<J, Δ=J-J' exactly. The scalar quartet (S=3/2) is forbidden in these transition matrix elements. Isolated-trimer inelastic residues are Wχ=3/16, WD=9/16, and a D elastic weight 9/16. At J'=J the local doublets are exactly degenerate, so the isolated triangle does **not** define a nonzero oscillation frequency. In a 27-spin lattice, overlapping ABC triangles strongly invalidate independent pseudospin descriptions and lift/disperse these collective freedoms.

Exact local operator identities, independently verified by explicit 8×8 matrices:

```
P²=P
χ²=(3/16)P
(D†D + DD†)/2=(9/8)P
D†D - DD†=3sqrt(3) χ
[χ,D]=-(sqrt(3)/2)D
D²=0
```

In a selected chirality state, the integrated D-vs-D† spectral asymmetry (including ground-manifold ω=0 contribution) measures χ exactly. From the previously reported |<χ_CCW>|≈0.01849477979 at the isotropic K/K' pair, |ΔW_D|≈0.09610169481. It vanishes after the symmetric K/K' average, not necessarily in each selected source. This is a *theoretical prediction from the pre-existing static data*, not a newly resolved dynamical spectrum.

## 3. Recomputed 27spin static targets (j=1)

The original triangle Codex stage provided 54-face complex Gram matrices, physically unfolded 27 momenta, and a projector onto the known orbital K/K' ground pair. We independently reconstructed the two physical face channels per Q, their spectra, and the *positive-energy* M0 by subtracting that ground pair. Full data: `results/27spin_static_mode_targets.csv`, reproducible via `code/static_mode_summary.py`.

| Channel at physical primitive Q | largest M0 positive energy | Known ground-pair weight on same eigenmode |
|---|---:|---:|
| D at sixfold star Q=(8/9,7/9) | 2.75762542347 | ~0 |
| χ at Q=(0,0) (physical Γ) | 0.36816876477 | 0 |
| D at physical K=(1/3,2/3) | 0.40807821491 | ~0 *for leading mode* |
| χ at physical K | 0.09702886380 | ~0 *for leading mode* |

For an **orthogonal subleading** D face channel at physical K, the known ground pair has M0=0.120894263834; for subleading χ Γ, M0=0.0184710714683. These must be removed as zero-energy contributions, not mistaken for scalar soft modes. Original report independently showed max dimer and chirality static covariance at small finite physical q and physical Γ respectively. Neither quantity is a finite-frequency peak position.

At r=1 the physical triangular K folds to the **three-site cell Γ**. The former is never interchangeable with physical Γ. For r<1, 27 primitive translations cease to be symmetries; retain six cell-basis face channels and only nine cell momenta.

## 4. Dynamical questions that can now be posed precisely

1. Does χ(physical Γ), a spin scalar/A2-type mode, have a finite-frequency isolated low-energy pole, or mostly broad high-energy weight? Its leading static M0 cannot be explained by the known K/K' ground transition.
2. Do dimer modes at physical Q≈(8/9,7/9) concentrate their unusually large inelastic M0=2.75763 in a small number of converged S=1/2 excitations?
3. Does a nonzero complex cross-spectral residue link D and χ at generic momenta? At Γ reflection symmetry can prohibit some cross channels. Do not enforce zero at a generic Q without checking its little group.
4. How do these scalar modes evolve from r=.9 to r=1, while the local isolated-trimer flavor gap J-J' closes? A peak energy moving to zero is *not* predicted for the coupled lattice merely by the isolated-trimer formula.
5. Does a selected K source exhibit D/D† integrated circular asymmetry |ΔW|≈0.09610, and does this disappear in the fully symmetry-restored ensemble as required?

## 5. Scope of validation

The supplied standalone local script recalculates every 8×8 operator matrix and inelastic spectral moment M0,M1,M2, with max local algebra error ≤~7e-16. Static 27spin script verifies momentum unfolding and known ground-pair deflation from existing exact spin correlators. The **only** 27-spin results represented as verified here are the supplied static covariance matrices and derived M0; no 27-spin Hamiltonian spectrum is claimed as newly computed. The 27-spin full scalar dynamical response is the next Codex task, with a concrete reproducible execution specification alongside this report.
