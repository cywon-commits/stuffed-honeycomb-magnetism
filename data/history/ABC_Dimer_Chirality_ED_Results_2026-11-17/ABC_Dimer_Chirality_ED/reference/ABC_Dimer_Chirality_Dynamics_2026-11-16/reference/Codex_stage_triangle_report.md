# 27-spin ABC-triangle correlations and the isotropic endpoint

## Established: ground sectors and exact benchmarks

All ground searches use Nup=14, Sz=1/2 (20,058,300 states before symmetry reduction). Every half-integer spin multiplet is represented. Cold generic starts avoid locking to an endpoint spin. The three cell-momentum stars have multiplicities 1,6,2, with compatible C6/C2/C3 reductions. Eight representative symmetry families cover their characters; conjugacy and the verified corner reflection cover the omitted equivalent characters. Uncomputed families are explicitly bounded, never assigned Ritz eigenvalues.

| Jprime/J | Ground energy/J | Cell momentum | S | C6 character | reflection | direct H residual/J |
|---:|---:|---|---:|---:|---:|---:|
| 0.6 | -11.686276944520 | Gamma | 1.5 | 0 | +1 | 6.62e-12 |
| 0.7 | -12.384408843833 | Gamma | 0.5 | 3 | -1 | 4.00e-12 |
| 0.8 | -13.290782890307 | Gamma | 0.5 | 3 | -1 | 9.37e-14 |
| 0.9 | -14.205410071922 | Gamma | 0.5 | 3 | -1 | 2.77e-14 |

Together with the verified 0.5 S=5/2 state, these samples establish ground-sector changes somewhere in [0.5,0.6] and [0.6,0.7]. At 0.7 the competing Gamma, character-0, S=3/2 minimum lies 0.001555624145J higher. No additional crossing scan was performed: narrow intervening sectors are not excluded. These are finite-cluster level changes, not thermodynamic transitions. Dimer and scalar-chirality probes are spin scalars and cannot directly connect crossing states of different total spin.

All 13 supplied ground/competitor reference rows reproduce with zero difference at the stored precision; input snapshot hashes match. Explicit 8x8 algebra and independent full-bit-basis 12-spin and 27-spin/two-up tests validate signs, Fourier normalization, overlapping operator order and RDMs. All 54 elementary faces have one AB and two hub bonds. Ordered face pairs split into 54 identical, 162 shared-edge, 486 shared-vertex and 2,214 disjoint pairs.

## Established: local doublet density and basis dependence

The local basis is |a,b,c> with computational index a+2b+4c and 1=up. We use orthonormal doublet flavors: AB-singlet times the C spin, and AB-triplet coupled to C to total spin 1/2, with consistent lowering phases. The three cyclic singlet-pairing vectors have off-diagonal overlaps -1/2, Gram eigenvalues 0,1.5,1.5, and sum to zero. Their singlet probabilities are not exclusive dimer occupations.

The complete fixed-M 8x8 reduced state is retained before averaging. The spin-multiplet twirl is reconstructed exactly from the three bonds and scalar chirality; its 2x2 spin-traced flavor matrix preserves chirality coherences. The restored ensemble includes spin-multiplet, point-group and conjugate-partner averaging, plus primitive-translation cosets at isotropy.

| j | doublet weight | quartet weight | conditional flavor purity, restored | physical-K magnetic S-dot |
|---:|---:|---:|---:|---:|
| 0.2 | 0.769963253 | 0.230036747 | 0.686867036 | 2.244681414 |
| 0.5 | 0.831990149 | 0.168009851 | 0.571242243 | 3.259990500 |
| 0.6 | 0.856833570 | 0.143166430 | 0.527663217 | 3.142081228 |
| 0.7 | 0.872483798 | 0.127516202 | 0.503997030 | 2.970813370 |
| 0.8 | 0.873161693 | 0.126838307 | 0.503344266 | 2.976779954 |
| 0.9 | 0.873412381 | 0.126587619 | 0.502940441 | 2.966798291 |
| 1 | 0.873480824 | 0.126519176 | 0.500000000 | 2.948982074 |

At isotropy the restored flavor state is I/2, not a pure coherent three-pairing state. The selected physical-K state has CCW face chirality alternating between approximately +/-0.01849477979; its conjugate partner reverses it. The restored mean is zero. Conditional flavor purity is 0.5011955295 in that selected state, versus 0.5 after restoration.

A different eigenbasis inside the same ground pair makes the endpoint caution concrete. The character-3 cell-Gamma combination has AB=-0.218284262864, BC=CA=-0.170968486810, zero chirality to numerical precision, and conditional flavor purity 0.5026082707. It has the same energy and doublet weight. Its bond values continue the anisotropic 0.9 trend, whereas the primitive-K or fully restored state has all bonds equal to -0.186740412162. Thus endpoint changes in local bond favoritism depend on restoring the enlarged isotropic symmetry and on the degenerate ground-state basis; they are not by themselves a new phase boundary.

The 54 faces overlap; no independent tensor-product lattice of 54 trimer flavors is derived. At isotropy the averaged doublet weight is already fixed by the energy, P=1/2-2E/81. The exact local identities remain chi^2=(3/16)P_doublet and the symmetrized D norm=(9/8)P_doublet. Their growth adds no independent resonance evidence. At isotropy they give 0.1637776546 and 0.9826659274. Raw selected-state RDMs, chirality, orthonormal flavor matrices and the alternative endpoint basis are saved separately.

## Established: genuine intertriangle fluctuations

All connected matrices are obtained as exact operator-seed Gram matrices, retaining complex off-diagonal elements. This evaluates the ordered products even for overlapping faces, without Wick factorization or assuming commutation. Selected disjoint chirality products are real to numerical precision, while selected overlapping products retain their allowed complex parts. Full covariances are PSD and satisfy Cauchy-Schwarz; category-masked matrices are not separately physical covariance matrices.

At isotropy the following decomposition projects each geometric category onto the leading full physical mode:

| field | representative primitive Q | leading equal-time eigenvalue | identical | shared edge | shared vertex | disjoint |
|---|---|---:|---:|---:|---:|---:|
| D | [0.8888888888888888, 0.7777777777777777] | 2.7576254235 | 0.9826659274 | 1.2364438984 | 0.2213195471 | 0.3171960506 |
| chi | [0.0, 0.0] | 0.3681687648 | 0.1637776546 | 0.0656706545 | 0.0399892745 | 0.0987311812 |

The dimer maximum belongs to the smallest nonzero primitive-momentum shell (one representative is equivalent to (-1/9,-2/9)), whereas chirality peaks at primitive Gamma. Disjoint pairs contribute about 11.5% of the leading dimer eigenvalue and 26.8% of the leading chirality eigenvalue. This confirms independent four-/six-spin spatial information beyond the local projector identities, but also shows why shared-site contributions must not be called long-range resonance.

The magnetic structure factor uses S(Q)=sum_ij exp[-2pi i Q.(r_j-r_i)] <Si.Sj>/27. Its maximum is 2.948982073595 at physical K=(1/3,2/3) and Kprime=(2/3,1/3), consistent with strong 120-degree correlations. Dimer/chirality and magnetic magnitudes have different operator powers and normalizations and are not directly interchangeable measures of order. Primitive K folds to cell Gamma; physical Gamma is distinct. No sector is labeled M.

Six face channels are unfolded into two geometric up/down-face channels at each of 27 allowed primitive momenta. A translation-covariant geometric anchor removes the color-dependent dimer phase. Unitarity, local permutation phases, and vanishing cross-fold blocks independently check this unfolding.

## What the ground-manifold and ensemble checks change

Total ensemble-connected covariance is decomposed into the average of within-selected-state connected covariances plus covariance of image-state means. This split depends on the chosen ensemble decomposition, so it is not a basis-independent quantum/classical classification. Both components and the selected source data are supplied.

A separate static projection onto the verified K/Kprime ground pair identifies zero-energy orbital contributions. They occupy subleading dimer channels at physical +/-K (eigenvalue about 0.1208943) and a subleading chirality channel at Gamma (about 0.0184711). The leading maxima quoted above are unchanged by removing this known ground-pair contribution. The remaining positive covariance is not interpreted as a resolved dynamical resonance spectrum.

## Suggested, and not yet established

**Suggested:** the doublet-flavor language is useful for organizing frustrated local states, and the nonzero disjoint covariances establish finite-cluster valence-bond and chirality fluctuations beyond trivial bond equivalence.

**Not established:** a distinct coherent trimer-resonance phase, chiral long-range order, spontaneous thermodynamic symmetry breaking, a tunneling frequency, or a lifetime. Local flavor mixedness does not rule out many-body coherence, but does not prove it. Finite-size scaling or a carefully ground-manifold-deflated scalar dynamical calculation would be needed. Optional dynamical Task D was not started.

## Validation and measured resources

Maximum 27-spin local-identity error: 3.031e-14; maximum RDM trace/Hermiticity/projector error: 4.541e-14; minimum full connected-covariance eigenvalue: -1.122e-14. Cauchy-Schwarz violations are zero at recorded precision. Every new ground passes direct H and S-squared residuals below 1e-8; actual values are in the tables. The reused isotropic checkpoint retains H residual 8.30e-9 and S-squared residual 1.72e-8, documented rather than relabeled as a tighter solve.

Measured new symmetry-block solves: 14; iterative H calls: 1742; summed solve walls: 32.09 minutes; summed seven-state observable walls: 15.12 minutes. These times overlap and exclude separately uninstrumented validation/projection/plotting, so are not elapsed end-to-end time. Hardware quota is two CPUs and 8 GiB RAM. Detailed costs, conditional exclusion evidence and raw numerical precision are preserved.
