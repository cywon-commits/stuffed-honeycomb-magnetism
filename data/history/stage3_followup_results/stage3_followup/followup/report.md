# Stage 3 follow-up: weak-coupling spectral matrix

All data use N=27, Nup=14, M=Sz=1/2, J=1, source cell momentum (1/3,2/3), probe (0,1/3), and forward target (1/3,0). The reverse probe around the same source targets (1/3,1/3); the correct conjugate source at (2/3,1/3) with negative probe instead targets (2/3,0). These operations are distinguished in partner_source_validation.json. Fourier normalization is 1/sqrt(9)=1/3, with negative Fourier exponent and no extra 2pi spectral prefactor.

## Additional converged poles at j=0.2375

| omega/J | S | F | hub weight | n/h intensity ratio | H residual/J |
|---:|---:|---:|---:|---:|---:|
| 0.001831371839 | 1.5 | +1 | 0.0739343634089 | 0.00960709 | 6.42e-13 |
| 0.021040229692 | 0.5 | +1 | 0.0359241424688 | 0.00944763 | 1.05e-11 |
| 0.022681526628 | 1.5 | -1 | 0.0817954672824 | 0.00980586 | 1.31e-11 |
| 0.022730224050 | 1.5 | +1 | 0.0442031643238 | 0.0095336 | 5.60e-13 |
| 0.041099346210 | 0.5 | -1 | 1.3552088776e-06 | 0.0820887 | 1.18e-11 |

The original two verified S=1/2 states are explicitly projected out of the initial seed and shifted by rank-one deflation in their respective reflection blocks. Subsequently verified states are added to the deflation set. Shift 64J exceeds the full bond spectral-width bound. ARPACK uses ncv=64 and tol=2e-12; every reported eigenpair passes a direct original-H residual, Casimir and orthogonality check. Both compatible reflection irreps carry seed weight and are searched. Probe-generated seeds include the allowed S=1/2 and 3/2 subspaces. This is not a proof of complete enumeration of all levels in either block.

## Cumulative hub weight and explicit unresolved bounds

| j | cutoff/J | verified weight | verified fraction | lower bound on total below cutoff | upper bound |
|---:|---:|---:|---:|---:|---:|
| 0.2375 | 0.01 | 0.1652994890 | 49.93350% | 0.1652994890 | 0.1654864854 |
| 0.2375 | 0.04 | 0.3272222631 | 98.84696% | 0.3272222631 | 0.3274159358 |
| 0.2375 | 0.05 | 0.3272236183 | 98.84737% | 0.3272236183 | 0.3274195947 |
| 0.2375 | 0.1 | 0.3272236183 | 98.84737% | 0.3272236183 | 0.3274317408 |
| 0.2 | 0.01 | 0.1657091809 | 49.94692% | 0.1657091809 | 0.1658087381 |
| 0.2 | 0.04 | 0.3292012341 | 99.22557% | 0.3292012341 | 0.3293043302 |
| 0.2 | 0.05 | 0.3292012341 | 99.22557% | 0.3292012341 | 0.3293055512 |
| 0.2 | 0.1 | 0.3292012341 | 99.22557% | 0.3292012341 | 0.3293119894 |

Verified weights are independent pole residues. Bounds are inferred from positivity, not extra measured states. If unresolved moments are r0,r1,r2, the first-moment bound is W_unresolved(<T)>=max(0,r0-r1/T). The upper bound also uses bounded-support moment inequalities and the one-sided Chebyshev inequality when the remaining centroid exceeds the cutoff. A stronger bounded-band lower bound is max(0,r0-[(L+T)r1-r2]/(LT)), where L=sum_bonds(J_b)/4-E0 bounds every excitation from above. The polynomial x(L+T-x)/(LT) majorizes the tail indicator on [0,L]. Bounds are conditioned on the verified numerical ground state and operator moments, with a displayed 1e-8 roundoff margin; they are not formal interval-arithmetic eigenvalue certificates.

Additional verified hub weight below 0.1J is 0.235858492693; the original moment requirement was at least 0.1391992364980938. Remaining matrix moments, PSD checks, and partial sums are in spectral_remainders.json. No fictitious pole energy is assigned to unresolved weight.

For an unresolved centroid mu>T, the one-sided Chebyshev (Cantelli) inequality additionally gives W_unresolved(<T)<=r0*variance/[variance+(mu-T)^2], where variance=r2/r0-mu^2. This makes the residual low-window weight much smaller than the total unresolved mass.

| j | verified fraction of total hub M0 below 0.1J | upper bound on total fraction below 0.1J | guaranteed recovered fraction of window weight |
|---:|---:|---:|---:|
| 0.2375 | 98.847367% | 98.910239% | 99.936435% |
| 0.2 | 99.225574% | 99.258960% | 99.966365% |

This bounds the missing hub weight in the requested window; it is not a proof that every individual eigenstate in that window is listed.

## Controlled coupling comparison

All energies and full complex 3x3 residues are exported in verified_poles.csv and complex_residues.csv, including energies divided by Jprime squared. Branches are matched by spin/reflection labels and actual eigenvector overlaps, not sorted index. matched_branch_overlaps.json records the evidence; only matched available states are compared. The earlier two S=1/2 branch overlaps exceed 0.9986. Missing counterpart states remain unresolved.

| j | omega/J | S | F | hub weight | rim Neel weight | Re R_nh | Im R_nh |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.2 | 0.000098187481 | 0.5 | +1 | 0.0094431293 | 6.3549271e-05 | -0.00067087852 | -0.00038733189 |
| 0.2 | 0.000108524474 | 0.5 | -1 | 0.082303353 | 0.00053524636 | -0.0057479934 | -0.0033186055 |
| 0.2 | 0.001074787802 | 1.5 | +1 | 0.073962698 | 0.00048063267 | -0.0051634937 | -0.0029811445 |
| 0.2 | 0.014697166970 | 0.5 | +1 | 0.036294533 | 0.00023271639 | -0.0025168928 | -0.0014531288 |
| 0.2 | 0.015675507400 | 1.5 | -1 | 0.08229353 | 0.00054213945 | -0.0057845421 | -0.0033397069 |
| 0.2 | 0.015708364677 | 1.5 | +1 | 0.044903991 | 0.00029030947 | -0.0031268259 | -0.0018052738 |
| 0.2375 | 0.000204512224 | 0.5 | +1 | 0.0095479137 | 9.6561523e-05 | -0.00083154725 | -0.00048009403 |
| 0.2375 | 0.000235370037 | 0.5 | -1 | 0.081817212 | 0.00078656004 | -0.0069473457 | -0.0040110519 |
| 0.2375 | 0.001831371839 | 1.5 | +1 | 0.073934363 | 0.00071029421 | -0.0062758555 | -0.0036233669 |
| 0.2375 | 0.021040229692 | 0.5 | +1 | 0.035924142 | 0.00033939795 | -0.0030239767 | -0.0017458938 |
| 0.2375 | 0.022681526628 | 1.5 | -1 | 0.081795467 | 0.00080207463 | -0.0070145956 | -0.0040498787 |
| 0.2375 | 0.022730224050 | 1.5 | +1 | 0.044203164 | 0.0004214151 | -0.0037377682 | -0.0021580015 |
| 0.2375 | 0.041099346210 | 0.5 | -1 | 1.3552089e-06 | 1.1124734e-07 | -3.362626e-07 | -1.941413e-07 |

The bright newly resolved hub modes have finite complex rim cross residues, while their rim Neel intensities remain roughly one percent of hub intensities at 0.2375. A much weaker pole can have a larger relative rim fraction without carrying appreciable absolute rim weight. These are coupled-system spectra, not an independently extracted irreducible hub susceptibility. The matched very-soft S=1/2 energies do not collapse exactly in Jprime-squared units. The nearly unchanged integrated hub M0 and increased M1 do not by themselves specify individual pole shifts. The very small unresolved high-energy fraction nevertheless carries much of M1 and M2; its centroid is an energy mean, not an extra pole.

## Validation and status

M0=B-dagger B, M1=B-dagger(H-E0)B, and M2=[(H-E0)B]-dagger[(H-E0)B] are computed directly. The full 4x4 n,m,h,Phi moments are provided in four_channel_exact_moments.json; Phi is a known momentum-dependent linear combination of n and m, so this four-channel matrix is rank deficient by construction. Zeroth moments are independently checked against real-space Gzz and symmetrized first moments against the bond f-sum. Fixed-M zz quantities are used throughout, not scalar-dot/3 for arbitrary spin. The target has no ground-multiplet elastic contribution. The 0.2 original-H and total-spin source check is in matched_j0.2.json.

The referenced stage3_cross_moments_direct_27spin.csv was not supplied in the available attachments. We compare all supplied diagonal benchmarks and independently reconstruct the full cross-moment matrices instead. The supplied 0.2 hub M0,M1,M2 agree to better than 6e-15. Previously verified 0.2375 moments are reused. Small-complex-sector tests compare physically deflated eigenpairs with a dense eigenspectrum and validate the tail-bound polynomial. Full residue, remainder and block-Hankel positivity are checked to numerical precision.

| Data product | Status | Use |
|---|---|---|
| verified_poles.csv | Direct original-H and S-squared checked | Physical pole energies/residues |
| Prior 72-dimensional quadrature | Approximate nodes, not eigenpairs | Prior moment-complete visualization only; never a deflation input |
| spectral_remainders.json | Exact operator moments minus verified sums | Remaining spectral mass and positivity bounds |

CPU quota is two cores and RAM cap 8 GiB. Each new solve records H calls, elapsed solve seconds, ncv and process high-water RSS in its JSON; setup, shared-source preparation and matched-comparison timing are separately scoped. Per-process wall times overlap and must not be added as end-to-end elapsed time. Each iterative H callback has a 1200-second cutoff; final verification and checkpoint I/O may extend the recorded elapsed time beyond it; a timeout record is never a converged eigenvalue. No complete spectrum is claimed. Finite-size stick spectra and cumulative-weight plots use no broadening and imply no thermodynamic magnon damping, velocity, phase transition or self-energy.

Measured accepted new work: 9 eigenpairs, 1916 iterative H calls, 94.53 summed solve-wall minutes (overlapping jobs), peak process RSS 1588.8 MiB. These sums exclude source preparation, final validation and the stopped duplicate worker.

Unresolved moment ledger (centroids are means, not additional poles):

| j | unresolved fraction of hub M0 | unresolved fraction of hub M1 | unresolved centroid / J |
|---:|---:|---:|---:|
| 0.2 | 0.774426% | 62.175843% | 1.675048980 |
| 0.2375 | 1.152633% | 62.532282% | 1.650087740 |

Next targeted run, if needed: deflate the delivered verified vectors in both reflection sectors, use separate Casimir-projected S=1/2 and 3/2 starts, and stop according to the explicit residual weight bounds rather than a chosen count of Ritz values. Match newly found states across couplings using their delivered vectors and quantum labels.
