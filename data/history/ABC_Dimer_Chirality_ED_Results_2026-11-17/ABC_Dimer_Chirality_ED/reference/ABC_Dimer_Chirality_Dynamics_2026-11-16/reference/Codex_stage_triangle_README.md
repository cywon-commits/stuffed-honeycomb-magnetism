# 27-spin ABC-triangle ED results

Read report.md first. Task A, local-RDM Task B and static higher-point Task C are complete at the specified couplings. Optional dynamical Task D was not started. No thermodynamic phase or coherent tunneling frequency is asserted.

## Reproduction

From the archive root:

    OPENBLAS_NUM_THREADS=1 python triangle/run_triangle.py

This reproduces supplied reference data, derived tables, validation, figures and report from cached numerical arrays, without an eigensolver. Add --tests for independent explicit 12-spin and 27-spin/two-up operator tests, 8x8 identities and phase/unfolding checks.

    OPENBLAS_NUM_THREADS=1 python triangle/run_triangle.py --tests --compute

The expensive command reruns missing symmetry eigenvectors and all seven native-observable calculations, reusing provided ground states. Add --fresh-ground to move the delivered new eigenpair records/checkpoints to a backup directory and recompute them. Endpoints 0.2, 0.5 and 1 remain the original verified inputs. Only ground-state checkpoints are distributed; nonwinning symmetry-block vectors can be regenerated. Their measured energies/residuals and the explicit concavity exclusions are retained. H is matrix free throughout; dense matrices appear only for local 8x8 algebra and tiny validation spaces.

Python, NumPy, SciPy, Matplotlib and g++/OpenMP are required. Versions/hardware are in environment.json and requirements.txt. C++ shared libraries compile locally on demand and are not packaged as machine-specific binaries. The two-CPU, 8-GiB quota and measured timing scope appear in computational_costs.json.

## Data and normalization

* triangle_extended_ground.csv: 15 ground couplings (11 prior + 4 new), with measured/bounded sector evidence in triangle_sector_comparison.csv. Previous competing branches are separate.
* triangle_rdm.csv and triangle_rdm_matrices.npz: all 54 faces at seven couplings, selected fixed-M and symmetry-restored ensembles. Local basis index is a+2b+4c, with spin up=1; site order is ABC even when the physical CCW order differs. Translation-related RDMs are exactly equal and expanded from six orientation representatives.
* Orthonormal doublet flavor 0 is the AB singlet times C spin. Flavor 1 is the AB triplet coupled to C to S=1/2. For M=+1/2 these are (|up down up>-|down up up>)/sqrt(2), and sqrt(2/3)|up up down>-(|up down up>+|down up up>)/sqrt(6); negative-M partners use consistent lowering. The normalized 2x2 flavor density is the local doublet projection traced over its spin, divided by its weight. The three cyclic singlet-pairing vectors are nonorthogonal and linearly dependent.
* The exact local spin twirl is I8/8+(2/3)sum_bonds <B_b>B_b+(4/3)<chi_ABC>chi_ABC. Full selected RDMs are retained before this twirl, including chirality and doublet/quartet coherences.
* intertriangle_dimer_correlators.csv and intertriangle_chirality_correlators.csv contain 20,412 ordered pairs each, including identical, edge-sharing, vertex-sharing and disjoint classes. Complex selected/ensemble values and the image-ensemble mean covariance are retained. representative_disjoint_pairs_j1.csv identifies nearest and several maximal-distance pairs. Category-masked matrices need not be PSD; the complete connected matrices are PSD.
* D_mean fields use a CCW face anchored at A; D_ABC fields retain the uploaded fixed-label convention. For odd orientation mu, physical CCW order is A,C,B, so D_CCW=omega^2 D_ABC^dagger and chi_CCW=-chi_ABC. face_geometry.csv records all sites, lifted physical vertices and centroids.
* The cell Fourier operator is sum_R exp[-2pi i k.R] O_(R,mu)/sqrt(9). Momentum-resolved matrices, eigenvectors and decomposition are in momentum_resolved_matrices.json and raw_observables_j*.npz. Native output channel order is D(mu=0..5), Ddagger(mu=0..5), chi(mu=0..5). Full real-space joint matrices use three blocks of 54 faces in that order.
* At isotropy, a geometric first-vertex gauge makes D translation covariant before primitive unfolding. Each of three primitive folds has two up/down-face channels, each normalized over 27 faces. The 6x6 unfolding is unitary and cross-fold blocks vanish. Physical K=(1/3,2/3) folds to cell Gamma; primitive Gamma is different. No exact M points exist here.
* Magnetic S-dot(Q)=sum_ij exp[-2pi i Q.(r_j-r_i)]<Si.Sj>/27, including three spin components. Its absolute scale is not directly interchangeable with four-/six-spin fields.

Raw per-state 8x8 matrices and full connected covariance are preserved. Isotropic_ground_manifold files separately project the verified K/Kprime pair and retain a C6-character-3 ground-basis alternative. Ensemble restoration and ground-basis choice must not be mistaken for a phase transition. Full connected covariance equals averaged within-image covariance plus covariance of image means; that decomposition is ensemble-dependent. No artificial broadening or physical linewidth is used.

## Limitations

All requested grid points are symmetry screened. Coarse spin-sector changes are bracketed by the sampled couplings, not an exhaustive crossing search. Exclusion bounds are conditional on numerically verified lowest roots, not formal interval arithmetic. The reused isotropic checkpoint has measured H residual 8.30e-9 and S-squared residual 1.72e-8; new ground states are tighter. Local identities, RDM positivity and covariance checks are much tighter than 1e-8. Stronger dynamical resonance/phase claims require further evidence and size scaling.
