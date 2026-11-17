# Static correlation results (Stages 1–2)

Completed coupling ratios: [0.2, 0.5, 1.0]. All results use the 27-site torus, Sz=+1/2 and J=1. Primary scalar matrices use an explicitly stated SU(2) spin-multiplet and exact point-group average; selected-state raw/connected matrices are also provided. No dynamical spectra were calculated.

| Jprime/J | E/J | cell momentum representative | S | rim m_N squared | hub S(Gamma) | 120-degree intensity | energy residual |
|---:|---:|---|---:|---:|---:|---:|---:|
| 0.2 | -10.185949128477 | (1/3,2/3) | 0.5 | 0.24532220 | 0.09517246 | 2.24468141 | 3.28e-14 |
| 0.5 | -11.083587020439 | (0,0) | 2.5 | 0.22257082 | 2.53725956 | 3.25999050 | 3.16e-14 |
| 1 | -15.125973385092 | (0,0) | 0.5 | 0.16383234 | 1.97524731 | 2.94898207 | 8.30e-09 |


From isotropy to Jprime/J=0.2, the rim Neel squared amplitude changes from 0.16383234 to 0.24532220; hub S(Gamma) changes from 1.97524731 to 0.09517246. At Jprime/J=0.5 the corresponding values are 0.22257082 and 2.53725956. These are equal-time finite-size correlations, with all definitions and state averages held fixed.

- Jprime/J=0.2: hub peak S=0.99531166 at cell momentum numerator(s) [[0, 1], [0, 2], [1, 0], [1, 1], [2, 0], [2, 2]] over 3; rim NN=-0.37033550, hub-rim bond=-0.01730469, all-distinct hub-pair mean=-0.08185344; uniform rim/hub Gamma cross=-0.01208820.
- Jprime/J=0.5: hub peak S=2.53725956 at cell momentum numerator(s) [[0, 0]] over 3; rim NN=-0.32302122, hub-rim bond=-0.08748200, all-distinct hub-pair mean=0.22340745; uniform rim/hub Gamma cross=-1.08303935.
- Jprime/J=1: hub peak S=1.97524731 at cell momentum numerator(s) [[0, 0]] over 3; rim NN=-0.18674041, hub-rim bond=-0.18674041, all-distinct hub-pair mean=0.15315591; uniform rim/hub Gamma cross=-1.94746953.


The intermediate point has a genuine S=5/2 ground multiplet, verified both by S-squared and by global-spin raising norms 8, 5, then zero at the highest weight; all raised partners retain the ground energy. Its selected Sz=+1/2 state has oppositely signed rim and hub polarizations. These are partially compensated finite-cluster moments, not proof of a thermodynamic ferrimagnetic phase. Both endpoints have S=1/2.

At isotropy, the uniform hub structure factor is large despite total S=1/2. The rim uniform component compensates the hub moment; hub Gamma weight alone is not a net-ferromagnetism diagnostic. The physical 120-degree probe includes its required 1/3 sublattice normalization. Away from isotropy it is a Fourier diagnostic that can also receive rim-Neel weight, not a unique identifier of 120-degree order.

Inversion enforces zero rim-staggered/hub cross correlation at Gamma in the averaged ensemble. This selection rule does not establish absence of a hub coupling. The nonzero uniform rim/hub cross channel and the full complex matrices at finite momenta provide the relevant complementary information. Correlation eigenspaces, including inversion-resolved Gamma modes and hub weights, are saved; arbitrary vectors inside degenerate eigenspaces are not compared as unique modes.

Ground-sector searches at each unequal coupling cover all nine cell momenta using exact point-group equivalences and every compatible cyclic character. Direct H residuals are in the table; total S is verified independently by applying S-squared. Small-cluster tests independently validate z,+,- seeds, H acting on the seeds, scalar Gram matrices and the reduced spectra. Production data pass positivity, Hermiticity, energy reconstruction, onsite sum rules, fixed-Sz conservation and inversion tests. CPU/RAM, H-product counts and wall times are recorded in sector_convergence.csv and technical_report.md.

These are finite-size equal-time correlations. Three couplings on one torus cannot establish thermodynamic long-range order, a transition, order destruction or damping. No exact M point exists on this torus.

For the next dynamical stage, retain Jprime/J=1 as a validation anchor and prioritize the contrast between 0.5 (spin 5/2, uniform hub peak) and 0.2 (spin 1/2, finite-momentum hub peak). Start at cell Gamma with the full n,m,h channel matrix, testing the conserved inelastic m+h null channel; also target the allowed hub-peak momentum stars listed above. At isotropy resolve physical K/Kprime using basis phases even though both fold to cell Gamma. Separate the entire degenerate zero-energy ground manifold before interpreting inelastic weight. Large static weight alone does not prove a low-frequency peak. Stage 3 remains unstarted.

See SPEC_REVIEW.md for conventions, technical_report.md for validation/resources, static_summary.csv for key numbers and static_j*_C.csv for every 3x3 matrix at all nine cell momenta.
