"""Independent N=27 reference energies (no symmetry reduction) for high S^z sectors,
used as unit tests for the three-coupling translation-reduced 27-spin solver."""
import csv, time, numpy as np
from ed_sz import geometry
from ed3 import sectors, solve
from simplex import couplings
g = geometry.catalog()['N27']; N = g['N']
pts = {'T1_mirror_r0.6': (1.0, 0.6, 0.6), 'T2_eps0.2_phi30': tuple(np.round(couplings(0.2, np.radians(30)), 12)),
       'T3_dice_JCA0': (1.0, 1.0, 0.0), 'T4_wedge_j0.2_d0.08': (1.0, 0.24, 0.16), 'T5_triangle': (1.0, 1.0, 1.0)}
secs = {}
from ed_sz import Sector, states_with
split = {'AB': [(a, b) for a, b, t in g['bonds'] if t == 0],
         'CA': [(a, b) for a, b, t in g['bonds'] if t == 1 and b % 3 == 0],
         'CB': [(a, b) for a, b, t in g['bonds'] if t == 1 and b % 3 == 1]}
rows = []
for nup in range(20, 28):          # S^z = 13/2 ... 27/2
    t = time.time(); s = Sector.__new__(Sector); s.N = N; s.nup = nup; s.s = states_with(N, nup)
    s.H = {k: s._build(v) for k, v in split.items()}
    for name, J in pts.items():
        e = solve(s, J)[0][0]
        rows.append(dict(point=name, J_AB=J[0], J_BC=J[1], J_CA=J[2], Sz=f'{nup - 13.5:.1f}', n_up=nup, E_min=f'{e:.12f}'))
    print(nup, len(s.s), round(time.time() - t, 1), 's', flush=True)
with open('../../data/field_dice/reference_N27_highSz.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
