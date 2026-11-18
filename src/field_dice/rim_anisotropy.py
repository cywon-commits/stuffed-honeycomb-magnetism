"""ED check of the order-by-disorder mechanism: rim (A-B) bond correlations, longitudinal (zz)
vs transverse (xx+yy), in the lowest state of S^z = N/6 - 1, N/6, N/6 + 1 on the honeycomb mirror."""
import numpy as np, csv
from ed3 import sectors, solve
from ed_sz import geometry
rows = []
for name in ['K18', 'C18']:
    g = geometry.catalog()[name]; secs = sectors(g)
    ab = [(a, b) for a, b, t in g['bonds'] if t == 0]
    for r in (0.4, 0.6, 0.8, 1.0):
        J = (1.0, r, r)
        for m in (2, 3, 4):
            s = secs[m]; _, v = solve(s, J); v = v[:, 0]
            p = np.abs(v) ** 2
            zz = np.mean([p @ np.where(((s.s >> a) & 1) == ((s.s >> b) & 1), 0.25, -0.25) for a, b in ab])
            tot = v @ (s.H['AB'] @ v) / len(ab)
            xy = tot - zz
            rows.append((name, r, m, zz, xy))
            print(f"{name} r={r:.1f} Sz={m}: zz={zz:+.4f} xx+yy={xy:+.4f} zz/(xx+yy)={zz/xy:5.2f}", flush=True)
with open('../../data/field_dice/rim_anisotropy_18.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['shape', 'r', 'Sz', 'zz_AB', 'xxyy_AB']); w.writerows(rows)
