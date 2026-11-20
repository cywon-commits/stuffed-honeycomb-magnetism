"""PD wedge: for J_AB=1, J_BC=j+d/2, J_CA=j-d/2 find the smallest d where ground S jumps (18 spins)."""
import numpy as np, csv
from ed3 import sectors, solve
from ed_sz import geometry
shapes = geometry.catalog()
rows = []
for name in ['K18', 'C18']:
    secs = sectors(shapes[name])
    def S_of(j, d):
        e0 = {m: solve(s, (1.0, j + d / 2, j - d / 2))[0][0] for m, s in secs.items()}
        mn = min(e0.values()); return max(m for m in e0 if e0[m] < mn + 1e-8)
    for j in [0.08, 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25, 0.3]:
        S0 = S_of(j, 0.0)
        lo, hi = 0.0, min(2 * j, 0.4) * 0.999
        if S_of(j, hi) == S0: rows.append((name, j, S0, None, None)); print(name, j, S0, 'no jump'); continue
        for _ in range(14):
            m = (lo + hi) / 2
            if S_of(j, m) == S0: lo = m
            else: hi = m
        Sj = S_of(j, hi)
        rows.append((name, j, S0, (lo + hi) / 2, Sj))
        print(f"{name} j={j:.3f} S0={S0} d_c={(lo+hi)/2:.5f} -> S={Sj}  classical wedge j^2={j*j:.4f} (exact {j*j - 0*j:.4f})", flush=True)
with open('../../data/field_dice/pd_wedge_18.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['shape', 'j', 'S_mirror', 'delta_c', 'S_after']); w.writerows(rows)
