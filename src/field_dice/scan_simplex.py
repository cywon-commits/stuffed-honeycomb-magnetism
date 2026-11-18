"""18-spin ED over the S3 fundamental domain (eps, phi): ground spin, plateau edges,
h_sat, frustration index f = E0 / sum_k J_k E_hex, honeycomb bond satisfactions s_k."""
import csv, time, numpy as np, scipy.sparse.linalg as sla
from ed_sz import geometry, Sector, states_with
from simplex import couplings, eps_max, h_sat_formula, f_classical
OUT = '../../data/field_dice/'
from ed3 import sectors, solve
pts = []
for phd in range(0, 61, 10):
    ph = np.radians(phd); em = eps_max(ph)
    for e in list(np.arange(0, em - 1e-9, 0.1)) + [0.97 * em]:
        pts.append((round(e, 4), phd))
rows = []
shapes = geometry.catalog()
for name in ['K18', 'C18']:
    secs = sectors(shapes[name]); t0 = time.time()
    ehex = min(solve(s, (1, 0, 0))[0][0] for s in secs.values())
    for e, phd in pts:
        j = couplings(e, np.radians(phd)); E = {}; V = {}
        for m, s in secs.items():
            w, v = solve(s, j, k=3 if m == 3 else 1); E[m] = w; V[m] = v[:, 0]
        e0 = {m: E[m][0] for m in E}; emin = min(e0.values())
        S = max(m for m in e0 if e0[m] < emin + 1e-8)
        mg = min(e0, key=e0.get); v = V[mg]
        sk = [float(v @ (secs[mg].H[b] @ v)) / ehex for b in ['AB', 'CB', 'CA']]
        lo = max((e0[3] - e0[m]) / (3 - m) for m in range(3)); hi = min((e0[m] - e0[3]) / (m - 3) for m in range(4, 10))
        rows.append(dict(shape=name, eps=e, phi_deg=phd, J_AB=j[0], J_BC=j[1], J_CA=j[2], S=S,
            f=emin / (sum(j) * ehex), f_cl=f_classical(*j), s_AB=sk[0], s_BC=sk[1], s_CA=sk[2],
            h_lo=max(lo, 0), h_hi=hi, width=hi - max(lo, 0), h_sat=e0[9] - e0[8], h_sat_formula=h_sat_formula(*j),
            plat_gap1=E[3][1] - E[3][0], plat_gap2=E[3][2] - E[3][0]))
    print(name, len(pts), 'points', round(time.time() - t0, 1), 's', flush=True)
with open(OUT + 'simplex_scan_18.csv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
