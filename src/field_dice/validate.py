import csv, numpy as np, time
from ed_sz import geometry, Sector
shapes = geometry.catalog()
ref = list(csv.DictReader(open('../../data/workspaces/finite_size_18_27_36pilot/gs_18_scan.csv')))
for name in ['K18','C18']:
    g = shapes[name]; t=time.time()
    secs = {n: Sector(g, n) for n in range(9, 19)}
    print(name, 'sectors built', round(time.time()-t,1),'s')
    for r in [0.2, 0.4, 0.5, 0.6, 0.7, 1.0]:
        E = {n-9: secs[n].lowest(1.0, r)[0][0] for n in secs}
        e0 = min(E.values()); S = max(m for m in E if E[m] < e0 + 1e-8)
        row = [x for x in ref if x['shape']==name and abs(float(x['r'])-r)<1e-9][0]
        print(f"  r={r}: E0={e0:.12f} ref={float(row['energy']):.12f} diff={e0-float(row['energy']):.1e} S={S} refS={row['S']}")
