import csv, numpy as np
from lswt3 import zero_point
Ehc, Ehq = zero_point((1, -1, 1), (1.0, 0.0, 0.0))
print(f"honeycomb per cell (J_AB=1): classical {Ehc:.4f}, LSWT {Ehq:.4f}  (per site of honeycomb: {Ehq/2:.4f}; known LSWT -0.5258)")
R = list(csv.DictReader(open('../../data/field_dice/simplex_scan_18.csv')))
print("phi eps  S | f_ED   f_cl   f_LSWT | LSWT/cl  ED/cl")
out = []
for x in R:
    if x['shape'] != 'K18' or x['S'] != '3': continue
    J = (float(x['J_AB']), float(x['J_BC']), float(x['J_CA']))
    k = int(np.argmin(J)); sig = [(1, 1, -1), (-1, 1, 1), (1, -1, 1)][k]
    Ec, Eq = zero_point(sig, J, n=24)
    fl = Eq / (sum(J) * Ehq); fc = Ec / (sum(J) * Ehc)
    out.append((int(x['phi_deg']), float(x['eps']), float(x['f']), fc, fl))
    print(f"{x['phi_deg']:>3} {float(x['eps']):.3f} | {float(x['f']):.3f} {fc:.3f} {fl:.3f} | {fl/fc:.3f} {float(x['f'])/fc:.3f}")
