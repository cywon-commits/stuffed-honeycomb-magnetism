"""Bulk partial-disorder wedge from the Seifert-Vojta effective hub model (PRB 99, 155156, Table I).
Physical per-pair couplings = 2 x table totals (verified: hub_eff.py reproduces 2 j^{xy,0} with the
exact sum rule sum_R K^xx = 9 chi_perp). Criterion: hub saturation field of the XXZ model in a
longitudinal field equals the static field 3 m_N (J_BC - J_CA) from the rim Neel order."""
import numpy as np
S = 0.5
jxy0 = [-0.153987, 0.027187, -0.016334, 0.001881, -0.002153, 0.000684, 0.000021, -0.000385, 0.000150]
Sjxy1 = [0.0745365, -0.0131602, 0.0079060, -0.0009108, 0.0010418, -0.0003309, -0.0000099, 0.0001868, -0.0000716]
Sjz = [-0.0114266, 0.0235538, -0.0005096, 0.0022996, 0.0007775, 0.0008679, 0.0006570, 0.0004147, 0.0003788]
jxy = [2 * (a + b / S) for a, b in zip(jxy0, Sjxy1)]          # per unordered pair, S = 1/2
jz = [2 * b / S for b in Sjz]
# triangular-lattice shells (d^2 in units of the hub spacing) and their neighbour vectors
shells = [1, 3, 4, 7, 9, 12, 13, 16, 19]
vecs = {d2: [(u, v) for u in range(-6, 7) for v in range(-6, 7) if u * u + v * v + u * v == d2] for d2 in shells}
def J_q(c, q):
    return sum(cn * np.cos(2 * np.pi * (q[0] * u + q[1] * v)) for cn, d2 in zip(c, shells) for (u, v) in vecs[d2])
print("neighbour counts:", [len(vecs[d]) for d in shells])
g = np.linspace(0, 1, 121, endpoint=False)
Jxy_min = min(J_q(jxy, (a, b)) for a in g for b in g)
Jz0 = J_q(jz, (0, 0))
hsat = S * (Jz0 - Jxy_min)
print(f"sum_j Jz = {Jz0:.4f}, min_q Jxy(q) = {Jxy_min:.4f}, hub h_sat = {hsat:.4f} J'^2/J")
for mN, lab in [(0.2677, 'QMC'), (0.2418, 'LSWT')]:
    c = hsat / (3 * mN)
    print(f"m_N = {mN} ({lab}): Delta J_c = {c:.3f} J'^2/J ->", ", ".join(f"j={j}: {c*j*j:.4f}" for j in [0.15, 0.2, 0.3]))
# classical-order (S^0 only) transverse coupling for comparison
print("S^0 transverse NN per pair:", 2 * jxy0[0], "| S=1/2 total:", round(jxy[0], 6), "| cancellation", round(1 - jxy[0] / (2 * jxy0[0]), 3))
