"""Finite-size theory of the honeycomb-mirror plateau: rim Neel rotor with an order-by-disorder easy axis.
H = L^2/(2I) - D n_z^2,  I = chi N (chi = 1/9, classical, per site),  D = N K (K = LSWT selection energy per site).
Sector L_z = m. Plateau width of the cluster W_N = 2 (E_1 - E_0) (h_eff symmetric about the plateau centre).
Limits: D -> 0: W_N = 2/I (tower); D -> inf: W_N -> 2 sqrt(2D/I) = 2 sqrt(2K/chi) (bulk spin flop).
W_N * I and q = <n_z^2>_{m=0} depend only on x = D I = chi K N^2.
"""
import numpy as np, csv

def levels(x, m, lmax=80):
    """Lowest eigenvalue/eigenvector of L^2/2 - x n_z^2 (units I = 1, D = x) in sector L_z = m."""
    ls = np.arange(abs(m), lmax + 1)
    a = lambda l: np.sqrt((l * l - m * m) / ((2 * l + 1) * (2 * l - 1.0)))
    C = np.zeros((len(ls), len(ls)))
    for i in range(1, len(ls)):
        C[i, i - 1] = C[i - 1, i] = a(ls[i])
    H = np.diag(ls * (ls + 1) / 2.0) - x * C @ C
    w, v = np.linalg.eigh(H)
    q = v[:, 0] @ (C @ C) @ v[:, 0]
    return w[0], q

def G(x):
    """W_N * I and q(m=0)."""
    e0, q = levels(x, 0); e1, _ = levels(x, 1)
    return 2 * (e1 - e0), q

if __name__ == '__main__':
    for x in (0, 0.1, 0.3, 1, 3, 10, 30, 100):
        g, q = G(x); print(f"x={x:6.1f} W*I={g:.4f} (limits 2, {2*np.sqrt(2*x):.3f}) q={q:.4f}")
