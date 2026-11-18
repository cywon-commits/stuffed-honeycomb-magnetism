"""S3 fundamental domain of the three-coupling model.

J_k = Jbar [1 + 2 eps cos(phi - 2 pi k/3)],  k = 0 (AB), 1 (BC), 2 (CA), Jbar = 1.
phi = 0: J_AB largest, J_BC = J_CA  (honeycomb-direction mirror, the old r < 1 line)
phi = pi/3: J_AB = J_BC > J_CA      (dice-direction mirror, equivalent to the old r > 1 line; hub = B)
eps_max(phi) = -1 / (2 cos(phi + 2 pi/3))  (J_CA = 0: bipartite Lieb-Mattis edge)
"""
import numpy as np
from scipy.optimize import minimize

RIM = [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1))]           # A-B
HUB = [(2, 0, (0, 0)), (2, 1, (0, 0)), (2, 0, (1, 0)),
       (2, 1, (1, -1)), (2, 0, (1, -1)), (2, 1, (0, -1))]           # C-A, C-B


def couplings(eps, phi):
    return tuple(1 + 2 * eps * np.cos(phi - 2 * np.pi * k / 3) for k in range(3))


def eps_max(phi):
    return -1 / (2 * np.cos(phi + 2 * np.pi / 3))


def to_eps_phi(jab, jbc, jca):
    jb = (jab + jbc + jca) / 3
    d = np.array([jab, jbc, jca]) / jb - 1
    x = (d[0] - 0.5 * d[1] - 0.5 * d[2]) / 3
    y = (np.sqrt(3) / 2 * (d[1] - d[2])) / 3
    return jb, np.hypot(x, y), np.arctan2(y, x)


def coupling_matrix(k, jab, jbc, jca):
    M = np.zeros((3, 3), complex)
    for i, j, d in RIM:
        ph = np.exp(2j * np.pi * (k[0] * d[0] + k[1] * d[1])); M[i, j] += jab * ph; M[j, i] += jab * np.conj(ph)
    for i, j, d in HUB:
        c = jca if j == 0 else jbc
        ph = np.exp(2j * np.pi * (k[0] * d[0] + k[1] * d[1])); M[i, j] += c * ph; M[j, i] += c * np.conj(ph)
    return M


def h_sat_numeric(jab, jbc, jca, S=0.5, n=48):
    D = np.diag(coupling_matrix((0, 0), jab, jbc, jca).real.sum(1))
    g = np.arange(n) / n
    w = min(np.linalg.eigvalsh(coupling_matrix((a, b), jab, jbc, jca) - D).min() for a in g for b in g)
    return -S * w


def h_sat_formula(jab, jbc, jca, S=0.5):
    """3 S lambda_max(weighted Laplacian of the A-B-C triangle) = 9 S Jbar (1 + eps)."""
    s = jab + jbc + jca
    return 3 * S * (s + np.sqrt(jab**2 + jbc**2 + jca**2 - jab * jbc - jbc * jca - jca * jab))


def classical_cell_energy(v, jab, jbc, jca, h, S=0.5):
    a, b, c = [np.array([np.sin(v[2*i]) * np.cos(v[2*i+1]), np.sin(v[2*i]) * np.sin(v[2*i+1]), np.cos(v[2*i])]) for i in range(3)]
    return 3 * S * S * (jab * a @ b + jbc * b @ c + jca * c @ a) - h * S * (a[2] + b[2] + c[2])


def classical_state(jab, jbc, jca, h, S=0.5, starts=24, seed=0):
    rng = np.random.default_rng(seed); best = None
    for _ in range(starts):
        r = minimize(classical_cell_energy, rng.uniform(0, np.pi, 6) * [1, 2, 1, 2, 1, 2],
                     args=(jab, jbc, jca, h, S), method='BFGS', options={'gtol': 1e-11})
        if best is None or r.fun < best.fun - 1e-12:
            best = r
    v = best.x
    mz = [np.cos(v[0]), np.cos(v[2]), np.cos(v[4])]
    if h == 0:   # moment magnitude of the net cell spin (direction free)
        vec = sum(np.array([np.sin(v[2*i]) * np.cos(v[2*i+1]), np.sin(v[2*i]) * np.sin(v[2*i+1]), np.cos(v[2*i])]) for i in range(3))
        return best.fun, np.linalg.norm(vec) / 3, mz
    return best.fun, np.mean(mz), mz


def f_classical(jab, jbc, jca):
    """E_cl / sum_k J_k E_hex,cl ; canted formula if feasible else collinear Lieb-Mattis."""
    js = sorted([jab, jbc, jca]); j3, j2, j1 = js[0], js[1], js[2]
    if j3 > 0 and j3 >= j1 * j2 / (j1 + j2):
        e = -0.5 * (jab * jbc / jca + jbc * jca / jab + jca * jab / jbc)
    else:
        e = j3 - j1 - j2           # weakest bond parallel, two strongest antiparallel
    return e / -(jab + jbc + jca)


if __name__ == "__main__":
    print("h_sat: numeric over k vs Laplacian formula vs 9S(1+eps)")
    worst = 0
    for phi in np.linspace(0, np.pi / 3, 7):
        for eps in np.linspace(0, eps_max(phi), 6):
            j = couplings(eps, phi); hn = h_sat_numeric(*j); hf = h_sat_formula(*j)
            worst = max(worst, abs(hn - hf), abs(hf - 4.5 * (1 + eps)))
    print("  max deviation over 42 points:", worst)
    print("classical f and M at h = 0 along phi (eps = 0.3 and eps_max):")
    for phi in np.linspace(0, np.pi / 3, 4):
        for eps in [0.3, eps_max(phi) * 0.999]:
            j = couplings(eps, phi); E, M, _ = classical_state(*j, 0)
            print(f"  phi={np.degrees(phi):4.0f} eps={eps:.3f} J={np.round(j,3)} f_cl={f_classical(*j):.4f} "
                  f"E_num/E_formula={E / (-(sum(j)) * f_classical(*j) * 0.25):.6f} M/Ms={M:.4f}")
