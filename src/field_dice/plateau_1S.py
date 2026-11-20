"""1/S width of the 1/3 plateau for the three-coupling model (J_AB, J_BC, J_CA).

Method (canonical / Legendre form of Chubukov-Golosov):
  e(m) = e_cl(m) + dE(m) + O(S^0),  h(m) = de/dm,
where dE(m) is the LSWT zero-point energy of the classical canted state with magnetization m.
If the classical window is zero (honeycomb mirror phi = 0, and the triangular point), e_cl'(m)
is continuous at m = S/3, and the plateau edges at order S^0 are
  h_c1 = h_cl + lim_{m -> S/3^-} dE'(m),   h_c2 = h_cl + lim_{m -> S/3^+} dE'(m).
Energies and m are per site. Spin S is a parameter; results are reported as the S^0 coefficients.
"""
import sys
import numpy as np
from scipy.optimize import minimize
from simplex import coupling_matrix

Z = np.array([0.0, 0.0, 1.0])


def vecs(v):
    return [np.array([np.sin(v[2*i]) * np.cos(v[2*i+1]), np.sin(v[2*i]) * np.sin(v[2*i+1]), np.cos(v[2*i])]) for i in range(3)]


def e_cell(v, J, h, S):
    a, b, c = vecs(v)
    jab, jbc, jca = J
    return 3 * S * S * (jab * a @ b + jbc * b @ c + jca * c @ a) - h * S * (a[2] + b[2] + c[2])


def classical(J, h, S, seeds=40, rng=np.random.default_rng(1)):
    best = None
    for _ in range(seeds):
        x0 = rng.uniform(0, np.pi, 6); x0[1::2] *= 2
        r = minimize(e_cell, x0, args=(J, h, S), method='BFGS', options={'gtol': 1e-12})
        if best is None or r.fun < best.fun - 1e-12:
            best = r
    return vecs(best.x), best.fun


def frame(e3):
    t = np.array([1.0, 0, 0]) if abs(e3[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = t - (t @ e3) * e3; e1 /= np.linalg.norm(e1)
    return e1, np.cross(e3, e1)


def blocks(spins, J, h, S, k):
    Jk = coupling_matrix(k, *J); J0 = coupling_matrix((0, 0), *J).real
    u = []
    for e3 in spins:
        e1, e2 = frame(e3); u.append(e1 + 1j * e2)
    U = np.array(u); E3 = np.array(spins)
    uu_bar = U @ U.conj().T          # u_s . conj(u_t)
    uu = U @ U.T                     # u_s . u_t
    loc = h * E3[:, 2] - S * (J0 * (E3 @ E3.T)).sum(1)
    A = 0.5 * S * Jk * uu_bar + np.diag(loc)
    B = 0.5 * S * Jk * uu
    return A, B


ETA = np.diag([1, 1, 1, -1, -1, -1])


def zero_point(spins, J, h, S, n=48):
    g = (np.arange(n) + 0.5) / n; corr = 0.0; wmin = np.inf
    for a in g:
        for b in g:
            A, B = blocks(spins, J, h, S, (a, b)); Am, _ = blocks(spins, J, h, S, (-a, -b))
            M = np.block([[A, B], [B.conj().T, Am.conj()]])
            w, v = np.linalg.eig(ETA @ M)
            nr = np.real(np.einsum('ij,ii,ij->j', v.conj(), ETA, v))
            om = np.real(w[nr > 1e-9])
            if len(om) != 3:            # instability or degeneracy: fall back to sorting
                om = np.sort(np.real(w))[3:]
            wmin = min(wmin, om.min())
            corr += 0.5 * (om.sum() - np.real(np.trace(A)))
    return corr / n**2 / 3, wmin      # per site


def branch(J, S, hs, n=48):
    out = []
    for h in hs:
        sp, ecl = classical(J, h, S)
        m = S * sum(s[2] for s in sp) / 3
        dE, wmin = zero_point(sp, J, h, S, n)
        out.append((h, m, ecl / 3 + h * m, dE, wmin))
    return np.array(out)


def edges(J, S=0.5, hc=None, d=(0.02, 0.04, 0.06, 0.08), n=48):
    """Return (h_c1, h_c2) at order S^0 by linear extrapolation of dE'(m) to m = S/3 from each side."""
    res = {}
    for side, sgn in (('lo', -1), ('hi', +1)):
        hs = hc + sgn * S * np.array(d)
        br = branch(J, S, hs, n)
        # slopes dE/dm between consecutive points, attributed to midpoints, extrapolated to m = S/3
        m, dE = br[:, 1], br[:, 3]
        mid = (m[1:] + m[:-1]) / 2; sl = np.diff(dE) / np.diff(m)
        p = np.polyfit(mid - S / 3, sl, 1)
        res[side] = (np.polyval(p, 0.0), br)
    return hc + res['lo'][0], hc + res['hi'][0], res


if __name__ == '__main__':
    S = float(sys.argv[1]) if len(sys.argv) > 1 else 0.5
    for r in [float(x) for x in sys.argv[2:]] or [1.0]:
        J = (1.0, r, r)
        h1, h2, res = edges(J, S, hc=3 * S * r)
        print(f"r={r:.3f} S={S}: h_c1={h1:.4f} h_c2={h2:.4f} width={h2-h1:.4f} "
              f"(classical point {3*S*r:.4f}); min omega lo/hi {res['lo'][1][:,4].min():.3g} {res['hi'][1][:,4].min():.3g}")
