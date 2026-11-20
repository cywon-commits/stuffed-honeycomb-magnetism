"""Bulk linear-spin-wave zero-point gains g = E_LSWT/E_cl - 1 (S = 1/2) for the 2D pieces
(honeycomb, square, Lieb, dice) and the triangular 120-degree state, for the resonance-efficiency
formula f_P = (1 + g_tri) / [2 (1 + g_P)]. Collinear pieces: generic Bogoliubov on a unit cell."""
import numpy as np, sys, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).parent))
S = 0.5


def collinear(sig, bonds, n=96):
    """sig: +-1 per basis site; bonds: (i, j, (R1, R2)) antiparallel, J = 1. Returns (E_cl, E_LSWT) per site."""
    nb = len(sig); z = np.zeros(nb)
    for i, j, _ in bonds:
        z[i] += 1; z[j] += 1
    Ecl = -S * S * len(bonds) / nb
    ETA = np.diag([1] * nb + [-1] * nb); corr = 0.0
    g = (np.arange(n) + 0.5) / n
    for k1 in g:
        for k2 in g:
            A = np.diag(S * z).astype(complex); B = np.zeros((nb, nb), complex)
            for i, j, R in bonds:
                ph = np.exp(2j * np.pi * (k1 * R[0] + k2 * R[1]))
                B[i, j] += S * ph; B[j, i] += S * ph
            Am = A.copy(); M = np.block([[A, B], [B.conj().T, Am.conj()]])
            w, v = np.linalg.eig(ETA @ M)
            nr = np.real(np.einsum('ij,ii,ij->j', v.conj(), ETA, v))
            om = np.real(w[nr > 1e-12])
            if len(om) != nb:
                om = np.sort(np.abs(np.real(w)))[nb:]
            corr += 0.5 * (om.sum() - np.real(np.trace(A)))
    return Ecl, Ecl + corr / n**2 / nb


PIECES = {
    # honeycomb: A(0) B(1); B at A + (1/3,1/3)-type positions: bonds to cells 0, -a1, -a2
    "honeycomb": ([1, -1], [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1))]),
    "square": ([1, -1], [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1)), (0, 1, (-1, -1))]),
    # Lieb: corner 0, edge 1 (x-bond), edge 2 (y-bond)
    "Lieb": ([1, -1, -1], [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 2, (0, 0)), (0, 2, (0, -1))]),
    # dice: hub 0 (6 bonds), rims 1, 2 (3 bonds each), triangular Bravais lattice
    "dice": ([1, -1, -1], [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1)),
                           (0, 2, (0, 0)), (0, 2, (1, 0)), (0, 2, (0, 1))]),
}

if __name__ == "__main__":
    from plateau_1S import zero_point
    t = [np.array([np.sin(a), 0, np.cos(a)]) for a in (0, 2 * np.pi / 3, 4 * np.pi / 3)]
    dE, _ = zero_point(t, (1.0, 1.0, 1.0), 0.0, S, n=96)
    Ecl_t = -1.5 * S * S; g_tri = (Ecl_t + dE) / Ecl_t - 1
    print(f"triangle: E_cl={Ecl_t:.4f} E_LSWT={Ecl_t+dE:.4f} g={g_tri:.4f}")
    for name, (sig, bonds) in PIECES.items():
        Ec, El = collinear(sig, bonds)
        g = El / Ec - 1
        print(f"{name:9s}: E_cl/site={Ec:.4f} E_LSWT/site={El:.4f} g={g:.4f} f_P(LSWT)={(1+g_tri)/(2*(1+g)):.4f}")
