"""Linear spin-wave (LSWT) and one-magnon results for the stuffed-honeycomb
Heisenberg antiferromagnet in a field.

H = J sum_rim S_i.S_j + J' sum_{rim,hub} S_i.T_R - h sum S^z,  J, J' >= 0.

Sublattices: A=0, B=1 (honeycomb rim), H=2 (hub).  Bond rules follow
src/workspaces/finite_size_18_27_36pilot/geometry.py:
  rim:  A_R - B_R, A_R - B_(R-a1), A_R - B_(R-a2)
  hub:  H_R - {A_R, B_R, A_(R+a1), B_(R+a1-a2), A_(R+a1-a2), B_(R-a2)}
Momentum k is given in cell reciprocal coordinates (k1, k2), phase exp(2 pi i k.dR).
"""
import numpy as np

RIM = [(0, 1, (0, 0)), (0, 1, (-1, 0)), (0, 1, (0, -1))]          # (i, j, R_j - R_i)
HUB = [(2, 0, (0, 0)), (2, 1, (0, 0)), (2, 0, (1, 0)),
       (2, 1, (1, -1)), (2, 0, (1, -1)), (2, 1, (0, -1))]


def coupling(k, J, Jp):
    """3x3 Hermitian Fourier exchange matrix J_ab(k) (each bond counted once per pair)."""
    M = np.zeros((3, 3), complex)
    for (i, j, d), c in [(b, J) for b in RIM] + [(b, Jp) for b in HUB]:
        ph = np.exp(2j * np.pi * (k[0] * d[0] + k[1] * d[1]))
        M[i, j] += c * ph
        M[j, i] += c * np.conj(ph)
    return M


def site_sum(J, Jp):
    """sum_j J_ij for each sublattice (A, B, H)."""
    return np.array([3 * J + 3 * Jp, 3 * J + 3 * Jp, 6 * Jp])


def kgrid(n=60):
    g = np.arange(n) / n
    pts = [(a, b) for a in g for b in g]
    pts += [(0, 0), (1 / 3, 2 / 3), (2 / 3, 1 / 3), (0.5, 0), (0, 0.5), (0.5, 0.5)]
    return pts


def h_sat(J, Jp, S=0.5, n=60):
    """Exact saturation field: one-magnon instability of the fully polarized state.
    omega(k) = h + S * eig[J(k) - diag(sum_j J_ij)];  h_sat = -S * min_k eig."""
    D = np.diag(site_sum(J, Jp))
    best = (np.inf, None)
    for k in kgrid(n):
        w = np.linalg.eigvalsh(coupling(k, J, Jp) - D).min()
        if w < best[0]:
            best = (w, k)
    return -S * best[0], best[1]


def collinear_stable(sig, h, J, Jp, S=0.5, n=24, tol=1e-9):
    """Classical (LSWT-level) local stability of a collinear state with sublattice
    signs sig = (sA, sB, sH) along the field.  The quadratic boson Hamiltonian
    H2 = sum_k [a^+ A a + (1/2)(a B a + h.c.)] must be positive semidefinite.
    A_ii = h s_i - S sum_j J_ij s_i s_j ;  parallel pair -> hopping S J_ij (A),
    antiparallel pair -> pairing S J_ij (B)."""
    sig = np.array(sig, float)
    for k in kgrid(n):
        Jk = coupling(k, J, Jp)
        Jm = coupling([-k[0], -k[1]], J, Jp)
        par = np.outer(sig, sig) > 0
        A = np.where(par, S * Jk, 0)
        B = np.where(~par, S * Jk, 0)
        loc = np.array([h * sig[i] - S * sum(
            (coupling((0, 0), J, Jp)[i, j]) * sig[i] * sig[j] for j in range(3)) for i in range(3)])
        A = A + np.diag(loc)
        Am = np.where(par, S * Jm, 0) + np.diag(loc)
        M = np.block([[A, B], [B.conj().T, Am.conj()]])
        if np.linalg.eigvalsh((M + M.conj().T) / 2).min() < -tol:
            return False
    return True


def window(sig, J, Jp, S=0.5, hmax=6.0, dh=0.0025):
    """Field interval in which the collinear state sig is locally stable (LSWT)."""
    hs = np.arange(0, hmax, dh)
    ok = [h for h in hs if collinear_stable(sig, h, J, Jp, S)]
    return (min(ok), max(ok)) if ok else None


if __name__ == "__main__":
    print("r, h_sat, k_min, h_sat formula S*max(6+3r, 9r)")
    for r in [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0]:
        hs, k = h_sat(1.0, r)
        print(f"{r:5.2f} {hs:8.4f} {np.round(k, 3)} {0.5 * max(6 + 3 * r, 9 * r):8.4f}")
