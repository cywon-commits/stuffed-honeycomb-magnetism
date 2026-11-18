"""LSWT for collinear states of the three-coupling model (J_AB, J_BC, J_CA), with field h.
Sublattice order (A, B, C); couplings via simplex.coupling_matrix (C = hub index 2)."""
import numpy as np
from simplex import coupling_matrix

S = 0.5


def blocks(sig, h, J, k):
    """A(k), B(k) of H2 = sum_k [a^+ A a + 1/2 (a B a + h.c.)] for collinear signs sig."""
    sig = np.array(sig, float)
    Jk = coupling_matrix(k, *J); J0 = coupling_matrix((0, 0), *J).real
    par = np.outer(sig, sig) > 0
    loc = h * sig - S * (J0 * np.outer(sig, sig)).sum(1)
    A = np.where(par, S * Jk, 0) + np.diag(loc)
    B = np.where(~par, S * Jk, 0)
    return A, B


def stable(sig, h, J, n=18, tol=1e-9):
    g = np.arange(n) / n
    for a in g:
        for b in g:
            A, B = blocks(sig, h, J, (a, b))
            Am, _ = blocks(sig, h, J, (-a, -b))
            M = np.block([[A, B], [B.conj().T, Am.conj()]])
            if np.linalg.eigvalsh((M + M.conj().T) / 2).min() < -tol:
                return False
    return True


def window(sig, J, hmax=12, n=18):
    hs = np.arange(0, hmax, 0.05); st = [stable(sig, h, J, n) for h in hs]
    if not any(st):
        return None
    def edge(lo, hi):
        for _ in range(28):
            m = (lo + hi) / 2
            if stable(sig, m, J, n): lo = m
            else: hi = m
        return (lo + hi) / 2
    i = st.index(True); j = len(st) - 1 - st[::-1].index(True)
    a = 0.0 if i == 0 else edge(hs[i], hs[i - 1])
    return a, edge(hs[j], hs[j + 1])


def zero_point(sig, J, h=0.0, n=36):
    """LSWT ground energy per cell: E_cl + (1/2) sum_k [sum_n omega_n(k) - tr A(k)] / Nk.
    Positive-norm Bogoliubov eigenvalues of sigma3 M are the magnon energies."""
    sig = np.array(sig, float)
    J0 = coupling_matrix((0, 0), *J).real
    Ecl = 0.5 * S * S * (J0 * np.outer(sig, sig)).sum() - h * S * sig.sum()
    g = (np.arange(n) + 0.5) / n; corr = 0.0
    eta = np.diag([1, 1, 1, -1, -1, -1])
    for a in g:
        for b in g:
            A, B = blocks(sig, h, J, (a, b)); Am, _ = blocks(sig, h, J, (-a, -b))
            M = np.block([[A, B], [B.conj().T, Am.conj()]])
            w, v = np.linalg.eig(eta @ M)
            norms = np.real(np.einsum('ij,ii,ij->j', v.conj(), eta, v))
            om = np.real(w[norms > 1e-9])
            corr += 0.5 * (om.sum() - np.real(np.trace(A)))
    return Ecl, Ecl + corr / n ** 2
