"""Three-coupling ED helpers: sectors with H_AB, H_CB, H_CA split; lowest levels."""
import numpy as np, scipy.sparse.linalg as sla
from ed_sz import Sector, states_with
def sectors(g):
    N = g['N']; split = {'AB': [(a, b) for a, b, t in g['bonds'] if t == 0],
        'CA': [(a, b) for a, b, t in g['bonds'] if t == 1 and b % 3 == 0],
        'CB': [(a, b) for a, b, t in g['bonds'] if t == 1 and b % 3 == 1]}
    secs = {}
    for n in range(N // 2, N + 1):
        s = Sector.__new__(Sector); s.N = N; s.nup = n; s.s = states_with(N, n)
        s.H = {k: s._build(v) for k, v in split.items()}; secs[n - N // 2] = s
    return secs
def solve(s, j, k=1):
    Hm = j[0] * s.H['AB'] + j[1] * s.H['CB'] + j[2] * s.H['CA']
    if Hm.shape[0] < 400:
        w, v = np.linalg.eigh(Hm.toarray()); return w[:k], v[:, :k]
    w, v = sla.eigsh(Hm, k=k, which='SA', tol=1e-11, ncv=30); o = np.argsort(w); return w[o], v[:, o]
