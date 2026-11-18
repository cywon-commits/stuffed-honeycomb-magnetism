"""Exact diagonalization of the stuffed-honeycomb Heisenberg model in every S^z sector.

Geometry (sites, bonds, coordinates) is imported from the verified finite-size
workspace (geometry.py).  H = J H_rim + J' H_hub is built once per sector as two
sparse matrices, so a coupling scan costs one Lanczos per point and sector.
Site order: 3*cell + s, s = 0 (A), 1 (B), 2 (H).
"""
import sys, pathlib, itertools
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy.special import comb

ROOT = pathlib.Path(__file__).resolve().parents[1] / "workspaces" / "finite_size_18_27_36pilot"
sys.path.insert(0, str(ROOT))
import geometry  # noqa: E402


def states_with(N, nup):
    out = np.fromiter((sum(1 << i for i in c) for c in itertools.combinations(range(N), nup)),
                      dtype=np.int64, count=int(comb(N, nup, exact=True)))
    out.sort()
    return out


class Sector:
    def __init__(self, g, nup):
        self.N = g["N"]; self.nup = nup
        self.s = states_with(self.N, nup)
        self.H = {}
        for t in (0, 1):
            pairs = [(a, b) for a, b, tt in g["bonds"] if tt == t]
            self.H[t] = self._build(pairs)

    def _build(self, pairs):
        s = self.s; n = len(s); diag = np.zeros(n)
        rows, cols, vals = [], [], []
        for a, b in pairs:
            ba = (s >> a) & 1; bb = (s >> b) & 1
            same = ba == bb
            diag += np.where(same, 0.25, -0.25)
            idx = np.nonzero(~same)[0]
            t = s[idx] ^ ((1 << a) | (1 << b))
            j = np.searchsorted(s, t)
            rows.append(idx); cols.append(j); vals.append(np.full(len(idx), 0.5))
        rows.append(np.arange(n)); cols.append(np.arange(n)); vals.append(diag)
        return sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n))

    def ham(self, J, Jp):
        return J * self.H[0] + Jp * self.H[1]

    def lowest(self, J, Jp, k=1):
        Hm = self.ham(J, Jp); n = Hm.shape[0]
        if n <= 400:
            w, v = np.linalg.eigh(Hm.toarray()); return w[:k], v[:, :k]
        w, v = sla.eigsh(Hm, k=k, which="SA", tol=1e-11, ncv=max(20, 3 * k + 10))
        o = np.argsort(w); return w[o], v[:, o]

    def perm_matrix(self, perm):
        """Matrix of a site permutation (sites i -> perm[i]) in this sector."""
        s = self.s; t = np.zeros_like(s)
        for i, p in enumerate(perm):
            t |= ((s >> i) & 1) << p
        j = np.searchsorted(self.s, t)
        assert np.all(self.s[j] == t)
        return sp.csr_matrix((np.ones(len(s)), (j, np.arange(len(s)))), shape=(len(s), len(s)))


def lattice_maps(g):
    """Hub-centred inversion C2 (A <-> B) and the primitive translation t (A->B->H)
    as site permutations, found from primitive triangular coordinates modulo the torus."""
    T = np.array(g["T"]); P = T @ np.array([[1, 1], [2, -1]])     # torus in primitive coords
    c = np.array(g["coords"]); Pinv = np.linalg.inv(P.T)
    def find(x):
        for i, y in enumerate(c):
            d = Pinv @ (x - y)
            if np.max(abs(d - np.rint(d))) < 1e-9:
                return i
        raise ValueError(x)
    hub0 = c[2]
    inv = [find(2 * hub0 - x) for x in c]
    shift = [find(x + np.array([1, 0])) for x in c]
    return np.array(inv), np.array(shift)


def check_symmetry(g, perm, kinds):
    E = {(tuple(sorted((a, b))), t) for a, b, t in g["bonds"]}
    Em = {(tuple(sorted((perm[a], perm[b]))), t) for a, b, t in g["bonds"]}
    if kinds == "typed":
        return E == Em
    return {e for e, _ in E} == {e for e, _ in Em}


if __name__ == "__main__":
    shapes = geometry.catalog()
    for name in ["K18", "C18"]:
        g = shapes[name]; inv, sh = lattice_maps(g)
        print(name, "inversion keeps typed bonds:", check_symmetry(g, inv, "typed"),
              "| inversion swaps A,B:", all((inv[i] % 3) == [1, 0, 2][i % 3] for i in range(g["N"])),
              "| t keeps untyped graph:", check_symmetry(g, sh, "untyped"),
              "| t cycles A->?:", sorted({(i % 3, sh[i] % 3) for i in range(g["N"])}))
