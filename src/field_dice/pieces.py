"""Bipartite pieces of the sqrt3 x sqrt3 triangular lattice: choose k sublattice pairs (AB, BC, CA) and
d bond directions; the k x d grid gives dimers, trimers, chains, zigzag chains, honeycomb, Lieb, square,
dice (and the frustrated triangle for k = d = 3). For a uniform cover each bond lies in
m = C(2,k-1) C(2,d-1) pieces, so H_tri = (1/m) sum_P H_P and E_tri >= (1/m) sum_P E_P.
Resonance efficiency f_P = E_tri / [(1/m) sum_P E_P]; classically 1/2 for every bipartite piece."""
import itertools, csv, sys, pathlib
from math import comb
import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from ed_sz import geometry

PAIR = {frozenset((0, 1)): "AB", frozenset((1, 2)): "BC", frozenset((0, 2)): "CA"}
DIRS = {(1, 0): 1, (-1, 0): 1, (0, 1): 2, (0, -1): 2, (-1, 1): 3, (1, -1): 3}


def labelled_bonds(g):
    P = np.array(g["T"]) @ np.array([[1, 1], [2, -1]])        # torus periods in primitive coordinates
    out = []
    for a, b, _ in g["bonds"]:
        d0 = g["coords"][b] - g["coords"][a]
        best = None
        for n1 in range(-7, 8):
            for n2 in range(-7, 8):
                d = tuple(int(v) for v in d0 + n1 * P[0] + n2 * P[1])
                if d in DIRS:
                    best = d
        assert best is not None, (a, b, d0)
        out.append((a, b, PAIR[frozenset((a % 3, b % 3))], DIRS[best]))
    return out


def bipartite(N, bonds):
    col = [-1] * N; adj = [[] for _ in range(N)]
    for a, b in bonds:
        adj[a].append(b); adj[b].append(a)
    for s in range(N):
        if col[s] < 0:
            col[s] = 0; st = [s]
            while st:
                u = st.pop()
                for v in adj[u]:
                    if col[v] < 0:
                        col[v] = 1 - col[u]; st.append(v)
                    elif col[v] == col[u]:
                        return False
    return True


def ground_energy(N, bonds):
    states = np.array([s for s in range(1 << N) if bin(s).count("1") == N // 2], dtype=np.int64)
    idx = {int(s): i for i, s in enumerate(states)}
    rows, cols, vals = [], [], []
    diag = np.zeros(len(states))
    for a, b in bonds:
        ba = (states >> a) & 1; bb = (states >> b) & 1
        diag += np.where(ba == bb, 0.25, -0.25)
        flip = np.nonzero(ba != bb)[0]
        new = states[flip] ^ ((1 << a) | (1 << b))
        rows += [idx[int(t)] for t in new]; cols += list(flip); vals += [0.5] * len(flip)
    H = sp.csr_matrix((vals, (rows, cols)), shape=(len(states),) * 2) + sp.diags(diag)
    return eigsh(H, k=1, which="SA", tol=1e-12)[0][0]


NAMES = {(1, 1): "dimers", (1, 2): "zigzag chains", (1, 3): "honeycomb", (2, 1): "trimers",
         (2, 2): "Lieb", (2, 3): "dice", (3, 1): "chains", (3, 2): "square", (3, 3): "triangle"}

if __name__ == "__main__":
    rows = []
    for name in ("K18", "C18"):
        g = geometry.catalog()[name]; N = g["N"]; LB = labelled_bonds(g)
        Etri = ground_energy(N, [(a, b) for a, b, _, _ in LB])
        for k in (1, 2, 3):
            for d in (1, 2, 3):
                m = comb(2, k - 1) * comb(2, d - 1)
                Es, bip, nb = [], True, []
                for ps in itertools.combinations(("AB", "BC", "CA"), k):
                    for ds in itertools.combinations((1, 2, 3), d):
                        bl = [(a, b) for a, b, p, mu in LB if p in ps and mu in ds]
                        bip &= bipartite(N, bl); nb.append(len(bl))
                        Es.append(ground_energy(N, bl) if (k, d) != (3, 3) else Etri)
                bound = sum(Es) / m
                g_P = np.mean([e / (-0.25 * n) for e, n in zip(Es, nb)]) - 1
                f = Etri / bound
                rows.append((name, NAMES[k, d], k, d, len(Es), m, nb[0], bip, Etri, bound, f, g_P))
                print(f"{name} {NAMES[k,d]:14s} pieces={len(Es)} m={m} bonds/piece={nb[0]:2d} bipartite={bip} "
                      f"E_P(mean)={np.mean(Es):9.5f} bound={bound:9.5f} f={f:.4f} g_P={g_P:.3f}", flush=True)
    out = pathlib.Path(__file__).resolve().parents[2] / "data" / "field_dice" / "pieces_18.csv"
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["shape", "piece", "k_pairs", "d_dirs", "n_pieces", "m", "bonds_per_piece",
                                        "bipartite", "E_tri", "bound", "f", "g_P"]); w.writerows(rows)
