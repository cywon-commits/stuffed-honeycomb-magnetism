"""Symmetry front end for the symmetry-adapted stuffed-honeycomb ED solver (ed_sym).

Builds the cluster from the project's reference geometry module
(src/workspaces/finite_size_18_27_36pilot/geometry.py), finds the full
bond-type- and hub-preserving automorphism group of the cluster, splits the
Hilbert space into complete sets of one-dimensional symmetry sectors
(momentum k x 1D irreps of an abelian part of the little group x spin flip),
counts sector dimensions with Burnside's lemma, and writes plain-text sector
specifications for the C++ solver.

Conventions
-----------
* Spin-1/2 sites, bit = 1 means spin up.  Bits are reordered so that the
  N_R rim sites occupy bits 0..N_R-1 and the N_H hub sites occupy the top
  bits N_R..N-1 (required by the solver's hub-sublattice representative search).
* A group element is a site permutation p (site i -> p[i]) optionally combined
  with a global spin flip.  Operators compose as op(x)op(y) = op(x o y) with
  (x o y)[i] = x[y[i]].
* Sector states are P|r> with P = sum_g chi(g)^* g.  chi restricted to the
  translations is exp(2 pi i k.R) with R the translation vector in cell
  coordinates of geometry.make(); k and -k give identical spectra.
"""
import itertools, json, math, pathlib, sys
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
GEOM_DIR = HERE.parent / "workspaces" / "finite_size_18_27_36pilot"
sys.path.insert(0, str(GEOM_DIR))
import geometry as G  # noqa: E402

CLUSTERS = {
    "N12": [[2, 0], [0, 2]],
    "K18": [[3, 0], [2, 2]],
    "C18": [[1, 2], [0, 6]],
    "N27": [[3, 0], [0, 3]],
    "N36": [[2, 2], [-2, 4]],
}


def compose(x, y):
    return tuple(x[i] for i in y)


def inverse(p):
    q = [0] * len(p)
    for i, j in enumerate(p):
        q[j] = i
    return tuple(q)


class Cluster:
    def __init__(self, name):
        self.name = name
        self.T = CLUSTERS[name]
        g = G.make(self.T)
        assert g["valid"], name
        self.N = g["N"]
        self.Ncell = g["Ncell"]
        self.cells = [tuple(map(int, c)) for c in g["cells"]]
        self.momenta = [tuple(map(float, q)) for q in g["momenta"]]
        self.bonds = [(int(a), int(b), int(t)) for a, b, t in g["bonds"]]  # t=0 rim-rim (J), t=1 hub-rim (J')
        self.is_hub = [i % 3 == 2 for i in range(self.N)]
        self.trans = {tuple(map(int, g["cells"][k])): tuple(map(int, g["perms"][k])) for k in range(self.Ncell)}
        self.trans_of_perm = {p: R for R, p in self.trans.items()}
        self.group = self._automorphisms()
        # bit order: rim sites first, hub sites on top
        rims = [i for i in range(self.N) if not self.is_hub[i]]
        hubs = [i for i in range(self.N) if self.is_hub[i]]
        self.bit_of_site = {s: b for b, s in enumerate(rims + hubs)}
        self.NR, self.NH = len(rims), len(hubs)

    # ---------------- automorphism group -----------------
    def _automorphisms(self):
        N = self.N
        adj = [dict() for _ in range(N)]
        for a, b, t in self.bonds:
            adj[a][b] = t
            adj[b][a] = t
        order = [0]
        seen = {0}
        for v in order:
            for w in sorted(adj[v]):
                if w not in seen:
                    seen.add(w)
                    order.append(w)
        out = []
        m, inv = {}, {}

        def bt():
            if len(m) == N:
                out.append(tuple(m[i] for i in range(N)))
                return
            i = order[len(m)]
            for j in range(N):
                if j in inv or self.is_hub[j] != self.is_hub[i]:
                    continue
                if all(adj[j].get(m[k]) == t for k, t in adj[i].items() if k in m) and \
                   all(adj[i].get(inv[l]) == t for l, t in adj[j].items() if l in inv):
                    m[i] = j
                    inv[j] = i
                    bt()
                    del m[i], inv[j]
        bt()
        # sanity: translations are automorphisms, group closed
        S = set(out)
        assert all(p in S for p in self.trans.values())
        assert all(compose(a, b) in S for a in out[:5] for b in out)
        return out

    # ---------------- sectors -----------------
    def chi_k(self, k, R):
        return np.exp(2j * np.pi * (k[0] * R[0] + k[1] * R[1]))

    def little_group(self, k):
        L = []
        for g in self.group:
            gi = inverse(g)
            ok = True
            for R, t in self.trans.items():
                R2 = self.trans_of_perm[compose(compose(g, t), gi)]
                if abs(self.chi_k(k, R2) - self.chi_k(k, R)) > 1e-9:
                    ok = False
                    break
            if ok:
                L.append(g)
        return L

    def _is_trans(self, p):
        return p in self.trans_of_perm

    def _quot_order(self, a):
        x, n = a, 1
        while not self._is_trans(x):
            x = compose(x, a)
            n += 1
        return n

    def _closure(self, gens):
        S = set(self.trans.values())
        frontier = list(S)
        gens = list(gens) + list(self.trans.values())
        while frontier:
            nxt = []
            for x in frontier:
                for y in gens:
                    z = compose(x, y)
                    if z not in S:
                        S.add(z)
                        nxt.append(z)
            frontier = nxt
        return S

    def abelian_part(self, L):
        """Return (subgroup S of L containing T whose irreps with fixed k are all 1D, label)."""
        abel = all(self._is_trans(compose(compose(a, b), compose(inverse(a), inverse(b)))) for a in L for b in L)
        if abel:
            return sorted(L), "full little group (abelian point part)"
        a = max(L, key=self._quot_order)
        return sorted(self._closure([a])), f"cyclic point part of order {self._quot_order(a)}"

    def characters(self, S, k):
        """All 1D characters of group S (list of perms, contains T) restricted to chi_k on T."""
        Sset = set(S)
        H = set(self.trans.values())
        gens = []
        while len(H) < len(Sset):
            a = max((x for x in Sset if x not in H), key=self._quot_order)
            gens.append(a)
            H = self._closure(gens)
        Tgens = list(self.trans.items())
        cands = []
        for a in gens:
            n = self._quot_order(a)
            x = a
            for _ in range(n - 1):
                x = compose(x, a)
            base = self.chi_k(k, self.trans_of_perm[x]) ** (1.0 / n)
            cands.append([base * np.exp(2j * np.pi * m / n) for m in range(n)])
        ident = tuple(range(self.N))
        chars = []
        for choice in itertools.product(*cands):
            phi = {ident: 1.0 + 0j}
            q = [ident]
            ok = True
            stepgens = [(t, self.chi_k(k, R)) for R, t in Tgens] + list(zip(gens, choice))
            while q and ok:
                x = q.pop()
                for y, py in stepgens:
                    z = compose(x, y)
                    v = phi[x] * py
                    if z in phi:
                        if abs(phi[z] - v) > 1e-8:
                            ok = False
                            break
                    else:
                        phi[z] = v
                        q.append(z)
            if ok:
                assert len(phi) == len(Sset)
                chars.append(phi)
        return chars

    def sectors(self, nup, use_flip=None):
        """Complete list of sectors for the given number of up spins."""
        if use_flip is None:
            use_flip = (2 * nup == self.N)
        out = []
        for ik, k in enumerate(self.momenta):
            L = self.little_group(k)
            S, how = self.abelian_part(L)
            for ic, phi in enumerate(self.characters(S, k)):
                for z in ([1, -1] if use_flip else [None]):
                    out.append(dict(k=k, ik=ik, ichar=ic, z=z, chi=phi, little_order=len(L) // self.Ncell,
                                    point_part=how, group_order=len(S) * (2 if z else 1)))
        return out

    # ---------------- Burnside dimension -----------------
    def burnside_dim(self, sec, nup):
        N = self.N
        tot = 0j
        for p, c in sec["chi"].items():
            cyc = cycles(p)
            # plain element: invariant states = unions of cycles with nup total
            poly = np.zeros(N + 1)
            poly[0] = 1
            for L in cyc:
                new = poly.copy()
                new[L:] += poly[:N + 1 - L]
                poly = new
            fixed = poly[nup]
            tot += np.conj(c) * fixed
            if sec["z"] is not None:
                # flip o p : fixed iff every cycle is even (alternating pattern, 2 choices)
                fixed_f = 2 ** len(cyc) if all(L % 2 == 0 for L in cyc) else 0
                tot += np.conj(c * sec["z"]) * fixed_f
        dim = tot / sec["group_order"]
        assert abs(dim.imag) < 1e-6 and abs(dim.real - round(dim.real)) < 1e-6, dim
        return int(round(dim.real))

    # ---------------- export -----------------
    def write_spec(self, sec, nup, path):
        b = self.bit_of_site
        lines = [f"# stuffed-honeycomb ED sector spec: cluster {self.name} T={self.T}",
                 f"# k={list(sec['k'])} ichar={sec['ichar']} z={sec['z']} point_part={sec['point_part']}",
                 f"N {self.N} NH {self.NH} NUP {nup}",
                 f"BONDS {len(self.bonds)}"]
        for a, c, t in self.bonds:
            lines.append(f"{b[a]} {b[c]} {t}")
        elems = []
        for p, ch in sec["chi"].items():
            perm_bits = [0] * self.N
            for s in range(self.N):
                perm_bits[b[s]] = b[p[s]]
            elems.append((0, ch, perm_bits))
            if sec["z"] is not None:
                elems.append((1, ch * sec["z"], perm_bits))
        lines.append(f"ELEMENTS {len(elems)}")
        for f, ch, pb in elems:
            lines.append(f"{f} {ch.real:.17g} {ch.imag:.17g} " + " ".join(map(str, pb)))
        pathlib.Path(path).write_text("\n".join(lines) + "\n")


def cycles(p):
    seen = [False] * len(p)
    out = []
    for i in range(len(p)):
        if not seen[i]:
            L, j = 0, i
            while not seen[j]:
                seen[j] = True
                j = p[j]
                L += 1
            out.append(L)
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cluster", choices=CLUSTERS)
    ap.add_argument("--nup", type=int, help="number of up spins (default: Sz=0 or 1/2)")
    ap.add_argument("--outdir", default="specs")
    ap.add_argument("--no-write", action="store_true")
    a = ap.parse_args()
    c = Cluster(a.cluster)
    nup = a.nup if a.nup is not None else (c.N + 1) // 2
    secs = c.sectors(nup)
    out = pathlib.Path(a.outdir) / f"{c.name}_nup{nup}"
    out.mkdir(parents=True, exist_ok=True)
    rows = []
    tot = 0
    for i, s in enumerate(secs):
        d = c.burnside_dim(s, nup)
        tot += d
        tag = f"s{i:03d}"
        if not a.no_write:
            c.write_spec(s, nup, out / f"{tag}.spec")
        rows.append(dict(tag=tag, k=list(s["k"]), ichar=s["ichar"], z=s["z"], dim=d,
                         group_order=s["group_order"], little_point_order=s["little_order"], point_part=s["point_part"],
                         real=all(abs(v.imag) < 1e-12 for v in s["chi"].values())))
    full = math.comb(c.N, nup)
    summary = dict(cluster=c.name, T=c.T, N=c.N, NH=c.NH, nup=nup, automorphisms=len(c.group),
                   full_dim=full, sum_sector_dims=tot, complete=(tot == full), sectors=rows)
    (out / "sectors.json").write_text(json.dumps(summary, indent=1))
    print(f"{c.name}: |Aut|={len(c.group)} nup={nup} full={full} sum={tot} complete={tot == full} sectors={len(rows)}")
    for r in rows:
        print(f"  {r['tag']} k={np.round(r['k'], 4).tolist()} chi#{r['ichar']} z={r['z']} dim={r['dim']} real={r['real']} |S|={r['group_order']}")


if __name__ == "__main__":
    main()
