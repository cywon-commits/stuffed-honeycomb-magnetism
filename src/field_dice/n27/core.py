"""27-spin three-coupling ED driver (docs/guides/27SPIN_THREE_COUPLING_REQUEST_EN.md).

H = J_AB sum_{A-B} S.S + J_BC sum_{B-C} S.S + J_CA sum_{C-A} S.S - h sum S^z

Solver: src/ed36/ed_sym (symmetry-adapted matrix-free Lanczos, validated against the
verified 18/27-spin data) with bond types 0 = AB, 1 = hub-B (J_BC), 2 = hub-A (J_CA).
The symmetry group is the automorphism group of the coupling-labelled 27-site graph
that keeps the hub sublattice C fixed: cell translations (9) everywhere, plus point
operations only where the couplings allow them (mirror lines, J_AB = J_BC = J_CA).

For each (couplings, n_up) every sector of a complete 1D-sector set is solved in
pass 1 (eigenvalue only).  The lowest level of each cell momentum is then rebuilt
(pass 2) and its residual ||H psi - E psi|| and S^2 action are measured.  k and -k
have identical spectra (H is real), so -k rows are copied from k and marked.
Results are cached in a JSON-lines file so interrupted runs resume.
"""
import hashlib, json, os, pathlib, subprocess, sys, time

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent.parent
ED36 = ROOT / "src" / "ed36"
sys.path.insert(0, str(ED36))
from symmetry import Cluster  # noqa: E402

EXE = ED36 / "ed_sym"
WORK = pathlib.Path(os.environ.get("N27_WORK", "/tmp/n27_work"))
TOL = 1e-9          # Ritz residual target in pass 1 -> true residual ~1e-9 < 1e-8 required
MAXIT = 2000
N = 27


def kkey(k):
    return tuple(round(x % 1.0, 6) % 1.0 for x in k)


def neg(k):
    return kkey((-k[0], -k[1]))


def ckey(J):
    return ",".join(f"{x:.12f}" for x in J)


class Cache:
    def __init__(self, path):
        self.path = pathlib.Path(path)
        self.d = {}
        if self.path.exists():
            for l in self.path.read_text().splitlines():
                if l.strip():
                    r = json.loads(l)
                    self.d[r["key"]] = r

    def get(self, key):
        return self.d.get(key)

    def put(self, key, rec):
        rec = dict(rec, key=key)
        self.d[key] = rec
        with open(self.path, "a") as f:
            f.write(json.dumps(rec) + "\n")


def _run(args):
    t = time.time()
    out = subprocess.run([str(EXE)] + [str(a) for a in args], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[-2000:])
    d = {"wall_s": time.time() - t}
    for l in out.stdout.splitlines():
        p = l.split()
        if l.startswith("BASIS"):
            d["dim"] = int(p[2])
        elif l.startswith("PASS1"):
            d["iters"] = int(p[2]); d["converged"] = p[4] == "yes"; d["E_pass1"] = float(p[6]); d["ritz_resid"] = float(p[8])
        elif l.startswith("GROUND"):
            d["E"] = float(p[2]); d["H_residual"] = float(p[4])
        elif l.startswith("S2ACTION"):
            d["S"] = float(p[2]); d["S2"] = float(p[4]); d["S2_residual"] = float(p[6])
        elif l.startswith("DONE"):
            d["peak_rss_mib"] = float(p[-2])
        elif l.startswith("NOT CONVERGED"):
            d["converged"] = False
    return d


class Point:
    """One coupling point: symmetry group, sectors and pair file are built once."""

    def __init__(self, J):
        self.J = tuple(float(x) for x in J)
        self.c = Cluster("N27", couplings=self.J)
        self.tag = hashlib.sha1(ckey(self.J).encode()).hexdigest()[:12]
        self.dir = WORK / self.tag
        self.dir.mkdir(parents=True, exist_ok=True)
        self.pairs = self.dir / "pairs.txt"
        if not self.pairs.exists():
            self.c.write_pairs(self.pairs)
        self._sec = {}

    def sectors(self, nup):
        if nup not in self._sec:
            secs = self.c.sectors(nup)
            for i, s in enumerate(secs):
                s["spec"] = self.dir / f"n{nup}_s{i:03d}.spec"
                s["i"] = i
                s["dim"] = self.c.burnside_dim(s, nup)
                if not s["spec"].exists():
                    self.c.write_spec(s, nup, s["spec"])
            self._sec[nup] = secs
        return self._sec[nup]

    def jarg(self):
        return ",".join(repr(x) for x in self.J)


def lowest_per_k(cache, J, nup, ks=None, vectors=True, log=print):
    """Lowest level of every requested cell momentum (all momenta if ks is None).

    Returns {k: row}; rows for -k partners are copied (source='time-reversal of k').
    """
    P = Point(J)
    secs = P.sectors(nup)
    allk = sorted({kkey(s["k"]) for s in secs})
    want = allk if ks is None else sorted({kkey(k) for k in ks})
    solve = []
    for k in want:
        if neg(k) in solve:
            continue
        solve.append(k)
    out = {}
    for k in solve:
        best = None
        for s in secs:
            if kkey(s["k"]) != k or s["dim"] == 0:   # empty sectors carry no levels
                continue
            key = f"p1|{ckey(J)}|{nup}|{s['i']}"
            r = cache.get(key)
            if r is None:
                d = _run([s["spec"], "--J", P.jarg(), "--tol", TOL, "--maxit", MAXIT, "--novec",
                          "--ckpt", P.dir / f"n{nup}_s{s['i']:03d}"])
                r = dict(d, J=list(J), nup=nup, sector=s["i"], k=list(k), z=s["z"], ichar=s["ichar"])
                cache.put(key, r)
                log(f"  pass1 J={J} nup={nup} k={k} s{s['i']} dim={r.get('dim')} E={r.get('E_pass1')} it={r.get('iters')} conv={r.get('converged')} {r['wall_s']:.0f}s")
            if best is None or r["E_pass1"] < best["E_pass1"]:
                best = r
        if best is None:   # every sector of this momentum is empty (possible only at extreme S^z)
            continue
        row = dict(best)
        if vectors:
            key = f"p2|{ckey(J)}|{nup}|{best['sector']}"
            r2 = cache.get(key)
            if r2 is None:
                s = secs[best["sector"]]
                d = _run([s["spec"], "--J", P.jarg(), "--tol", TOL, "--maxit", MAXIT,
                          "--ckpt", P.dir / f"n{nup}_s{s['i']:03d}", "--pairs", P.pairs])
                r2 = dict(d, J=list(J), nup=nup, sector=best["sector"])
                cache.put(key, r2)
                log(f"  pass2 J={J} nup={nup} k={k} E={r2.get('E')} res={r2.get('H_residual')} S={r2.get('S')} S2res={r2.get('S2_residual')} {r2['wall_s']:.0f}s")
            row.update({x: r2.get(x) for x in ("E", "H_residual", "S", "S2", "S2_residual", "peak_rss_mib")})
        else:
            row["E"] = row["E_pass1"]
        row["source"] = "computed"
        out[k] = row
        if neg(k) in want and neg(k) != k:
            out[neg(k)] = dict(row, k=list(neg(k)), source=f"time-reversal partner of k={list(k)}")
    return out


def gs_spin_from_levels(Emin_by_sz, tol=1e-8):
    """Ground spin = largest S^z whose minimum equals the global minimum (within tol)."""
    E0 = min(Emin_by_sz.values())
    return max(sz for sz, e in Emin_by_sz.items() if e <= E0 + tol), E0


def nup_of(sz):
    return int(round(N / 2 + sz))


def sz_of(nup):
    return nup - N / 2
