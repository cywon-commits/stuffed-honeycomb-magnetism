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

    def __init__(self, J, cluster="N27"):
        self.J = tuple(float(x) for x in J)
        self.c = Cluster(cluster, couplings=self.J)
        self.tag = hashlib.sha1((ckey(self.J) if cluster == "N27" else cluster + "|" + ckey(self.J)).encode()).hexdigest()[:12]
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


def lowest_per_k(cache, J, nup, ks=None, vectors=True, log=print, cluster="N27"):
    """Lowest level of every requested cell momentum (all momenta if ks is None).

    Returns {k: row}; rows for -k partners are copied (source='time-reversal of k').
    """
    P = Point(J, cluster)
    ck = ckey(J) if cluster == "N27" else cluster + "|" + ckey(J)
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
            key = f"p1|{ck}|{nup}|{s['i']}"
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
            key = f"p2|{ck}|{nup}|{best['sector']}"
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


# ---------------------------------------------------------------------------------------
# Second request (Tasks E-G): bond-type resolved correlations and sublattice magnetization
# ---------------------------------------------------------------------------------------
BOND_NAMES = {0: "AB", 1: "CB", 2: "CA"}   # bond type -> label (1 = hub-B = J_BC, 2 = hub-A = J_CA)


def _obs_pairs(P):
    path = P.dir / "pairs_obs.txt"
    if not path.exists():
        P.c.write_pairs(path)
    lines = path.read_text().splitlines()
    npair = int(lines[0].split()[1])
    pairs = [tuple(map(int, l.split())) for l in lines[1:1 + npair]]
    sites = [tuple(map(int, l.split())) for l in lines[2 + npair:]]
    return path, pairs, sites


def _run_obs(P, spec, ckpt, pairs_path):
    out = subprocess.run([str(x) for x in (EXE, spec, "--J", P.jarg(), "--tol", TOL, "--maxit", MAXIT,
                          "--ckpt", ckpt, "--pairs", pairs_path)], capture_output=True, text=True)
    if out.returncode != 0:
        raise RuntimeError(out.stderr[-2000:])
    d = {"orb": {}, "orbzz": {}, "siteorb": {}, "htype": {}}
    for l in out.stdout.splitlines():
        p = l.split()
        if l.startswith("ORBIT "): d["orb"][int(p[1])] = float(p[2])
        elif l.startswith("ORBITZZ "): d["orbzz"][int(p[1])] = float(p[2])
        elif l.startswith("SITEORB "): d["siteorb"][int(p[1])] = float(p[2])
        elif l.startswith("HTYPE "): d["htype"][int(p[1])] = float(p[2])
        elif l.startswith("GROUND"): d["E"] = float(p[2]); d["H_residual"] = float(p[4])
        elif l.startswith("S2ACTION"): d["S"] = float(p[2]); d["S2"] = float(p[4]); d["S2_residual"] = float(p[6])
        elif l.startswith("BASIS"): d["dim"] = int(p[2])
        elif l.startswith("PASS1"): d["iters"] = int(p[2]); d["converged"] = p[4] == "yes"
    return d


def sector_observables(P, raw, pairs, sites):
    """Bond-type zz / transverse correlations per bond (/27) and sublattice m from orbit sums."""
    osize, zsize = {}, {}
    for _, _, o, i, j in pairs:
        osize[o] = osize.get(o, 0) + 1
    oid = {(min(i, j), max(i, j)): o for _, _, o, i, j in pairs}
    for _, so, i in sites:
        zsize[so] = zsize.get(so, 0) + 1
    so_of = {i: so for _, so, i in sites}
    nb = {0: 0, 1: 0, 2: 0}
    full = {0: 0.0, 1: 0.0, 2: 0.0}; zz = {0: 0.0, 1: 0.0, 2: 0.0}
    for a, b, t in P.c.bonds:
        o = oid[(min(a, b), max(a, b))]
        full[t] += raw["orb"][o] / osize[o]; zz[t] += raw["orbzz"][o] / osize[o]; nb[t] += 1
    res = {}
    for t, name in BOND_NAMES.items():
        res[f"zz_{name}"] = zz[t] / nb[t]
        res[f"t_{name}"] = (full[t] - zz[t]) / nb[t]
        res[f"Htype_{name}"] = raw["htype"][t] / nb[t]          # independent: separate H_t action
    for s, name in ((0, "A"), (1, "B"), (2, "C")):
        idx = [i for i in range(P.c.N) if i % 3 == s]
        res[f"m_{name}"] = sum(raw["siteorb"][so_of[i]] / zsize[so_of[i]] for i in idx) / len(idx)
    res["n_bonds_per_type"] = nb[0]
    return res


OBS_KEYS = ["zz_AB", "t_AB", "zz_CA", "t_CA", "zz_CB", "t_CB", "m_A", "m_B", "m_C",
            "Htype_AB", "Htype_CA", "Htype_CB"]


def lowest_multiplet_observables(cache, J, nup, ks, log=print, degtol=1e-8, cluster="N27"):
    """Lowest level over the momenta ks (all sectors there); observables averaged over every
    sector state within degtol of it (a -k partner counts as one more state with equal values)."""
    P = Point(J, cluster)
    ck = ckey(J) if cluster == "N27" else cluster + "|" + ckey(J)
    pairs_path, pairs, sites = _obs_pairs(P)
    lv = lowest_per_k(cache, J, nup, ks=ks, vectors=False, log=log, cluster=cluster)
    E0 = min(r["E"] for r in lv.values())
    secs = P.sectors(nup)
    want = sorted({kkey(k) for k in ks}) if ks is not None else sorted({kkey(s["k"]) for s in secs})
    members = []
    for k in want:
        for s in secs:
            if kkey(s["k"]) != k or s["dim"] == 0:
                continue
            src = k if f"p1|{ck}|{nup}|{s['i']}" in cache.d else neg(k)
            # pass-1 energy of this sector (k itself, or its -k partner by time reversal)
            if src == k:
                e = cache.get(f"p1|{ck}|{nup}|{s['i']}")["E_pass1"]
                if e <= E0 + degtol:
                    members.append((k, s, "computed"))
    # -k partners whose own sectors were not computed: count via time reversal of a computed member
    counted = {m[0] for m in members}
    extra = [(neg(k), s, "time-reversal") for k, s, _ in members if neg(k) in want and neg(k) not in counted and neg(k) != k]
    out_states = []
    for k, s, how in members:
        key = f"obs|{ck}|{nup}|{s['i']}"
        r = cache.get(key)
        if r is None:
            raw = _run_obs(P, s["spec"], P.dir / f"n{nup}_s{s['i']:03d}", pairs_path)
            r = dict(sector_observables(P, raw, pairs, sites), E=raw["E"], H_residual=raw["H_residual"], S=raw["S"],
                     S2_residual=raw["S2_residual"], dim=raw.get("dim"), iters=raw.get("iters"), converged=raw.get("converged"),
                     J=list(J), nup=nup, sector=s["i"], k=list(k))
            cache.put(key, r)
            log(f"  obs J={J} nup={nup} k={k} s{s['i']} E={r['E']:.12f} S={r['S']} zzAB={r['zz_AB']:.6f} tAB={r['t_AB']:.6f} mC={r['m_C']:.6f} {r.get('H_residual')}")
        out_states.append((k, r, how))
    for k, s, how in extra:
        r = cache.get(f"obs|{ck}|{nup}|{s['i']}")
        out_states.append((k, r, how))
    n = len(out_states)
    avg = {x: sum(r[x] for _, r, _ in out_states) / n for x in OBS_KEYS}
    best = min(out_states, key=lambda x: x[1]["E"])
    return dict(avg, E=E0, multiplet_dim=n, multiplet_k=sorted({tuple(k) for k, _, _ in out_states}),
                S=best[1]["S"], H_residual=max(r["H_residual"] for _, r, _ in out_states),
                S2_residual=max(r["S2_residual"] for _, r, _ in out_states),
                converged=all(r.get("converged", True) for _, r, _ in out_states),
                dim=best[1].get("dim"), iters=best[1].get("iters"), k=list(best[0]),
                spread={x: max(r[x] for _, r, _ in out_states) - min(r[x] for _, r, _ in out_states) for x in OBS_KEYS})
