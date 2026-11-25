"""36-spin Task I (docs/guides/36SPIN_의뢰서_I_KO.md): S^z ladder of the honeycomb-mirror plateau
at a third cluster size, J = (1, r, r), N = 36 torus T = [[2,2],[-2,4]] (plateau sector S^z = N/6 = 6).

Solver: src/ed36/ed_sym (symmetry-adapted Lanczos).  Sectors: cell momentum k with the 1D irreps of the
abelian part of its little group (C6 at Gamma, C3 at K), as in the solver's default.  Sectors whose spectra
are identical by symmetry are solved once:
  * time reversal: at k = -k, a character and its complex conjugate;
  * conjugation by an element of the full (hub-preserving) automorphism group that fixes k, e.g. the
    sigma_v mirrors of C3v at K, which exchange the C3 characters omega and omega-bar (2D irrep E).
The equivalence classes are computed from the characters; they are verified numerically at S^z = 8.
K' = -K is the time-reversal partner of K (identical spectrum).

Pass 1 (eigenvalue) runs for every class representative.  The lowest sector of each S^z is rebuilt (pass 2)
for ||H psi - E psi|| and S from S^2, and for S^z = 5, 6, 7 the bond-type zz / transverse correlations are
measured with the pair-orbit sums.  All solver results are cached (resumable).

Usage: N36_WORK=/scratch python3 task_i.py [U11 U13 U12 I1 I2 I3 I1b I4dims ...]
"""
import csv, json, math, os, pathlib, subprocess, sys, time
import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent.parent.parent
WORK = pathlib.Path(os.environ.get("N36_WORK", "/tmp/n36_work"))
os.environ["N27_WORK"] = str(WORK)                      # core.Point writes its specs below this directory
sys.path.insert(0, str(ROOT / "src" / "field_dice" / "n27"))
sys.path.insert(0, str(ROOT / "src" / "field_dice"))
import core                                               # noqa: E402
from core import Cache, Point, _run, _run_obs, _obs_pairs, sector_observables, kkey, neg, ckey, EXE, TOL, MAXIT  # noqa: E402
from symmetry import compose, inverse                    # noqa: E402

OUT = ROOT / "data" / "field_dice" / "n36"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
cache = Cache(WORK / "cache_n36.jsonl")
N = 36
PLATEAU = 6
GAMMA, K = (0.0, 0.0), kkey((1 / 3, 2 / 3))
DEGTOL = 1e-8          # multiplet (cross-sector) degeneracy
GAPTOL = 1e-6          # in-sector gap below which the lowest level is reported as degenerate
GHOST = 1e-9           # Ritz copies closer than this to E0 are Lanczos ghosts


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def nup_of(sz):
    return int(round(N / 2 + sz))


# --------------------------------------------------------------------------- sector equivalence classes
def classes(P, nup, ks):
    secs = [s for s in P.sectors(nup) if kkey(s["k"]) in ks and s["dim"] > 0]
    parent = {s["i"]: s["i"] for s in secs}

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    def union(a, b, why):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[max(ra, rb)] = min(ra, rb)
            reason[(min(ra, rb), max(ra, rb))] = why

    reason = {}
    close = lambda c1, c2: all(abs(c1[h] - c2.get(h, 1e9)) < 1e-9 for h in c1)
    for s in secs:
        k = kkey(s["k"])
        same = [t for t in secs if kkey(t["k"]) == k and t["i"] != s["i"]]
        if neg(k) == k:   # time reversal
            cc = {h: np.conj(v) for h, v in s["chi"].items()}
            for t in same:
                if close(cc, t["chi"]):
                    union(s["i"], t["i"], "time reversal (conjugate character)")
        for g in P.c.group:
            gi = inverse(g)
            img = {}
            ok = True
            for h in s["chi"]:
                hh = compose(compose(gi, h), g)
                if hh not in s["chi"]:
                    ok = False
                    break
                img[h] = s["chi"][hh]
            if not ok:
                continue
            for t in same:
                if close(img, t["chi"]):
                    union(s["i"], t["i"], "conjugation by a little-group element")
    out = {}
    for s in secs:
        out.setdefault(find(s["i"]), []).append(s)
    return out, secs


def char_label(s):
    """Characters of the point elements fixing hub site 2 (order, fixed sites) -> value."""
    lab = {}
    for p, c in s["chi"].items():
        if p[2] != 2:
            continue
        cyc = []
        seen = set()
        for i in range(len(p)):
            if i in seen:
                continue
            L, j = 0, i
            while j not in seen:
                seen.add(j); j = p[j]; L += 1
            cyc.append(L)
        order = math.lcm(*cyc)
        if order == 1:
            continue
        key = f"C{order}" if sum(1 for i in range(len(p)) if p[i] == i) <= 3 else f"o{order}"
        lab.setdefault(key, set()).add(complex(round(c.real, 3), round(c.imag, 3)))
    return ";".join(f"{k}:{'/'.join(str(v) for v in sorted(vs, key=lambda z: (z.real, z.imag)))}" for k, vs in sorted(lab.items()))


# --------------------------------------------------------------------------- solving
def pass1(P, s, nup):
    key = f"n36p1|{P.c.name}|{ckey(P.J)}|{nup}|{s['i']}"
    r = cache.get(key)
    if r is None:
        d = _run([s["spec"], "--J", P.jarg(), "--tol", TOL, "--maxit", MAXIT, "--novec",
                  "--ckpt", P.dir / f"n{nup}_s{s['i']:03d}"])
        r = dict(d, J=list(P.J), nup=nup, sector=s["i"], k=list(kkey(s["k"])), ichar=s["ichar"])
        cache.put(key, r)
        log(f"  pass1 J={P.J} Sz={nup - N // 2} k={kkey(s['k'])} s{s['i']} dim={r.get('dim')} E={r.get('E_pass1')} "
            f"it={r.get('iters')} conv={r.get('converged')} {r['wall_s']:.0f}s rss={r.get('peak_rss_mib')}MiB")
    return r


def pass2_obs(P, s, nup):
    key = f"n36obs|{P.c.name}|{ckey(P.J)}|{nup}|{s['i']}"
    r = cache.get(key)
    if r is None:
        pairs_path, pairs, sites = _obs_pairs(P)
        t = time.time()
        raw = _run_obs(P, s["spec"], P.dir / f"n{nup}_s{s['i']:03d}", pairs_path)
        r = dict(sector_observables(P, raw, pairs, sites), E=raw["E"], H_residual=raw["H_residual"], S=raw["S"],
                 S2_residual=raw["S2_residual"], dim=raw.get("dim"), J=list(P.J), nup=nup, sector=s["i"], wall_s=time.time() - t)
        cache.put(key, r)
        log(f"  pass2+obs J={P.J} Sz={nup - N // 2} s{s['i']} E={r['E']:.12f} res={r['H_residual']:.2e} S={r['S']} "
            f"zzAB={r['zz_AB']:.6f} tAB={r['t_AB']:.6f} {r['wall_s']:.0f}s")
    return r


def in_sector_gap(r):
    E0 = r["E_pass1"]
    higher = [x for x in r.get("ritz", []) if x > E0 + GHOST]
    return (higher[0] - E0) if higher else None


def ladder_level(J, sz, ks=(GAMMA, K), verify=False, observables=False):
    """Lowest level of S^z over the sectors at ks.  Returns the summary row and the per-sector rows."""
    P = Point(J, "N36")
    nup = nup_of(sz)
    cls, secs = classes(P, nup, set(ks))
    rows = []
    for rep, members in sorted(cls.items()):
        srep = next(s for s in members if s["i"] == rep)
        r = pass1(P, srep, nup)
        for m in members:
            rm = pass1(P, m, nup) if (verify and m["i"] != rep) else r
            src = "computed" if (m["i"] == rep or verify) else f"equivalent to s{rep}"
            rows.append(dict(sector=m["i"], k=kkey(m["k"]), ichar=m["ichar"], label=char_label(m), dim=m["dim"],
                             E=rm["E_pass1"], ritz_resid=rm.get("ritz_resid"), iters=rm.get("iters"), converged=rm.get("converged"),
                             gap=in_sector_gap(rm), source=src, rep=rep, wall_s=rm.get("wall_s"), rss=rm.get("peak_rss_mib")))
    if verify:
        worst = max(abs(x["E"] - next(y["E"] for y in rows if y["sector"] == x["rep"])) for x in rows)
        log(f"VERIFY equivalence classes J={J} Sz={sz}: max |E(member) - E(representative)| = {worst:.2e}")
    E0 = min(x["E"] for x in rows)
    deg = [x for x in rows if x["E"] <= E0 + DEGTOL]
    mult = sum(2 if x["k"] == K else 1 for x in deg)          # K' partner (time reversal) has the same level
    best = min(rows, key=lambda x: x["E"])
    srep = next(s for s in secs if s["i"] == best["rep"])
    o = pass2_obs(P, srep, nup)
    summ = dict(r=J[1], Sz=sz, E=E0, k=best["k"], sector=best["sector"], label=best["label"], multiplet_dim=mult,
                multiplet=[(x["k"], x["ichar"]) for x in deg], H_residual=o["H_residual"], S=o["S"], S2_residual=o["S2_residual"],
                gap=best["gap"], obs=o, converged=all(x["converged"] for x in rows))
    if best["k"] != GAMMA:
        log(f"NOTE r={J[1]} Sz={sz}: lowest level not at Gamma but at {best['k']}")
    return summ, rows


# --------------------------------------------------------------------------- output
COLS = ["task", "r", "Sz", "n_up", "k1", "k2", "irrep_index", "irrep_characters", "dim", "E", "residual", "residual_kind",
        "iterations", "converged", "in_sector_gap", "source", "S", "S2_residual",
        "zz_AB", "t_AB", "zz_CA", "t_CA", "zz_CB", "t_CB", "m_A", "m_B", "m_C", "wall_s", "peak_rss_mib"]


def _rewrite(path, cols, keyf, newrows):
    rows = list(csv.DictReader(open(path))) if path.exists() else []
    nk = {keyf(r) for r in newrows}
    rows = [r for r in rows if keyf(r) not in nk] + newrows
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader(); w.writerows(rows)


def write_rows(task, J, sz, summ, rows):
    out = []
    for x in rows:
        lowest = x["sector"] == summ["sector"]
        o = summ["obs"] if lowest else {}
        out.append({c: "" for c in COLS} | dict(
            task=task, r=J[1], Sz=sz, n_up=nup_of(sz), k1=f"{x['k'][0]:.6f}", k2=f"{x['k'][1]:.6f}", irrep_index=x["ichar"],
            irrep_characters=x["label"], dim=x["dim"], E=f"{x['E']:.12f}",
            residual=(o.get("H_residual") if lowest else x["ritz_resid"]),
            residual_kind=("||H psi - E psi|| (eigenvector)" if lowest else "Lanczos Ritz estimate"),
            iterations=x["iters"], converged=x["converged"], in_sector_gap=("" if x["gap"] is None else f"{x['gap']:.3e}"),
            source=x["source"], S=o.get("S", ""), S2_residual=o.get("S2_residual", ""),
            **{c: (f"{o[c]:.12f}" if c in o else "") for c in ("zz_AB", "t_AB", "zz_CA", "t_CA", "zz_CB", "t_CB", "m_A", "m_B", "m_C")},
            wall_s=f"{x['wall_s']:.0f}" if x["wall_s"] else "", peak_rss_mib=x["rss"] or ""))
    _rewrite(OUT / "task_I.csv", COLS, lambda r: (r["task"], str(r["r"]), str(r["Sz"]), r["k1"], r["k2"], str(r["irrep_index"])), out)


SUMCOLS = ["task", "r", "quantity", "value", "note"]


def summary(task, r, q, v, note=""):
    _rewrite(OUT / "summary_n36.csv", SUMCOLS, lambda x: (x["task"], str(x["r"]), x["quantity"]),
             [dict(task=task, r=r, quantity=q, value=v, note=note)])


# --------------------------------------------------------------------------- rotor analysis
def rotor_fit(E):
    """E: {Sz: E}.  With S^z = 4 use the symmetric gaps; otherwise solve E_m - E_6 = e_m(x)/I - h0 m, m = -1, 1, 2."""
    from rotor_aniso import levels
    from scipy.optimize import brentq
    e = lambda x, m: levels(x, m)[0] - levels(x, 0)[0]
    p = PLATEAU
    d1 = (E[p + 1] + E[p - 1]) / 2 - E[p]
    h0 = (E[p - 1] - E[p + 1]) / 2
    if (p - 2) in E:
        d2 = (E[p + 2] + E[p - 2]) / 2 - E[p]; how = "symmetric gaps (S^z = 4..8)"
    else:
        d2 = E[p + 2] - E[p] + 2 * h0; how = "E_m - E_6 = e_m(x)/I - h0 m with m = -1, 1, 2"
    R = d2 / d1
    g = lambda x: e(x, 2) / e(x, 1) - R
    x = brentq(g, 0.0, 200.0) if g(0.0) * g(200.0) < 0 else float("nan")
    I = e(x, 1) / d1 if x == x else float("nan")
    return dict(W=2 * d1, Delta1=d1, Delta2=d2, ratio=R, h0=h0, x=x, I=I, I_over_N9=I / (N / 9), D_over_N=(x / I) / N if x == x else float("nan"), how=how)


def task_I(task, r, szs=(6, 7, 8, 5), verify_sz=8):
    J = (1.0, r, r)
    E = {}
    for sz in szs:
        summ, rows = ladder_level(J, sz, verify=(sz == verify_sz))
        write_rows(task, J, sz, summ, rows)
        E[sz] = summ["E"]
        o = summ["obs"]
        gap = "n/a" if summ["gap"] is None else f"{summ['gap']:.2e}"
        summary(task, r, f"E_min(Sz={sz})", f"{summ['E']:.12f}",
                f"k={summ['k']} irrep={summ['label']} multiplet_dim={summ['multiplet_dim']} S={summ['S']} "
                f"residual={summ['H_residual']:.2e} in-sector gap={gap}")
        if sz in (5, 6, 7):
            R = o["zz_AB"] / o["t_AB"]
            summary(task, r, f"R=zz_AB/t_AB(Sz={sz})", R, f"Delta=R-1/2={R - 0.5:+.6f}")
        log(f"{task} r={r} Sz={sz}: E={summ['E']:.12f} k={summ['k']} S={summ['S']} mult={summ['multiplet_dim']}")
    if all(m in E for m in (5, 6, 7, 8)):
        f = rotor_fit(E)
        for q in ("W", "Delta1", "Delta2", "ratio", "h0", "x", "I", "I_over_N9", "D_over_N"):
            summary(task, r, q if q != "ratio" else "Delta2/Delta1", f[q], f["how"] if q in ("Delta2", "x") else "")
        log(f"{task} r={r}: W36={f['W']:.4f} Delta2/Delta1={f['ratio']:.3f} x={f['x']:.3f} I/(N/9)={f['I_over_N9']:.3f} D/N={f['D_over_N']:.4f}")


# --------------------------------------------------------------------------- unit tests
def save_test(name, res):
    p = OUT / "unit_tests.json"
    d = json.loads(p.read_text()) if p.exists() else {}
    d[name] = res
    p.write_text(json.dumps(d, indent=1, default=str))
    log(f"{name} {'PASS' if res['pass_'] else 'FAIL'} {json.dumps({k: v for k, v in res.items() if k != 'rows'}, default=str)[:400]}")


def U11():
    """N27, r = 0.6, S^z = 9/2, Gamma: -10.395873254626 (ed_sym, sectors solved here)."""
    P = Point((1.0, 0.6, 0.6), "N27")
    best = min(pass1(P, s, 18)["E_pass1"] for s in P.sectors(18) if kkey(s["k"]) == GAMMA and s["dim"] > 0)
    save_test("U11", dict(E=best, E_ref=-10.395873254626, diff=abs(best + 10.395873254626), pass_=abs(best + 10.395873254626) <= 1e-9))


def U12():
    t = time.time()
    out = subprocess.run([sys.executable, "-u", str(ROOT / "src" / "ed36" / "validate_reference.py")],
                         capture_output=True, text=True, cwd=str(ROOT / "src" / "ed36"))
    res = [json.loads(l) for l in out.stdout.splitlines() if l.startswith("{")]
    save_test("U12", dict(rows=res, n=len(res), all_pass=all(x.get("pass") for x in res), wall_s=time.time() - t,
                          pass_=out.returncode == 0 and len(res) == 10 and all(x.get("pass") for x in res)))


def U13():
    """N36, r = 0.8: E(S^z = 17, Gamma lowest) - E(S^z = 18) = -h_sat = -1.5 (2 + r) = -4.2."""
    J = (1.0, 0.8, 0.8)
    P = Point(J, "N36")
    e = {}
    for nup in (35, 36):
        e[nup] = min(pass1(P, s, nup)["E_pass1"] for s in P.sectors(nup) if kkey(s["k"]) == GAMMA and s["dim"] > 0)
    diff = e[35] - e[36]
    save_test("U13", dict(E17=e[35], E18=e[36], diff=diff, expected=-4.2, err=abs(diff + 4.2), pass_=abs(diff + 4.2) <= 1e-10))


def U14():
    """Pair-orbit sum rule zz_AB + t_AB = <H_AB>/N_AB for every N36 observable state (and for H_CA + H_CB)."""
    rows = [r for k, r in cache.d.items() if k.startswith("n36obs|")]
    wAB = max((abs(r["zz_AB"] + r["t_AB"] - r["Htype_AB"]) for r in rows), default=None)
    wH = max((abs(r["zz_CA"] + r["t_CA"] + r["zz_CB"] + r["t_CB"] - r["Htype_CA"] - r["Htype_CB"]) for r in rows), default=None)
    save_test("U14", dict(states=len(rows), max_diff_AB=wAB, max_diff_CA_plus_CB=wH,
                          pass_=bool(rows) and max(wAB, wH) <= 1e-10))


def I4dims():
    """Task I4 feasibility: Burnside dimensions at J = (1, 0.15 + d/2, 0.15 - d/2) (off the mirror line)."""
    from symmetry import Cluster
    J = (1.0, 0.16, 0.14)
    c = Cluster("N36", couplings=J)
    lines = [f"|G| = {len(c.group)} (hub-preserving automorphisms of the coupling-labelled graph)"]
    for sz in range(0, 7):
        secs = c.sectors(18 + sz)
        dims = [c.burnside_dim(s, 18 + sz) for s in secs]
        lines.append(f"Sz={sz}: {len(secs)} sectors, total {sum(dims)}, largest {max(dims)}, smallest {min(dims)}")
        summary("I4", 0.15, f"dims(Sz={sz})", max(dims), f"{len(secs)} sectors, total {sum(dims)}, |G|={len(c.group)}")
    log("I4 dims: " + " | ".join(lines))


TASKS = dict(U11=U11, U12=U12, U13=U13, U14=U14, I4dims=I4dims,
             I1=lambda: task_I("I1", 0.8), I2=lambda: task_I("I2", 0.6), I3=lambda: task_I("I3", 0.9),
             I1b=lambda: task_I("I1", 0.8, szs=(6, 7, 8, 5, 4)))

if __name__ == "__main__":
    for t in sys.argv[1:]:
        t0 = time.time()
        log(f"=== {t} start")
        TASKS[t]()
        log(f"=== {t} done in {time.time() - t0:.0f}s")
