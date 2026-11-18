"""Unit tests U1-U6 of docs/guides/27SPIN_THREE_COUPLING_REQUEST_EN.md (section 3).

Run:  N27_WORK=/scratch/dir python3 tests.py [U1 U2 ...]
Writes data/field_dice/n27/unit_tests.json.
"""
import csv, itertools, json, math, pathlib, subprocess, sys, tempfile, time
import numpy as np
from core import Cache, Point, lowest_per_k, gs_spin_from_levels, kkey, nup_of, sz_of, EXE, ROOT, WORK, Cluster

OUT = ROOT / "data" / "field_dice" / "n27"
OUT.mkdir(parents=True, exist_ok=True)
cache = Cache(WORK / "cache_tests.jsonl") if WORK.mkdir(parents=True, exist_ok=True) is None else None


def verified_geometry(L, jp=1.0):
    """Verbatim copy of geometry() in
    src/history/hub_ground_state_crossings_verified/hub_ground_state_crossings/ed.py
    (copied because that module imports scipy at import time)."""
    def sid(x, y, s): return 3 * ((x % L) * L + y % L) + s
    coords = []; b = []
    for x in range(L):
        for y in range(L):
            coords.extend([(x + 2 * y, x - y), (x + 2 * y + 1, x - y), (x + 2 * y, x - y + 1)])
            a = sid(x, y, 0)
            b.extend((a, sid(xx, yy, 1), 1.) for xx, yy in [(x, y), (x - 1, y), (x, y - 1)])
            b.extend((sid(x, y, 2), sid(xx, yy, s), jp) for xx, yy, s in [(x, y, 0), (x, y, 1), (x + 1, y, 0), (x + 1, y - 1, 1), (x + 1, y - 1, 0), (x, y - 1, 1)])
    return np.array(coords), b


def U1():
    """geometry3(3, 1, r, r) == geometry(3, r) bond for bond; J_CA on hub->A, J_BC on hub->B."""
    L = 3
    coords_v, _ = verified_geometry(L)
    key = lambda u, v: ((u + 2 * v) % (3 * L), (u - v) % (3 * L))
    look = {key(*c): i for i, c in enumerate(coords_v)}
    res = {}
    for jab, jbc, jca in [(1, 0.6, 0.6), (1.3, 0.7, 0.2), (1, 1, 0)]:
        c = Cluster("N27", couplings=(jab, jbc, jca))
        import geometry as G
        g = G.make(c.T)
        mp = [look[key(*map(int, xy))] for xy in g["coords"]]   # my site -> verified site index
        J = (jab, jbc, jca)
        mine = sorted((min(mp[a], mp[b]), max(mp[a], mp[b]), J[t]) for a, b, t in c.bonds)
        # verified geometry with the three-coupling assignment: hub bond to s=0 -> J_CA, s=1 -> J_BC
        _, vb = verified_geometry(L, 1.0)
        ver = []
        for a, b, w in vb:
            if a % 3 == 2:
                w = jca if b % 3 == 0 else jbc
            else:
                w = jab
            ver.append((min(a, b), max(a, b), w))
        ver = sorted(ver)
        same = len(mine) == len(ver) == 81 and all(x[:2] == y[:2] and abs(x[2] - y[2]) < 1e-15 for x, y in zip(mine, ver))
        res[str(J)] = same
        if (jab, jbc, jca) == (1, 0.6, 0.6):   # also equals geometry(3, r) with jp = r
            _, vr = verified_geometry(L, 0.6)
            vr = sorted((min(a, b), max(a, b), w) for a, b, w in vr)
            res["equals geometry(3, 0.6)"] = all(x[:2] == y[:2] and abs(x[2] - y[2]) < 1e-15 for x, y in zip(mine, vr))
    return dict(details=res, pass_=all(res.values()))


def U2():
    ref = [(0.275, -10.290051743981, 3.5), (0.35, -10.472017603235, 3.5), (0.425, -10.732527324629, 2.5)]
    rows = []
    for r, E, S in ref:
        J = (1.0, r, r)
        lv = lowest_per_k(cache, J, 14, ks=[(0, 0)], vectors=True)
        row = lv[(0.0, 0.0)]
        rows.append(dict(r=r, E=row["E"], E_ref=E, dE=abs(row["E"] - E), S=row["S"], S_ref=S, H_residual=row["H_residual"], S2_residual=row["S2_residual"]))
    ok = all(x["dE"] < 1e-9 and x["S"] == x["S_ref"] and x["H_residual"] < 1e-8 for x in rows)
    return dict(rows=rows, pass_=ok)


def ref_points():
    pts = {}
    for row in csv.DictReader(open(ROOT / "data" / "field_dice" / "reference_N27_highSz.csv")):
        pts.setdefault(row["point"], (tuple(float(row[x]) for x in ("J_AB", "J_BC", "J_CA")), {}))[1][int(row["n_up"])] = float(row["E_min"])
    return pts


def U3():
    rows = []
    for name, (J, levels) in ref_points().items():
        for nup, E in sorted(levels.items()):
            lv = lowest_per_k(cache, J, nup, vectors=False)
            Emin = min(r["E"] for r in lv.values())
            rows.append(dict(point=name, J=J, nup=nup, E=Emin, E_ref=E, dE=abs(Emin - E)))
    return dict(rows=rows, max_dE=max(x["dE"] for x in rows), pass_=all(x["dE"] < 1e-9 for x in rows))


def U4():
    rows = []
    S = 0.5
    for name, (J, _) in ref_points().items():
        e = {}
        for nup in (26, 27):
            lv = lowest_per_k(cache, J, nup, vectors=False)
            e[nup] = min(r["E"] for r in lv.values())
        hsat = e[27] - e[26]
        sJ = sum(J); s2 = sum(x * x for x in J); sx = J[0] * J[1] + J[1] * J[2] + J[2] * J[0]
        pred = 3 * S * (sJ + math.sqrt(max(s2 - sx, 0.0)))
        rows.append(dict(point=name, J=J, h_sat=hsat, formula=pred, diff=abs(hsat - pred)))
    return dict(rows=rows, pass_=all(x["diff"] < 1e-9 for x in rows))


def U5():
    rows = []
    for J in [(1.0, 1.0, 0.0), (1.0, 0.5, 0.0)]:
        em = {}
        for nup in range(14, 21):
            lv = lowest_per_k(cache, J, nup, vectors=(nup == 18))
            em[sz_of(nup)] = min(r["E"] for r in lv.values())
        S, E0 = gs_spin_from_levels(em)
        rows.append(dict(J=J, Emin_by_Sz={str(k): v for k, v in em.items()}, S=S, E0=E0))
    return dict(rows=rows, pass_=all(x["S"] == 4.5 for x in rows))


def U6():
    """Hermiticity on random vectors (off-mirror N27 sector) and, on the 12-site torus with three
    unequal couplings, equality of the union of all sector spectra with the full fixed-Sz spectrum
    (checks translation covariance and completeness of the reduced operators)."""
    out = {}
    J = (1.3, 0.7, 0.2)
    P = Point(J)
    s = P.sectors(14)[1]
    o = subprocess.run([str(EXE), str(s["spec"]), "--J", P.jarg(), "--mode", "bench", "--nbench", "1"], capture_output=True, text=True, check=True).stdout
    l = next(x for x in o.splitlines() if x.startswith("HERMITICITY"))
    a, b = float(l.split("=")[1].split()[0]), float(l.split("=")[2].strip(" )"))
    out["N27_hermiticity_relative"] = a / b
    c = Cluster("N12", couplings=J)
    worst = 0.0
    for nup in (6, 7):
        states = [st for st in range(1 << c.N) if bin(st).count("1") == nup]
        idx = {st: i for i, st in enumerate(states)}
        H = np.zeros((len(states),) * 2)
        for i, st in enumerate(states):
            for p, q, t in c.bonds:
                w = J[t]
                if ((st >> p) & 1) == ((st >> q) & 1):
                    H[i, i] += 0.25 * w
                else:
                    H[i, i] -= 0.25 * w
                    H[idx[st ^ (1 << p) ^ (1 << q)], i] += 0.5 * w
        ref = np.linalg.eigvalsh(H)
        ev = []
        with tempfile.TemporaryDirectory() as td:
            for i, sec in enumerate(c.sectors(nup)):
                sp = pathlib.Path(td) / f"s{i}.spec"
                c.write_spec(sec, nup, sp)
                txt = subprocess.run([str(EXE), str(sp), "--J", ",".join(map(repr, J)), "--mode", "dense"], capture_output=True, text=True, check=True).stdout.splitlines()
                k0 = next(j for j, x in enumerate(txt) if x.startswith("DENSE"))
                D = int(txt[k0].split()[1])
                M = np.zeros((D, D), complex)
                for x in txt[k0 + 1:k0 + 1 + D * D]:
                    p, q, re_, im_ = x.split()
                    M[int(p), int(q)] = float(re_) + 1j * float(im_)
                worst = max(worst, float(np.abs(M - M.conj().T).max()) if D else 0.0)
                ev.extend(np.linalg.eigvalsh((M + M.conj().T) / 2))
        ev = np.sort(ev)
        out[f"N12_nup{nup}_spectrum_maxdiff"] = float(np.abs(ev - ref).max()) if len(ev) == len(ref) else None
        out[f"N12_nup{nup}_dims"] = [len(ev), len(ref)]
    out["N12_max_hermiticity_error"] = worst
    ok = out["N27_hermiticity_relative"] < 1e-12 and worst < 1e-12 and all(
        out[f"N12_nup{n}_spectrum_maxdiff"] is not None and out[f"N12_nup{n}_spectrum_maxdiff"] < 1e-10 for n in (6, 7))
    return dict(details=out, pass_=ok)


TESTS = dict(U1=U1, U6=U6, U2=U2, U3=U3, U4=U4, U5=U5)


def U7():
    """Observables on the 12-site three-coupling cluster, J=(1,0.6,0.6), S^z=1 (n_up=7): reduced-basis
    multiplet averages vs. dense unreduced diagonalization (trace over the degenerate ground space)."""
    from core import lowest_multiplet_observables, BOND_NAMES
    J = (1.0, 0.6, 0.6)
    c = Cluster("N12", couplings=J)
    nup = 7
    states = [st for st in range(1 << c.N) if bin(st).count("1") == nup]
    idx = {st: i for i, st in enumerate(states)}
    D = len(states)
    H = np.zeros((D, D))
    ops = {t: (np.zeros((D, D)), np.zeros((D, D))) for t in range(3)}   # (zz, transverse) per bond type
    for i, st in enumerate(states):
        for p, q, t in c.bonds:
            same = ((st >> p) & 1) == ((st >> q) & 1)
            zz = 0.25 if same else -0.25
            H[i, i] += J[t] * zz
            ops[t][0][i, i] += zz
            if not same:
                j = idx[st ^ (1 << p) ^ (1 << q)]
                H[j, i] += 0.5 * J[t]
                ops[t][1][j, i] += 0.5
    w, v = np.linalg.eigh(H)
    deg = int(np.sum(w < w[0] + 1e-8))
    V = v[:, :deg]
    nb = {t: sum(1 for *_, tt in c.bonds if tt == t) for t in range(3)}
    ref = {}
    for t, name in BOND_NAMES.items():
        ref[f"zz_{name}"] = float(np.trace(V.T @ ops[t][0] @ V)) / deg / nb[t]
        ref[f"t_{name}"] = float(np.trace(V.T @ ops[t][1] @ V)) / deg / nb[t]
    for s, name in ((0, "A"), (1, "B"), (2, "C")):
        sites = [i for i in range(c.N) if i % 3 == s]
        mz = np.array([sum(((st >> i) & 1) - 0.5 for i in sites) / len(sites) for st in states])
        ref[f"m_{name}"] = float(np.sum((V ** 2).T @ mz)) / deg
    mine = lowest_multiplet_observables(cache, J, nup, ks=None, cluster="N12")
    diffs = {k: abs(mine[k] - ref[k]) for k in ref}
    return dict(E_dense=float(w[0]), E_reduced=mine["E"], degeneracy_dense=deg, multiplet_dim=mine["multiplet_dim"],
                ref=ref, mine={k: mine[k] for k in ref}, max_diff=max(diffs.values()), diffs=diffs,
                pass_=bool(max(diffs.values()) <= 1e-10 and abs(mine["E"] - w[0]) < 1e-10 and deg == mine["multiplet_dim"]))


def U8():
    J = (1.0, 0.6, 0.6)
    lv = lowest_per_k(cache, J, 18, ks=[(0, 0)], vectors=False)
    E = lv[(0.0, 0.0)]["E"]
    return dict(E=E, E_ref=-10.395873254626, diff=abs(E + 10.395873254626), pass_=abs(E + 10.395873254626) <= 1e-10)


def U9():
    """Sum rule zz_t + t_t = <H_t>/27 for every Task E state; <H_t> from a separate H_t action."""
    tc = Cache(WORK / "cache_tasks.jsonl")
    rows = [r for k, r in tc.d.items() if k.startswith("obs|")]
    worst = 0.0
    for r in rows:
        for name in ("AB", "CA", "CB"):
            worst = max(worst, abs(r[f"zz_{name}"] + r[f"t_{name}"] - r[f"Htype_{name}"]))
    return dict(states_checked=len(rows), max_diff=worst, pass_=bool(rows) and worst <= 1e-10)


TESTS.update(U7=U7, U8=U8, U9=U9)

if __name__ == "__main__":
    want = sys.argv[1:] or list(TESTS)
    path = OUT / "unit_tests.json"
    res = json.loads(path.read_text()) if path.exists() else {}
    for name in want:
        t = time.time()
        r = TESTS[name]()
        r["wall_s"] = time.time() - t
        res[name] = r
        print(name, "PASS" if r["pass_"] else "FAIL", f"{r['wall_s']:.0f}s", json.dumps({k: v for k, v in r.items() if k not in ('rows',)})[:600], flush=True)
        path.write_text(json.dumps(res, indent=1))
