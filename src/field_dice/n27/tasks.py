"""Production tasks A-D of docs/guides/27SPIN_THREE_COUPLING_REQUEST_EN.md (section 4).

Run:  N27_WORK=/scratch/dir python3 tasks.py A B D C
Each task appends rows to data/field_dice/n27/task_<X>.csv and derived quantities to
summary_n27.csv.  All solver results are cached (N27_WORK/cache_tasks.jsonl), so a rerun
after an interruption only computes what is missing.
"""
import cmath, csv, json, math, pathlib, sys, time
from core import Cache, lowest_per_k, gs_spin_from_levels, kkey, nup_of, sz_of, WORK, ROOT

OUT = ROOT / "data" / "field_dice" / "n27"
OUT.mkdir(parents=True, exist_ok=True)
WORK.mkdir(parents=True, exist_ok=True)
cache = Cache(WORK / "cache_tasks.jsonl")
S_SPIN = 0.5
ALLK = None
GK = [(0.0, 0.0), (1 / 3, 2 / 3), (2 / 3, 1 / 3)]
COLS = ["task", "J_AB", "J_BC", "J_CA", "eps", "phi_deg", "n_up", "Sz", "k1", "k2", "E", "S", "H_residual", "S2_residual",
        "converged", "lanczos_iterations", "sector_dim", "source"]


def eps_phi(J):
    jb = sum(J) / 3
    z = sum((J[k] / jb - 1) * cmath.exp(2j * math.pi * k / 3) for k in range(3)) / 3
    return abs(z), (math.degrees(cmath.phase(z)) % 360 if abs(z) > 1e-12 else 0.0)


def couplings(eps, phi_deg):
    return tuple(1 + 2 * eps * math.cos(math.radians(phi_deg) - 2 * math.pi * k / 3) for k in range(3))


def log(msg):
    print(time.strftime("%H:%M:%S"), msg, flush=True)


def _rewrite(path, cols, keyf, newrows):
    rows = list(csv.DictReader(open(path))) if path.exists() else []
    newkeys = {keyf(r) for r in newrows}
    rows = [r for r in rows if keyf(r) not in newkeys] + newrows
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def write_rows(task, J, nup, lv):
    """Replace (not append) the rows of this (task, J, n_up) so reruns stay duplicate-free."""
    eps, phi = eps_phi(J)
    rows = []
    for k, r in sorted(lv.items()):
        rows.append(dict(zip(COLS, [task, repr(J[0]), repr(J[1]), repr(J[2]), f"{eps:.12f}", f"{phi:.6f}", str(nup), str(sz_of(nup)),
                                    f"{k[0]:.6f}", f"{k[1]:.6f}", f"{r['E']:.12f}", str(r.get("S", "")), str(r.get("H_residual", "")),
                                    str(r.get("S2_residual", "")), str(r.get("converged")), str(r.get("iters")), str(r.get("dim")), str(r.get("source"))])))
    _rewrite(OUT / f"task_{task}.csv", COLS, lambda r: (r["task"], r["J_AB"], r["J_BC"], r["J_CA"], r["n_up"], r["k1"], r["k2"]), rows)


SUMCOLS = ["task", "point", "J_AB", "J_BC", "J_CA", "eps", "phi_deg", "quantity", "value", "note"]


def summary(rec):
    row = {k: str(rec.get(k, "")) for k in SUMCOLS}
    _rewrite(OUT / "summary_n27.csv", SUMCOLS, lambda r: (r["task"], r["point"], r["quantity"]), [row])


def emin(lv):
    good = [r for r in lv.values() if r.get("converged", True)]
    assert len(good) == len(lv), "unconverged level would enter a derived quantity"
    return min(r["E"] for r in good)


def maxwell(E, m0=4.5):
    lo = [(E[m0] - E[m]) / (m0 - m) for m in E if m < m0]
    hi = [(E[m] - E[m0]) / (m - m0) for m in E if m > m0]
    return (max(lo) if lo else None), (min(hi) if hi else None)


# ------------------------------------------------------------------ Task A
def task_A():
    restrict = None
    for r in (0.4, 0.6, 0.8):
        J = (1.0, r, r)
        E = {}
        for nup in range(14, 21):
            ks = ALLK
            if nup in (14, 15) and restrict:
                ks = GK
            lv = lowest_per_k(cache, J, nup, ks=ks, vectors=True, log=log)
            write_rows("A", J, nup, lv)
            E[sz_of(nup)] = emin(lv)
            if r == 0.4 and nup in (14, 15):
                best = min(lv.items(), key=lambda kv: kv[1]["E"])[0]
                ok = best in {kkey(k) for k in GK}
                restrict = ok if restrict is None else (restrict and ok)
                log(f"A r=0.4 nup={nup}: global minimum at k={best}; restrict later S^z=1/2,3/2 to Gamma/K/K': {restrict}")
        hlo, hhi = maxwell(E)
        Sgs, E0 = gs_spin_from_levels(E)
        eps, phi = eps_phi(J)
        band = [(1 / chi - 3 * r) / 9 for chi in (0.16, 0.14)]
        base = dict(task="A", point=f"r={r}", J_AB=J[0], J_BC=J[1], J_CA=J[2], eps=eps, phi_deg=phi)
        summary(dict(base, quantity="h_lo", value=hlo, note="Maxwell, S^z=1/2..13/2"))
        summary(dict(base, quantity="h_hi", value=hhi, note="Maxwell, S^z=1/2..13/2"))
        summary(dict(base, quantity="W_27", value=hhi - hlo, note=f"finite-size rotor-hub band {band[0]:.3f}-{band[1]:.3f}"))
        summary(dict(base, quantity="ground_spin_S", value=Sgs, note=f"E0={E0:.12f}"))
        log(f"A r={r}: h_lo={hlo:.6f} h_hi={hhi:.6f} W={hhi - hlo:.6f} band={band} S_gs={Sgs}")


# ------------------------------------------------------------------ Task B
def task_B():
    pts = [("eps0.1_phi30", couplings(0.1, 30)), ("eps0.2_phi30", couplings(0.2, 30)),
           ("dice", (1.0, 1.0, 0.0)), ("mirror_1.1_1.1_0.8", (1.1, 1.1, 0.8))]
    for name, J in pts:
        E = {}
        for nup in (17, 18, 19):
            lv = lowest_per_k(cache, J, nup, ks=ALLK, vectors=True, log=log)
            write_rows("B", J, nup, lv)
            E[sz_of(nup)] = emin(lv)
        hlo, hhi = maxwell(E)
        Jmin = min(J)
        o = sorted(J)
        rt = math.sqrt(max((o[1] - Jmin) * (o[2] - Jmin), 0.0))
        hm, hp = 3 * S_SPIN * (Jmin - rt), 3 * S_SPIN * (Jmin + rt)
        eps, phi = eps_phi(J)
        base = dict(task="B", point=name, J_AB=J[0], J_BC=J[1], J_CA=J[2], eps=eps, phi_deg=phi)
        summary(dict(base, quantity="h_lo", value=hlo, note="from S^z=7/2,9/2"))
        summary(dict(base, quantity="h_hi", value=hhi, note="from S^z=9/2,11/2"))
        summary(dict(base, quantity="h_minus_classical", value=hm))
        summary(dict(base, quantity="h_plus_classical", value=hp))
        summary(dict(base, quantity="h_hi/h_plus", value=(hhi / hp) if hp else None, note="N=18: 1.51-1.52 off-mirror"))
        summary(dict(base, quantity="h_lo-h_minus", value=hlo - hm, note="N=18: -0.08...-0.13"))
        log(f"B {name} J={J}: h_lo={hlo:.6f} h_hi={hhi:.6f} h-={hm:.6f} h+={hp:.6f} ratio={hhi / hp if hp else float('nan'):.4f}")


# ------------------------------------------------------------------ Task D
def task_D():
    Jhex = (1.0, 0.0, 0.0)
    lv = lowest_per_k(cache, Jhex, 18, ks=ALLK, vectors=True, log=log)
    write_rows("D", Jhex, 18, lv)
    Ehex = emin(lv)
    summary(dict(task="D", point="hexagon J=(1,0,0)", J_AB=1, J_BC=0, J_CA=0, quantity="E_hex", value=Ehex, note="S^z=9/2, all momenta"))
    for name, J in [("dice", (1.0, 1.0, 0.0)), ("eps0.560_phi30", couplings(0.560, 30))]:
        lv = lowest_per_k(cache, J, 18, ks=[(0, 0)], vectors=True, log=log)
        write_rows("D", J, 18, lv)
        E0 = lv[(0.0, 0.0)]["E"]
        f = E0 / (sum(J) * Ehex)
        eps, phi = eps_phi(J)
        base = dict(task="D", point=name, J_AB=J[0], J_BC=J[1], J_CA=J[2], eps=eps, phi_deg=phi)
        summary(dict(base, quantity="E0_Sz9/2_Gamma", value=E0))
        summary(dict(base, quantity="f", value=f, note="N=18: 0.82 (dice), 0.83 (phi=30 edge); LSWT 0.89"))
        log(f"D {name}: E0={E0:.12f} Ehex={Ehex:.12f} f={f:.6f}")


# ------------------------------------------------------------------ Task C
def ground_spin_jump(J, full=False):
    """True if the ground spin exceeds 1/2: E_min(S^z=3/2) == E_min(S^z=1/2) within 1e-8.
    S^z=1/2 restricted to Gamma, K, K' unless full=True."""
    lv1 = lowest_per_k(cache, J, 14, ks=(ALLK if full else GK), vectors=full, log=log)
    lv3 = lowest_per_k(cache, J, 15, ks=ALLK, vectors=full, log=log)
    e1, e3 = emin(lv1), emin(lv3)
    return e3 <= e1 + 1e-8, e1, e3, lv1, lv3


def task_C(js=(0.2, 0.15, 0.3)):
    for j in js:
        Jd = lambda d: (1.0, j + d / 2, j - d / 2)
        lo, hi = 0.0, 0.08
        jl, *_ = ground_spin_jump(Jd(lo))
        jh, *_ = ground_spin_jump(Jd(hi))
        log(f"C j={j}: S>1/2 at d=0? {jl}; at d=0.08? {jh}")
        base = dict(task="C", point=f"j={j}")
        if jl or not jh:
            summary(dict(base, quantity="d_c", value="", note=f"no single jump bracketed in [0,0.08]: S>1/2 at d=0: {jl}, at d=0.08: {jh}"))
            continue
        while hi - lo > 0.002 + 1e-12:
            mid = round((lo + hi) / 2, 10)
            jm, e1, e3, *_ = ground_spin_jump(Jd(mid))
            log(f"C j={j} d={mid:.6f}: E1/2={e1:.10f} E3/2={e3:.10f} jump={jm}")
            if jm:
                hi = mid
            else:
                lo = mid
        # final bracket and d=0: all momenta, vectors, S^z = 1/2 ... 9/2
        res = {}
        for d in (0.0, lo, hi):
            J = Jd(d)
            E = {}
            for nup in range(14, 19):
                lv = lowest_per_k(cache, J, nup, ks=ALLK, vectors=True, log=log)
                write_rows("C", J, nup, lv)
                E[sz_of(nup)] = emin(lv)
            res[d] = gs_spin_from_levels(E)
        summary(dict(base, J_AB=1, quantity="d_c", value=hi,
                     note=f"bracket [{lo:.6f},{hi:.6f}]; S(d=0)={res[0.0][0]}, S(lo)={res[lo][0]}, S(hi)={res[hi][0]}; "
                          f"d_c/j={hi / j:.4f}; predictions: finite-size 0.14j={0.14 * j:.4f}, bulk 0.13j^2={0.13 * j * j:.4f}"))
        summary(dict(base, quantity="S_after_jump", value=res[hi][0]))
        log(f"C j={j}: d_c in ({lo:.6f},{hi:.6f}], S after jump {res[hi][0]}, d_c/j={hi / j:.4f}")


TASKS = dict(A=task_A, B=task_B, D=task_D, C=task_C)

if __name__ == "__main__":
    for t in sys.argv[1:] or ["A", "B", "D", "C"]:
        t0 = time.time()
        log(f"=== task {t} start")
        TASKS[t]()
        log(f"=== task {t} done in {time.time() - t0:.0f}s")
