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
    phi = math.degrees(cmath.phase(z)) % 360 if abs(z) > 1e-12 else 0.0
    return abs(z), (0.0 if abs(phi - 360) < 1e-6 else phi)


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
                          f"d_c/j={hi / j:.4f}; predictions: finite-size 0.14j={0.14 * j:.4f}, bulk (Seifert-Vojta, corrected) 0.34-0.37j^2={0.34 * j * j:.4f}-{0.37 * j * j:.4f}"))
        summary(dict(base, quantity="S_after_jump", value=res[hi][0]))
        log(f"C j={j}: d_c in ({lo:.6f},{hi:.6f}], S after jump {res[hi][0]}, d_c/j={hi / j:.4f}")


TASKS = dict(A=task_A, B=task_B, D=task_D, C=task_C)


# ======================================================================= second request (E-G)
from core import lowest_multiplet_observables, OBS_KEYS   # noqa: E402

ECOLS = COLS + ["zz_AB", "t_AB", "zz_CA", "t_CA", "zz_CB", "t_CB", "m_A", "m_B", "m_C", "multiplet_dim", "multiplet_k"]


def ref18_anisotropy():
    out = {}
    p = ROOT / "data" / "field_dice" / "rim_anisotropy_18.csv"
    if p.exists():
        for r in csv.DictReader(open(p)):
            out.setdefault(r.get("r") or r.get("Jp") or "", []).append(r)
    return out


def task_E():
    for r in (0.6, 0.8, 0.4, 0.2, 1.0):          # decisive points first
        J = (1.0, r, r)
        ks = [(0, 0)] if r in (0.4, 0.6, 0.8) else GK
        for nup in (17, 18, 19):
            o = lowest_multiplet_observables(cache, J, nup, ks=ks, log=log)
            assert o["converged"], "unconverged level"
            eps, phi = eps_phi(J)
            row = dict(zip(COLS, ["E", repr(J[0]), repr(J[1]), repr(J[2]), f"{eps:.12f}", f"{phi:.6f}", str(nup), str(sz_of(nup)),
                                  f"{o['k'][0]:.6f}", f"{o['k'][1]:.6f}", f"{o['E']:.12f}", str(o["S"]), str(o["H_residual"]),
                                  str(o["S2_residual"]), str(o["converged"]), str(o["iters"]), str(o["dim"]), "computed (multiplet average)"]))
            row.update({x: f"{o[x]:.12f}" for x in ("zz_AB", "t_AB", "zz_CA", "t_CA", "zz_CB", "t_CB", "m_A", "m_B", "m_C")})
            row["multiplet_dim"] = str(o["multiplet_dim"])
            row["multiplet_k"] = ";".join(f"({a:.4f},{b:.4f})" for a, b in o["multiplet_k"])
            _rewrite(OUT / "task_E.csv", ECOLS, lambda x: (x["J_BC"], x["n_up"]), [row])
            R = o["zz_AB"] / o["t_AB"]
            base = dict(task="E", point=f"r={r},Sz={sz_of(nup)}", J_AB=J[0], J_BC=J[1], J_CA=J[2], eps=eps, phi_deg=phi)
            summary(dict(base, quantity="R=zz_AB/t_AB", value=R, note=f"multiplet_dim={o['multiplet_dim']}, k={o['multiplet_k']}"))
            summary(dict(base, quantity="Delta=R-0.5", value=R - 0.5))
            log(f"E r={r} Sz={sz_of(nup)}: E={o['E']:.12f} S={o['S']} R={R:.6f} Delta={R - 0.5:+.6f} m=({o['m_A']:.4f},{o['m_B']:.4f},{o['m_C']:.4f}) dim={o['multiplet_dim']}")


def task_G():
    pts = [("eps0.1_phi30", couplings(0.1, 30)), ("eps0.2_phi30", couplings(0.2, 30)), ("mirror_1.1_1.1_0.8", (1.1, 1.1, 0.8))]
    for name, J in pts:
        E = {}
        for nup in (14, 15, 16):
            lv = lowest_per_k(cache, J, nup, ks=GK, vectors=True, log=log)
            write_rows("G", J, nup, lv)
            E[sz_of(nup)] = emin(lv)
        for nup in (17, 18, 19):                  # Task B levels (all 9 momenta, cached)
            E[sz_of(nup)] = emin(lowest_per_k(cache, J, nup, ks=ALLK, vectors=True, log=log))
        lo = {m: (E[4.5] - E[m]) / (4.5 - m) for m in E if m < 4.5}
        m_lo = max(lo, key=lo.get)
        hlo_B = lo[3.5]
        o = sorted(J); Jmin = o[0]
        hm = 3 * S_SPIN * (Jmin - math.sqrt(max((o[1] - Jmin) * (o[2] - Jmin), 0.0)))
        eps, phi = eps_phi(J)
        base = dict(task="G", point=name, J_AB=J[0], J_BC=J[1], J_CA=J[2], eps=eps, phi_deg=phi)
        summary(dict(base, quantity="h_lo_full_Maxwell", value=lo[m_lo], note=f"set by S^z={m_lo}; candidates {json.dumps({str(k): round(v, 6) for k, v in lo.items()})}"))
        summary(dict(base, quantity="h_lo_sector", value=m_lo))
        summary(dict(base, quantity="h_lo-h_minus_corrected", value=lo[m_lo] - hm, note=f"Task B (S^z=7/2 only): {hlo_B - hm:.6f}; change {lo[m_lo] - hlo_B:+.3e}"))
        log(f"G {name}: full-Maxwell h_lo={lo[m_lo]:.6f} (S^z={m_lo}); Task B h_lo={hlo_B:.6f}; h_lo-h-={lo[m_lo] - hm:+.6f}")


def plateau_scan(task, r, full_k_low):
    J = (1.0, r, r)
    E, where = {}, {}
    for nup in range(14, 21):
        ks = ALLK if (nup >= 16 or full_k_low) else GK
        lv = lowest_per_k(cache, J, nup, ks=ks, vectors=True, log=log)
        write_rows(task, J, nup, lv)
        E[sz_of(nup)] = emin(lv)
        where[sz_of(nup)] = min(lv.items(), key=lambda kv: kv[1]["E"])[0]
    lo = {m: (E[4.5] - E[m]) / (4.5 - m) for m in E if m < 4.5}
    hi = {m: (E[m] - E[4.5]) / (m - 4.5) for m in E if m > 4.5}
    mlo, mhi = max(lo, key=lo.get), min(hi, key=hi.get)
    Sgs, E0 = gs_spin_from_levels(E)
    eps, phi = eps_phi(J)
    base = dict(task=task, point=f"r={r}", J_AB=J[0], J_BC=J[1], J_CA=J[2], eps=eps, phi_deg=phi)
    summary(dict(base, quantity="h_lo", value=lo[mlo], note=f"set by S^z={mlo}"))
    summary(dict(base, quantity="h_hi", value=hi[mhi], note=f"set by S^z={mhi}"))
    summary(dict(base, quantity="W_27", value=hi[mhi] - lo[mlo]))
    summary(dict(base, quantity="ground_spin_S", value=Sgs, note=f"E0={E0:.12f}; lowest-level momenta {json.dumps({str(k): list(v) for k, v in where.items()})}"))
    log(f"{task} r={r}: h_lo={lo[mlo]:.6f} (S^z={mlo}) h_hi={hi[mhi]:.6f} (S^z={mhi}) W={hi[mhi] - lo[mlo]:.6f} S_gs={Sgs}")


def task_F(rs=(0.2, 0.3, 0.5)):
    for r in rs:
        plateau_scan("F", r, full_k_low=(r == 0.2))


def task_F2():
    task_F(rs=(0.7, 0.9))


TASKS.update(E=task_E, G=task_G, F=task_F, F2=task_F2)


# ======================================================================= third request (H)
TOLDEG = 1e-8


def levels_up(J, ks_half=GK, cap=5.5, need_full_half=False):
    """Lowest level per S^z: S^z = 1/2 over ks_half (all momenta if need_full_half), then S^z = 3/2, 5/2, ...
    over all momenta until the level rises above the running minimum (or S^z = cap).
    Returns (S, E0, {Sz: E}, {Sz: (k_min, lv)})."""
    E, where = {}, {}
    lv = lowest_per_k(cache, J, 14, ks=(ALLK if need_full_half else ks_half), vectors=False, log=log)
    E[0.5] = emin(lv); where[0.5] = (min(lv.items(), key=lambda kv: kv[1]["E"])[0], lv)
    m = 1.5
    while m <= cap + 1e-9:
        lv = lowest_per_k(cache, J, nup_of(m), ks=ALLK, vectors=False, log=log)
        E[m] = emin(lv); where[m] = (min(lv.items(), key=lambda kv: kv[1]["E"])[0], lv)
        if E[m] > min(E.values()) + TOLDEG:
            break
        m += 1.0
    S, E0 = gs_spin_from_levels(E, TOLDEG)
    return S, E0, E, where


def confirm_S2(J, S, where):
    """Rebuild the lowest state of the S^z = S sector (its lowest momentum) and return S from S^2."""
    k = where[S][0]
    lv = lowest_per_k(cache, J, nup_of(S), ks=[k], vectors=True, log=log)
    r = lv[k]
    write_rows("H", J, nup_of(S), {k: r})
    return r["S"], r["S2_residual"], r["H_residual"]


def H1_point(j, res=0.0004, hi0=0.02):
    Jd = lambda d: (1.0, j + d / 2, j - d / 2)
    jump = lambda d: ground_spin_jump_all(Jd(d))
    base = dict(task="H1", point=f"j={j}")
    if jump(0.0)[0]:
        summary(dict(base, quantity="d_c", value="", note="S>1/2 already at d=0: not decidable (skipped)"))
        log(f"H1 j={j}: S>1/2 at d=0, skipped")
        return None
    hi = hi0
    while not jump(hi)[0]:
        hi *= 2
        if hi > 0.0801:
            summary(dict(base, quantity="d_c", value="", note="no jump up to d=0.08"))
            return None
    lo = 0.0
    while hi - lo > res + 1e-12:
        mid = round((lo + hi) / 2, 10)
        jm, e1, e3 = jump(mid)
        log(f"H1 j={j} d={mid:.6f}: E1/2={e1:.10f} E3/2={e3:.10f} jump={jm}")
        if jm:
            hi = mid
        else:
            lo = mid
    out = {}
    for d in (lo, hi):           # final bracket: S^z = 1/2, 3/2, 5/2, all momenta, vectors (S from S^2)
        J = Jd(d)
        E = {}
        for nup in (14, 15, 16):
            lv = lowest_per_k(cache, J, nup, ks=ALLK, vectors=True, log=log)
            write_rows("H", J, nup, lv)
            E[sz_of(nup)] = emin(lv)
        Sg, E0 = gs_spin_from_levels(E, TOLDEG)
        k0 = min(lowest_per_k(cache, J, nup_of(Sg), ks=ALLK, vectors=True, log=log).items(), key=lambda kv: kv[1]["E"])
        out[d] = (Sg, k0[1]["S"], E0)
    lin = (0.056 * j, 0.058 * j); quad = (0.28 * j * j, 0.39 * j * j)
    ov = lambda a, b: not (b[1] < a[0] or b[0] > a[1])
    br = (lo, hi)
    summary(dict(base, J_AB=1, quantity="d_c", value=hi,
                 note=f"bracket ({lo:.6f},{hi:.6f}]; S(lo)={out[lo][0]} (S^2: {out[lo][1]}), S(hi)={out[hi][0]} (S^2: {out[hi][1]}); "
                      f"d_c/j={lo / j:.4f}-{hi / j:.4f}; d_c/j^2={lo / j / j:.3f}-{hi / j / j:.3f}; "
                      f"linear 0.056-0.058j=[{lin[0]:.5f},{lin[1]:.5f}] overlap={ov(br, lin)}; "
                      f"quadratic 0.28-0.39j^2=[{quad[0]:.5f},{quad[1]:.5f}] overlap={ov(br, quad)}"))
    log(f"H1 j={j}: d_c in ({lo:.6f},{hi:.6f}]  S {out[lo][0]}->{out[hi][0]}  linear overlap={ov(br, lin)} quadratic overlap={ov(br, quad)}")
    return lo, hi


def ground_spin_jump_all(J):
    """Task H1 rule: S^z = 1/2 and 3/2 both over all 9 momenta; S > 1/2 iff E(3/2) <= E(1/2) + 1e-8."""
    e1 = emin(lowest_per_k(cache, J, 14, ks=ALLK, vectors=False, log=log))
    e3 = emin(lowest_per_k(cache, J, 15, ks=ALLK, vectors=False, log=log))
    return e3 <= e1 + TOLDEG, e1, e3


def H2_point(j, d_first, step=0.0025, res=0.0005):
    """Track the ground spin from the first jump upward; bisect every change of S to res."""
    Jd = lambda d: (1.0, j + d / 2, j - d / 2)
    memo = {}

    def S_at(d):
        d = round(d, 10)
        if d not in memo:
            S, E0, E, where = levels_up(Jd(d))
            memo[d] = (S, E0, E, where)
            log(f"H2 j={j} d={d:.6f}: S={S} E0={E0:.10f} levels={ {k: round(v, 8) for k, v in E.items()} }")
        return memo[d][0]

    dmax = 2 * j * j
    grid = [round(d_first + i * step, 10) for i in range(0, 200) if d_first + i * step <= dmax + 1e-12]
    thresholds = {}
    prev_d, prev_S = None, None
    for d in grid:
        S = S_at(d)
        if prev_S is not None and S > prev_S:
            for target in [prev_S + 1 + i for i in range(int(round(S - prev_S)))]:
                lo, hi = prev_d, d
                while hi - lo > res + 1e-12:
                    mid = round((lo + hi) / 2, 10)
                    if S_at(mid) >= target:
                        hi = mid
                    else:
                        lo = mid
                thresholds[target] = (lo, hi)
        elif prev_S is not None and S < prev_S:
            log(f"H2 j={j}: S decreased {prev_S}->{S} at d={d} (non-monotonic)")
            thresholds.setdefault("non_monotonic", []).append((prev_d, d, prev_S, S))
        prev_d, prev_S = d, S
        if S >= 4.5:
            break
    # S^2 confirmation at both ends of every bracket, rows for the energy levels at every evaluated d
    for d, (S, E0, E, where) in sorted(memo.items()):
        for m, (k, lv) in where.items():
            write_rows("H", Jd(d), nup_of(m), lv)
    base = dict(task="H2", point=f"j={j}", J_AB=1)
    for target, br in sorted((t, b) for t, b in thresholds.items() if t != "non_monotonic"):
        s_lo = confirm_S2(Jd(br[0]), memo[round(br[0], 10)][0], memo[round(br[0], 10)][3])
        s_hi = confirm_S2(Jd(br[1]), memo[round(br[1], 10)][0], memo[round(br[1], 10)][3])
        summary(dict(base, quantity=f"d(->S={target})", value=br[1],
                     note=f"bracket ({br[0]:.6f},{br[1]:.6f}]; S^2 check lo={s_lo[0]}, hi={s_hi[0]}; d/j^2={br[1] / j / j:.3f}"))
        log(f"H2 j={j}: d(->S={target}) in ({br[0]:.6f},{br[1]:.6f}]  S^2 check {s_lo[0]}/{s_hi[0]}")
    reached = [t for t in thresholds if t != "non_monotonic" and t >= 4.5]
    if reached:
        br = thresholds[4.5]
        summary(dict(base, quantity="d_full", value=br[1],
                     note=f"S=9/2 (N/6) reached in ({br[0]:.6f},{br[1]:.6f}]; d_full/j^2={br[0] / j / j:.3f}-{br[1] / j / j:.3f}; "
                          f"theory 0.34-0.37j^2=[{0.34 * j * j:.5f},{0.37 * j * j:.5f}]; classical j^2={j * j:.5f}"))
    else:
        last = max(memo)
        summary(dict(base, quantity="d_full", value="", note=f"S=9/2 not reached up to d={last:.6f} (2j^2={dmax:.5f}); S there={memo[last][0]}"))
    if "non_monotonic" in thresholds:
        summary(dict(base, quantity="non_monotonic", value=json.dumps(thresholds["non_monotonic"])))
    return thresholds


def task_H():
    dc = {}
    for j in (0.10, 0.125):
        r = H1_point(j)
        if r:
            dc[j] = r[1]
    # first-jump upper bracket ends from task C for j = 0.15, 0.2
    dc.setdefault(0.15, 0.00875); dc.setdefault(0.2, 0.0125)
    for j in (0.10, 0.15, 0.20):
        if j in dc:
            H2_point(j, dc[j])
    H1_point(0.175)


TASKS.update(H=task_H, H1=lambda: [H1_point(j) for j in (0.10, 0.125)], H175=lambda: H1_point(0.175))

if __name__ == "__main__":
    for t in sys.argv[1:] or ["A", "B", "D", "C"]:
        t0 = time.time()
        log(f"=== task {t} start")
        TASKS[t]()
        log(f"=== task {t} done in {time.time() - t0:.0f}s")
