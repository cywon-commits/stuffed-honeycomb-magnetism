"""Ground-state observables from ed_sym pair-orbit output.

Reconstructs G_ij = <S_i.S_j> (G_ii = 3/4) from the orbit sums printed by
`ed_sym ... --pairs pairs.txt` (valid because the state lies in a 1D sector,
so G is constant on every orbit of the symmetry group), then evaluates the
same quantities, with the same definitions, as
src/workspaces/finite_size_18_27_36pilot/observables.py::values().
"""
import json, math, pathlib, re, sys
import numpy as np
from symmetry import Cluster


def parse_log(text):
    d = {"orbits": {}}
    for l in text.splitlines():
        if l.startswith("ORBIT "):
            _, k, v = l.split()
            d["orbits"][int(k)] = float(v)
        elif l.startswith("GROUND"):
            p = l.split()
            d["E0"] = float(p[2]); d["residual"] = float(p[4])
        elif l.startswith("LOADED"):
            d["E0"] = float(l.split()[2])
        elif l.startswith("S2 "):
            d["S2_solver"] = float(l.split()[1])
        elif l.startswith("PASS1"):
            p = l.split(); d["lanczos_iterations"] = int(p[2]); d["ritz_resid"] = float(p[8])
        elif l.startswith("BASIS"):
            d["dim"] = int(l.split()[2])
    return d


def correlation_matrix(c, pairs_path, orbit_sums):
    lines = pathlib.Path(pairs_path).read_text().splitlines()[1:]
    pairs = [tuple(map(int, l.split())) for l in lines]
    size = {}
    for _, _, o, i, j in pairs:
        size[o] = size.get(o, 0) + 1
    G = np.eye(c.N) * 0.75
    for _, _, o, i, j in pairs:
        G[i, j] = G[j, i] = orbit_sums[o] / size[o]
    return G


def values(c, G):
    N, nh = c.N, c.Ncell
    nr = 2 * nh
    h = np.arange(2, N, 3)
    rim = np.array([i for i in range(N) if i % 3 != 2])
    sgn = np.array([1 if i % 3 == 0 else (-1 if i % 3 == 1 else 0) for i in range(N)])
    H2 = float(G[np.ix_(h, h)].sum()); R2 = float(G[np.ix_(rim, rim)].sum()); HR = float(G[np.ix_(h, rim)].sum())
    S2 = float(G.sum())
    S = 0.5 * (math.sqrt(1 + 4 * S2) - 1)
    Sround = round(2 * S) / 2
    neel = float(sgn @ G @ sgn / nr ** 2)
    cells = np.array(c.cells, float)
    hs = []
    for q in c.momenta:
        phase = np.exp(-2j * np.pi * (cells @ np.array(q))) / np.sqrt(nh)
        hs.append(float(np.vdot(phase, G[np.ix_(h, h)] @ phase).real))
    peak = max(hs)
    peakq = [list(q) for q, s in zip(c.momenta, hs) if abs(s - peak) < 1e-8]
    return {"S_from_S2": S, "S": Sround, "S2": S2, "m_tot": Sround / N, "m_tot_rms": math.sqrt(Sround * (Sround + 1)) / N,
            "H2": H2, "R2": R2, "HR": HR, "spin_sum_rule_error": abs(H2 + R2 + 2 * HR - Sround * (Sround + 1)),
            "p_H": (H2 - 0.75 * nh) / (nh * (nh - 1) / 4), "m_R_uniform2": R2 / nr ** 2, "m_N2": neel,
            "legacy_Neel_FGamma": nr ** 2 * neel / nh, "hub_S_Gamma": hs[0], "hub_peak": peak, "hub_peak_q": peakq,
            "hub_structure_factor": {json.dumps([round(x, 6) for x in q]): s for q, s in zip(c.momenta, hs)},
            "energy_per_site_note": "E0/N given separately"}


def analyze(cluster, pairs_path, log_text):
    c = Cluster(cluster)
    d = parse_log(log_text)
    G = correlation_matrix(c, pairs_path, d["orbits"])
    v = values(c, G)
    # energy from G as an independent check of the correlations
    r = None
    m = re.search(r"^SPEC .* r=([0-9.eE+-]+)", log_text, re.M)
    if m:
        r = float(m.group(1))
        Eg = sum((1.0 if t == 0 else r) * G[a, b] for a, b, t in c.bonds)
        v["energy_from_G"] = Eg
    v.update({k: d[k] for k in ("E0", "residual", "lanczos_iterations", "ritz_resid", "dim", "S2_solver") if k in d})
    if "E0" in d:
        v["E0_per_site"] = d["E0"] / c.N
        if "energy_from_G" in v:
            v["energy_from_G_error"] = abs(v["energy_from_G"] - d["E0"])
    v["r"] = r
    return v


if __name__ == "__main__":
    cluster, pairs, log = sys.argv[1:4]
    print(json.dumps(analyze(cluster, pairs, pathlib.Path(log).read_text()), indent=1))
