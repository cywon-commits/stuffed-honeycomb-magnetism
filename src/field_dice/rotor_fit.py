"""Fit the rim-rotor parameters I and D separately from the S^z ladder around the plateau sector.
Symmetric gaps Delta_m = (E_{N/6+m} + E_{N/6-m})/2 - E_{N/6} (m = 1, 2) cancel the field offset of the
plateau centre; for H = L^2/2I - D n_z^2 the ratio Delta_2/Delta_1 depends on x = D I only (3 at x = 0)."""
import csv, sys, pathlib
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(pathlib.Path(__file__).parent))
from rotor_aniso import levels

D = pathlib.Path(__file__).resolve().parents[2] / "data" / "field_dice"


def gaps(x):
    e0, e1, e2 = (levels(x, m)[0] for m in (0, 1, 2))
    return e1 - e0, e2 - e0


def fit(d1, d2):
    R = d2 / d1
    g = lambda x: gaps(x)[1] / gaps(x)[0] - R
    if g(0.0) * g(80.0) > 0:
        return None
    x = brentq(g, 0.0, 80.0)
    I = gaps(x)[0] / d1
    return x, I, x / I


def ladders():
    out = {}
    for f in ("n27/task_A.csv", "n27/task_F.csv"):
        for x in csv.DictReader(open(D / f)):
            if float(x["J_AB"]) != 1.0 or x["J_BC"] != x["J_CA"]:
                continue
            k = (27, round(float(x["J_BC"]), 3)); out.setdefault(k, {})
            sz = float(x["Sz"]); out[k][sz] = min(out[k].get(sz, 1e9), float(x["E"]))
    for x in csv.DictReader(open(D / "sz_scan_18.csv")):
        k = (x["shape"], round(float(x["Jp"]), 3))
        out[k] = {float(m): float(x[f"E{m}"]) for m in range(10)}
    return out


if __name__ == "__main__":
    K = {round(float(x["r"]), 2): float(x["K_S0.5"]) for x in csv.DictReader(open(D / "obd_flop_lswt.csv"))}
    rows = []
    for (c, r), e in sorted(ladders().items(), key=lambda t: (str(t[0][0]), t[0][1])):
        N = 27 if c == 27 else 18; p = N / 6
        if not all(p + m in e for m in (-2, -1, 0, 1, 2)) or r < 0.15 or r > 0.95:
            continue
        d1 = (e[p + 1] + e[p - 1]) / 2 - e[p]; d2 = (e[p + 2] + e[p - 2]) / 2 - e[p]
        res = fit(d1, d2)
        x, I, Dp = res if res else (np.nan,) * 3
        Kl = K.get(round(r, 2), np.nan)
        rows.append((c, r, d1, d2, d2 / d1, x, I, I / (N / 9), Dp / N, Kl))
        print(f"{str(c):>4} r={r:.2f} D1={d1:.4f} D2={d2:.4f} ratio={d2/d1:.3f} x={x:.3f} I/(N/9)={I/(N/9):.3f} K_eff={Dp/N:.4f} K_LSWT={Kl:.4f}")
    with open(D / "rotor_fit_ladder.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["cluster", "r", "Delta1", "Delta2", "ratio", "x", "I", "I_over_N9", "K_eff", "K_LSWT"]); w.writerows(rows)
