"""Figures 2-5 of the manuscript from the data in data/field_dice (run from manuscript/figures)."""
import csv, re, sys, pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

ROOT = pathlib.Path(__file__).resolve().parents[2]
D = ROOT / "data" / "field_dice"
sys.path.insert(0, str(ROOT / "src" / "field_dice"))
from simplex import couplings, eps_max            # noqa: E402
from rotor_aniso import G                          # noqa: E402

# validated categorical order (light mode); grays for theory/reference lines
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, MUTED, GRID = "#222222", "#6b6b6b", "#d9d9d9"
plt.rcParams.update({"font.size": 8, "axes.labelsize": 8, "legend.fontsize": 7, "xtick.labelsize": 7,
                     "ytick.labelsize": 7, "axes.edgecolor": MUTED, "axes.linewidth": 0.6,
                     "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK,
                     "lines.linewidth": 1.4, "legend.frameon": False, "savefig.dpi": 300})
S = 0.5
rows = lambda f: list(csv.DictReader(open(f)))


def label(ax, s):
    ax.text(0.02, 0.97, s, transform=ax.transAxes, va="top", ha="left", fontsize=9, color=INK, weight="bold")


def save(fig, name):
    for ext in ("pdf", "png"):
        fig.savefig(pathlib.Path(__file__).with_name(f"{name}.{ext}"), bbox_inches="tight")


# ---------------------------------------------------------------- Fig. 2
def eps_collinear(phi):
    """Smallest eps where 1/J_CA = 1/J_AB + 1/J_BC (classical non-collinear/collinear boundary)."""
    f = lambda e: (lambda J: 1 / J[2] - 1 / J[0] - 1 / J[1])(couplings(e, phi))
    lo, hi = 1e-9, eps_max(phi) * (1 - 1e-9)
    if f(hi) < 0:
        return np.nan
    for _ in range(80):
        m = (lo + hi) / 2
        lo, hi = (m, hi) if f(m) < 0 else (lo, m)
    return hi


def fig2():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.9), gridspec_kw={"width_ratios": [1.15, 1]})
    ph = np.linspace(0, np.pi / 3, 200)
    em = eps_max(ph)
    a.plot(em * np.cos(ph), em * np.sin(ph), color=MUTED, lw=0.8)
    a.plot([0, em[0]], [0, 0], color=MUTED, lw=0.8); a.plot([0, em[-1] * np.cos(np.pi / 3)], [0, em[-1] * np.sin(np.pi / 3)], color=MUTED, lw=0.8)
    ec = np.array([eps_collinear(p) for p in ph])
    a.plot(ec * np.cos(ph), ec * np.sin(ph), color=INK, lw=1.4, label=r"classical boundary, Eq. (7)")
    R = [x for x in rows(D / "simplex_scan_18.csv") if x["shape"] == "K18"]
    cmap = plt.get_cmap("Blues")
    for x in R:
        e, p = float(x["eps"]), np.radians(float(x["phi_deg"]))
        s = int(x["S"]) / 3
        a.scatter(e * np.cos(p), e * np.sin(p), s=26, c=[cmap(0.15 + 0.85 * s)], edgecolors="white", linewidths=0.5, zorder=3)
    sm = plt.cm.ScalarMappable(cmap=matplotlib.colors.LinearSegmentedColormap.from_list("b", [cmap(0.15), cmap(1.0)]),
                               norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, ax=a, fraction=0.04, pad=0.02); cb.set_label(r"$S_{\rm tot}/(N/6)$, $N=18$", color=INK)
    cb.outline.set_linewidth(0.4)
    a.annotate("triangle", (0, 0), (0.05, -0.07), fontsize=7, color=MUTED)
    a.annotate("honeycomb", (em[0], 0), (em[0] - 0.33, -0.07), fontsize=7, color=MUTED)
    a.annotate("dice", (em[-1] / 2, em[-1] * np.sin(np.pi / 3)), (0.16, 0.45), fontsize=7, color=MUTED)
    a.text(0.66, 0.235, "Lieb–Mattis edge", fontsize=7, color=MUTED, rotation=-26, ha="center")
    a.set_aspect("equal"); a.set_xlim(-0.08, 1.08); a.set_ylim(-0.1, 0.55)
    a.set_xlabel(r"$\varepsilon\cos\varphi$"); a.set_ylabel(r"$\varepsilon\sin\varphi$")
    a.legend(loc="upper right", handlelength=1.5)
    label(a, "(a)")

    W = rows(D / "pd_wedge_18.csv")
    for shp, col, mk, lab in (("K18", C2, "s", r"$N=18$ K18 ($0\to3=N/6$)"), ("C18", C4, "^", r"$N=18$ C18 ($0\to1$)")):
        pts = [(float(x["j"]), float(x["delta_c"])) for x in W if x["shape"] == shp and x["delta_c"]]
        j, d = np.array(pts).T
        b.plot(j, d, mk, color=col, ms=5, mec="white", mew=0.5, label=lab)
    first, full = [], []
    for x in rows(D / "n27" / "summary_n27.csv"):
        if x["task"] in ("C", "H1") and x["quantity"] == "d_c" and x["value"]:
            first.append(x)
        if x["task"] == "H2" and x["quantity"] == "d_full":
            full.append(x)
    for lst, mfc, lab in ((first, "white", r"$N=27$ first jump ($1/2\to3/2$)"), (full, C3, r"$N=27$ ferrimagnet ($3/2\to9/2$)")):
        for n, x in enumerate(lst):
            jj = float(x["point"].split("=")[1]); lo, hi = map(float, re.search(r"[\[(]([\d.]+),\s*([\d.]+)\]", x["note"]).groups())
            b.errorbar(jj, (lo + hi) / 2, yerr=(hi - lo) / 2, fmt="D", color=C3, mfc=mfc, ms=5, mew=1, capsize=2, lw=1,
                       label=lab if n == 0 else None)
    jg = np.linspace(0, 0.32, 100)
    b.plot(jg, 0.164 * jg, ls="--", color=C3, lw=0.9)
    b.plot(jg, 0.209 * jg, ls="--", color=C2, lw=0.9)
    b.fill_between(jg, 0.34 * jg**2, 0.37 * jg**2, color=INK, alpha=0.25, lw=0, label=r"hub model $0.34$–$0.37\,J'^2/J$")
    b.plot(jg, jg**2, ls=":", color=INK, lw=1, label=r"classical $J'^2/J$")
    b.set_xlim(0, 0.32); b.set_ylim(0, 0.072)
    b.set_xlabel(r"$J'=(J_{BC}+J_{CA})/2$  ($J_{AB}=J=1$)"); b.set_ylabel(r"threshold $|J_{BC}-J_{CA}|_c$")
    b.grid(color=GRID, lw=0.4); b.legend(loc="upper left", bbox_to_anchor=(0.06, 1.0), fontsize=6.3)
    label(b, "(b)")
    fig.tight_layout(); save(fig, "fig2_zero_field")


# ---------------------------------------------------------------- Fig. 4 (plateau edges)
def hpm(J):
    Js = sorted(J); jm, j1, j2 = Js[0], Js[1], Js[2]
    r = np.sqrt((j1 - jm) * (j2 - jm)); return 3 * S * (jm - r), 3 * S * (jm + r)


def fig3():
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 2.7))
    R = [x for x in rows(D / "simplex_scan_18.csv") if int(x["phi_deg"]) >= 10 and float(x["eps"]) > 0]
    B = {}
    for x in rows(D / "n27" / "summary_n27.csv"):
        if x["task"] == "B": B.setdefault(x["point"], {"J": (float(x["J_AB"]), float(x["J_BC"]), float(x["J_CA"]))})[x["quantity"]] = float(x["value"])
    for shp, col, mk in (("K18", C1, "o"), ("C18", C2, "s")):
        hp = [hpm((float(x["J_AB"]), float(x["J_BC"]), float(x["J_CA"]))) for x in R if x["shape"] == shp]
        hi = [float(x["h_hi"]) for x in R if x["shape"] == shp]; lo = [float(x["h_lo"]) for x in R if x["shape"] == shp]
        hp = np.array(hp)
        a.plot(hp[:, 1], hi, mk, color=col, ms=4, mec="white", mew=0.4, label=f"$N=18$ ({shp})")
        m = hp[:, 0] > 0
        b.plot(hp[m, 0], np.array(lo)[m], mk, color=col, ms=4, mec="white", mew=0.4, label=f"$N=18$ ({shp})")
    hp27 = np.array([[v["h_minus_classical"], v["h_plus_classical"], v["h_lo"], v["h_hi"]] for v in B.values()])
    a.plot(hp27[:, 1], hp27[:, 3], "D", color=C3, ms=5, mec="white", mew=0.5, label="$N=27$")
    m = hp27[:, 0] > 0
    b.plot(hp27[m, 0], hp27[m, 2], "D", color=C3, ms=5, mec="white", mew=0.5, label="$N=27$")
    x = np.linspace(1.2, 2.25, 10)
    a.plot(x, x, ":", color=INK, lw=1, label=r"classical $h=h_+$")
    a.plot(x, 1.515 * x, "--", color=MUTED, lw=0.9); a.plot(x, 1.44 * x, "--", color=MUTED, lw=0.9)
    a.text(1.21, 1.515 * 1.21 + 0.08, r"$1.51\,h_+$", fontsize=7, color=MUTED)
    a.text(1.30, 1.44 * 1.30 - 0.16, r"$1.44\,h_+$", fontsize=7, color=MUTED)
    a.set_xlabel(r"classical upper edge $h_+$"); a.set_ylabel(r"ED upper edge $h_{\rm hi}$")
    x = np.linspace(0, 1.5, 10)
    b.plot(x, x, ":", color=INK, lw=1, label=r"classical $h=h_-$"); b.plot(x, x - 0.11, "--", color=MUTED, lw=0.9, label=r"$h_- - 0.11$")
    b.set_xlabel(r"classical lower edge $h_-$"); b.set_ylabel(r"ED lower edge $h_{\rm lo}$")
    for ax, l in ((a, "(a)"), (b, "(b)")):
        ax.grid(color=GRID, lw=0.4); ax.legend(loc="lower right"); label(ax, l)
    fig.tight_layout(); save(fig, "fig4_plateau_edges")


# ---------------------------------------------------------------- Fig. 4
def fig4():
    fig, (a, b, c) = plt.subplots(1, 3, figsize=(7.0, 2.5), gridspec_kw={"width_ratios": [1.35, 1, 1]})
    L = rows(D / "obd_flop_lswt.csv"); rl = np.array([float(x["r"]) for x in L]); Kl = np.array([float(x["K_S0.5"]) for x in L])
    rg = np.linspace(0.1, 0.9, 81); Kg = np.interp(rg, rl, Kl)
    chi = 1 / 9
    V = rows(D / "rotor_aniso_vs_ed.csv")
    r = np.array([float(x["r"]) for x in V])
    for N, col, mk, key in ((18, C2, "s", None), (27, C3, "D", "W27_ED")):
        I = chi * N
        a.plot(rg, [G(k * N * I)[0] / I for k in Kg], color=col, lw=1.2)
        if N == 18:
            a.plot(r, [float(x["W18_K18"]) for x in V], "s", color=col, ms=4, mec="white", mew=0.4)
            a.plot(r, [float(x["W18_C18"]) for x in V], "^", color=col, ms=4, mec="white", mew=0.4)
        else:
            a.plot(r, [float(x[key]) for x in V], mk, color=col, ms=5, mec="white", mew=0.5)
    a.plot(rg, 2 * np.sqrt(2 * Kg / chi), color=INK, lw=1.4)
    a.text(0.62, 1.22, r"bulk $2\sqrt{2K/\chi}$", fontsize=7, color=INK)
    a.text(0.12, 1.09, "$N=18$", fontsize=7, color=INK); a.text(0.12, 0.62, "$N=27$", fontsize=7, color=INK)
    a.set_xlim(0.1, 0.92); a.set_ylim(0, 1.4)
    a.set_xlabel(r"$r=J'/J$ (honeycomb mirror)"); a.set_ylabel(r"$1/3$-plateau width $W/J$")
    a.legend(handles=[Line2D([], [], color=MUTED, lw=1.2, label="rotor theory (no fit)"),
                      Line2D([], [], ls="", marker="s", color=C2, label="ED $N=18$ (K18; C18 $\\triangle$)"),
                      Line2D([], [], ls="", marker="D", color=C3, label="ED $N=27$")], loc="lower right")
    E27 = {}
    for x in rows(D / "n27" / "task_E.csv"):
        E27.setdefault((float(x["J_BC"]), x["Sz"]), []).append(x)
    rows18 = rows(D / "rim_anisotropy_18.csv")
    rr = [0.4, 0.6, 0.8]
    d18 = {s: [np.mean([float(x["zz_AB"]) / float(x["xxyy_AB"]) - 0.5 for x in rows18 if x["shape"] == s and float(x["r"]) == q and x["Sz"] == "3"]) for q in rr] for s in ("K18", "C18")}
    summ = {float(x["J_BC"]): float(x["value"]) for x in rows(D / "n27" / "summary_n27.csv") if x["task"] == "E" and x["quantity"].startswith("Delta") and x["point"].endswith("Sz=4.5")}
    if not summ:
        summ = {x: v for x, v in zip((0.2, 0.4, 0.6, 0.8, 1.0), (0.0078, 0.038, 0.1208, 0.3289, 0.1093))}
    r27 = sorted(k for k in summ if k < 1); d27 = [summ[k] for k in r27]
    b.plot(rr, d18["K18"], "s-", color=C2, ms=4, mec="white", mew=0.4, lw=1, label="$N=18$ (K18)")
    b.plot(rr, d18["C18"], "^-", color=C2, ms=4, mec="white", mew=0.4, lw=1, alpha=0.8, label="$N=18$ (C18)")
    b.plot(r27, d27, "D-", color=C3, ms=5, mec="white", mew=0.5, lw=1, label="$N=27$")
    b.set_xlabel(r"$r$"); b.set_ylabel(r"rim anisotropy $\Delta=\langle zz\rangle/\langle xx{+}yy\rangle-\frac{1}{2}$")
    b.set_xlim(0.15, 0.85); b.set_ylim(0, 0.36); b.legend(loc="upper left", bbox_to_anchor=(0.04, 0.9))
    F = [x for x in rows(D / "rotor_fit_ladder.csv") if x["cluster"] == "27" and x["K_eff"] not in ("", "nan") and float(x["r"]) >= 0.55]
    rf = [float(x["r"]) for x in F]
    c.plot(rf, [float(x["K_eff"]) / float(x["K_LSWT"]) for x in F], "D-", color=C3, ms=5, mec="white", mew=0.5, lw=1,
           label=r"$K_{\rm eff}/K_{\rm LSWT}$")
    c.plot(rf, [float(x["I_over_N9"]) for x in F], "o--", color=C1, ms=4, mec="white", mew=0.4, lw=1, label=r"$I/(N/9)$")
    c.axhline(1, color=INK, lw=0.9, ls=":")
    c.set_xlim(0.5, 0.95); c.set_ylim(0, 1.4); c.set_xlabel(r"$r$"); c.set_ylabel(r"rotor fit, $N=27$")
    c.legend(loc="lower right")
    for ax, l in ((a, "(a)"), (b, "(b)"), (c, "(c)")):
        ax.grid(color=GRID, lw=0.4); label(ax, l)
    fig.tight_layout(); save(fig, "fig5_mirror_plateau")


# ---------------------------------------------------------------- Fig. 5
def fig5():
    from lswt3 import zero_point
    from simplex import f_classical
    Ehc, Ehq = zero_point((1, -1, 1), (1.0, 0.0, 0.0))
    fig, ax = plt.subplots(figsize=(3.4, 2.6))
    R = [x for x in rows(D / "simplex_scan_18.csv") if x["shape"] == "K18"]
    for phi, mk in (("30", "s"), ("60", "o")):
        P = sorted([x for x in R if x["phi_deg"] == phi and x["S"] == "3"], key=lambda x: float(x["eps"]))
        e = [float(x["eps"]) for x in P]
        fl, fe = [], []
        for x in P:
            J = (float(x["J_AB"]), float(x["J_BC"]), float(x["J_CA"])); k = int(np.argmin(J))
            sig = [(1, 1, -1), (-1, 1, 1), (1, -1, 1)][k]
            Ec, Eq = zero_point(sig, J, n=24)
            fl.append((Eq / (sum(J) * Ehq)) / (Ec / (sum(J) * Ehc))); fe.append(float(x["f"]) / float(x["f_cl"]))
        ax.plot(e, fl, "-", marker=mk, color=C1, ms=4, mec="white", mew=0.4, lw=1.1)
        ax.plot(e, fe, "-", marker=mk, color=C2, ms=4, mec="white", mew=0.4, lw=1.1)
    for x in rows(D / "n27" / "summary_n27.csv"):
        if x["task"] == "D" and x["quantity"] == "f":
            mk = "o" if x["point"] == "dice" else "s"
            fc = f_classical(float(x["J_AB"]), float(x["J_BC"]), float(x["J_CA"]))
            ax.plot(float(x["eps"]), float(x["value"]) / fc, mk, color=C3, ms=6, mec="white", mew=0.5)
    ax.set_xlabel(r"$\varepsilon$ (collinear region)"); ax.set_ylabel(r"$f/f_{\rm cl}$")
    ax.set_ylim(0.8, 1.06); ax.grid(color=GRID, lw=0.4)
    ax.legend(handles=[Line2D([], [], color=C1, lw=1.1, label="LSWT"), Line2D([], [], color=C2, lw=1.1, label="ED $N=18$"),
                       Line2D([], [], ls="", marker="o", color=C3, label="ED $N=27$"),
                       Line2D([], [], ls="", marker="o", color=MUTED, label=r"$\varphi=60^\circ$ (dice mirror)"),
                       Line2D([], [], ls="", marker="s", color=MUTED, label=r"$\varphi=30^\circ$")], loc="upper right", ncol=2,
              columnspacing=0.8)
    fig.tight_layout(); save(fig, "fig3_resonance")


if __name__ == "__main__":
    for f in sys.argv[1:] or ["2", "3", "4", "5"]:
        globals()["fig" + f]()
