"""Fig. 1: (a) three-coupling triangular lattice, (b) coupling triangle with its S3 structure.
Run from the repository root: python3 manuscript/figures/fig1_model.py"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

plt.rcParams.update({"font.size": 8, "font.family": "serif", "mathtext.fontset": "cm"})
fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.8, 3.1), gridspec_kw={"width_ratios": [1, 1.15]})

# (a) lattice: sublattice index (n1 - n2) mod 3 -> A, B, C
a1, a2 = np.array([1.0, 0.0]), np.array([0.5, np.sqrt(3) / 2])
col = {0: "#1f77b4", 1: "#ff7f0e", 2: "#2ca02c"}          # A, B, C
style = {frozenset({0, 1}): ("k", "-", 1.6), frozenset({1, 2}): ("0.35", "--", 1.1), frozenset({2, 0}): ("0.55", ":", 1.1)}
pts = {}
for n1 in range(-4, 6):
    for n2 in range(-1, 5):
        p = n1 * a1 + n2 * a2
        if 0 <= p[0] <= 5.0 and -0.01 <= p[1] <= 3.5:
            pts[(n1, n2)] = p
for (n1, n2), p in pts.items():
    for d in [(1, 0), (0, 1), (-1, 1)]:
        q = (n1 + d[0], n2 + d[1])
        if q in pts:
            s, t = (n1 - n2) % 3, (q[0] - q[1]) % 3
            c, ls, lw = style[frozenset({s, t})]
            ax.plot(*zip(p, pts[q]), color=c, ls=ls, lw=lw, zorder=1)
for (n1, n2), p in pts.items():
    ax.scatter(*p, s=34, color=col[(n1 - n2) % 3], edgecolor="k", lw=0.5, zorder=2)
for lab, s, y in [("A", 0, 0), ("B", 1, 0), ("C", 2, 0)]:
    pass
ax.plot([], [], "k-", lw=1.6, label=r"$J_{AB}$"); ax.plot([], [], color="0.35", ls="--", label=r"$J_{BC}$")
ax.plot([], [], color="0.55", ls=":", label=r"$J_{CA}$")
for s, name in col.items():
    ax.scatter([], [], s=34, color=name, edgecolor="k", lw=0.5, label="ABC"[s])
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=6, frameon=False, handlelength=1.6, columnspacing=0.8)
ax.set_aspect("equal"); ax.axis("off"); ax.set_title("(a)", loc="left")

# (b) coupling triangle in the E plane: x = eps cos phi, y = eps sin phi
def J_of(x, y):
    return [1 + 2 * (x * np.cos(2 * np.pi * k / 3) + y * np.sin(2 * np.pi * k / 3)) for k in range(3)]
hexv = [np.array([np.cos(t), np.sin(t)]) for t in (0, 2 * np.pi / 3, 4 * np.pi / 3)]
bx.add_patch(Polygon(hexv, closed=True, fill=False, ec="k", lw=1.0))
X, Y = np.meshgrid(np.linspace(-0.6, 1.05, 500), np.linspace(-0.95, 0.95, 500))
J = np.array(J_of(X, Y)); inside = (J > 0).all(0)
with np.errstate(divide="ignore", invalid="ignore"):
    u = 1 / J; s = np.sort(u, axis=0)
    canted = inside & (s[2] <= s[0] + s[1])
bx.contourf(X, Y, np.where(inside, canted.astype(float), np.nan), levels=[-0.5, 0.5, 1.5],
            colors=["#dbe9f6", "#f6dcc4"], zorder=0)
for t in np.arange(6) * np.pi / 3:                          # mirror lines
    bx.plot([0, 1.0 * np.cos(t) if t % (2 * np.pi / 3) == 0 else 0.5 * np.cos(t)],
            [0, 1.0 * np.sin(t) if t % (2 * np.pi / 3) == 0 else 0.5 * np.sin(t)], color="0.5", lw=0.6, ls="-.")
fd = [np.array([0, 0]), np.array([1, 0]), 0.5 * np.array([np.cos(np.pi / 3), np.sin(np.pi / 3)])]
bx.add_patch(Polygon(fd, closed=True, fill=False, ec="#d62728", lw=1.4, zorder=3))
for v, lab in zip(hexv, ["honeycomb\n" r"($J_{AB}$ only)", "honeycomb\n" r"($J_{BC}$)", "honeycomb\n" r"($J_{CA}$)"]):
    bx.scatter(*v, s=26, marker="h", color="k", zorder=4)
for k in range(3):
    m = 0.5 * np.array([np.cos(np.pi / 3 + 2 * np.pi * k / 3), np.sin(np.pi / 3 + 2 * np.pi * k / 3)])
    bx.scatter(*m, s=22, marker="D", color="k", zorder=4)
bx.scatter(0, 0, s=30, marker="^", color="k", zorder=4)
bx.annotate("honeycomb", (1, 0), (0.78, -0.22), fontsize=7, ha="center")
bx.annotate("dice", fd[2], (0.33, 0.52), fontsize=7, ha="center")
bx.annotate("triangle", (0, 0), (-0.2, -0.13), fontsize=7, ha="center")
bx.text(0.55, 0.06, r"$\varphi=0$" "\n" r"($r<1$)", fontsize=6.5, ha="center")
bx.text(0.02, 0.27, r"$\varphi=\pi/3$" "\n" r"($r>1$)", fontsize=6.5, ha="right")
bx.text(0.72, 0.30, r"$J_{CA}=0$" "\n Lieb–Mattis", fontsize=6.5, ha="left", rotation=-30)
from matplotlib.patches import Patch
bx.legend(handles=[Patch(fc="#f6dcc4", label="non-collinear (1/J triangle inequality)"),
                   Patch(fc="#dbe9f6", label=r"collinear Lieb–Mattis, $M=M_s/3$"),
                   Patch(fc="none", ec="#d62728", label="fundamental domain")],
          loc="upper center", bbox_to_anchor=(0.5, 0.02), ncol=1, frameon=False, fontsize=6.5)
bx.set_xlim(-0.62, 1.12); bx.set_ylim(-1.2, 0.95); bx.set_aspect("equal"); bx.axis("off")
bx.set_title("(b)", loc="left")
fig.tight_layout()
for ext in ("pdf", "png"):
    fig.savefig(f"manuscript/figures/fig1_model.{ext}", dpi=300)
