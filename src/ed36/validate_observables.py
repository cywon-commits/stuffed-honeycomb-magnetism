"""Validate pair-orbit observables (S, <S_H^2>, <S_R^2>, <S_H.S_R>, Neel, hub S(q))
against the archived 18-spin ground-state scan (gs_18_scan.csv)."""
import csv, json, pathlib, subprocess, sys, tempfile
from symmetry import Cluster
from analyze import analyze

HERE = pathlib.Path(__file__).resolve().parent
REF = HERE.parent.parent / "data" / "workspaces" / "finite_size_18_27_36pilot" / "gs_18_scan.csv"
EXE = HERE / "ed_sym"
KEYS = ["energy", "S", "H2", "R2", "HR", "p_H", "m_R_uniform2", "m_N2", "hub_S_Gamma", "hub_peak"]


def run(spec, r, extra):
    return subprocess.run([str(EXE), str(spec), "--r", repr(r), "--tol", "1e-10", "--maxit", "600"] + extra,
                          capture_output=True, text=True, check=True).stdout


def main():
    ref = list(csv.DictReader(open(REF)))
    out = []
    for shape, r in [("K18", 0.3), ("C18", 0.3), ("K18", 0.5), ("C18", 0.5), ("K18", 1.0), ("C18", 1.0)]:
        row = next(x for x in ref if x["shape"] == shape and abs(float(x["r"]) - r) < 1e-12)
        c = Cluster(shape)
        with tempfile.TemporaryDirectory() as td:
            pairs = pathlib.Path(td) / "pairs.txt"
            c.write_pairs(pairs)
            best = None
            for i, s in enumerate(c.sectors(9)):
                p = pathlib.Path(td) / f"s{i}.spec"
                c.write_spec(s, 9, p)
                o = run(p, r, ["--novec"])
                E = float(next(l for l in o.splitlines() if l.startswith("PASS1")).split()[6])
                if best is None or E < best[0] - 1e-10:
                    best = (E, p, s)
            log = run(best[1], r, ["--pairs", str(pairs)])
            v = analyze(shape, pairs, log)
        mine = dict(v, energy=v["E0"])
        diffs = {k: abs(float(mine[k]) - float(row[k])) for k in KEYS}
        res = dict(cluster=shape, r=r, sector_k=list(best[2]["k"]), sector_z=best[2]["z"], diffs=diffs,
                   max_diff=max(diffs.values()), energy_from_G_error=float(v["energy_from_G_error"]),
                   S2_solver_vs_G=abs(v["S2_solver"] - v["S2"]), residual=v["residual"])
        res["pass"] = bool(res["max_diff"] < 1e-8 and res["energy_from_G_error"] < 1e-8)
        print(json.dumps(res), flush=True)
        out.append(res)
    (HERE / "results" / "validation_observables_18.json").write_text(json.dumps(out, indent=1))
    sys.exit(0 if all(o["pass"] for o in out) else 1)


if __name__ == "__main__":
    main()
