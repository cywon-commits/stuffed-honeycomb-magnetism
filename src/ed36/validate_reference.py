"""Compare ed_sym Lanczos ground energies per momentum with the archived
18-spin (K18, C18; NEW_18_SPIN_ED) and 27-spin (verified checkpoint) results.

Reference files (repository paths):
  data/workspaces/finite_size_18_27_36pilot/sector_completeness_18.csv
  data/workspaces/finite_size_18_27_36pilot/gs_18_scan.csv
  data/workspaces/finite_size_18_27_36pilot/reference27_results.json
"""
import csv, json, pathlib, subprocess, sys, tempfile, time
import numpy as np
from symmetry import Cluster

HERE = pathlib.Path(__file__).resolve().parent
REF = HERE.parent.parent / "data" / "workspaces" / "finite_size_18_27_36pilot"
EXE = HERE / "ed_sym"


def lanczos_E0(spec, r, vec=True):
    cmd = [str(EXE), str(spec), "--r", repr(r), "--tol", "1e-9", "--maxit", "500"] + ([] if vec else ["--novec"])
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    d = {}
    for l in out.splitlines():
        if l.startswith("PASS1"):
            d["E0"] = float(l.split()[6]); d["iters"] = int(l.split()[2])
        if l.startswith("GROUND"):
            d["E0"] = float(l.split()[2]); d["residual"] = float(l.split()[4])
        if l.startswith("BASIS"):
            d["dim"] = int(l.split()[2])
    return d


def kkey(k):
    return tuple(round(x % 1.0, 6) % 1.0 for x in k)


def sector_mins(c, nup, r, vec=True):
    per_k = {}
    with tempfile.TemporaryDirectory() as td:
        for i, s in enumerate(c.sectors(nup)):
            p = pathlib.Path(td) / f"s{i}.spec"
            c.write_spec(s, nup, p)
            d = lanczos_E0(p, r, vec)
            assert d["dim"] == c.burnside_dim(s, nup)
            kk = kkey(s["k"])
            per_k.setdefault(kk, []).append(d)
    return per_k


def main():
    results = []
    ref18 = list(csv.DictReader(open(REF / "sector_completeness_18.csv")))
    gs18 = list(csv.DictReader(open(REF / "gs_18_scan.csv")))
    for shape in ("K18", "C18"):
        c = Cluster(shape)
        for r in ("0.1", "0.3", "0.5", "1.0"):
            t = time.time()
            mine = sector_mins(c, 9, float(r))
            refk = {}
            for row in ref18:
                if row["shape"] == shape and row["r"] == r:
                    kk = kkey((float(row["k1"]), float(row["k2"])))
                    refk[kk] = min(refk.get(kk, 1e9), float(row["energy"]))
            diffs = []
            for kk, e in refk.items():
                mk = kk if kk in mine else kkey((-kk[0], -kk[1]))
                diffs.append(abs(min(d["E0"] for d in mine[mk]) - e))
            gsref = next(float(g["energy"]) for g in gs18 if g["shape"] == shape and g["r"] == r)
            E0 = min(d["E0"] for v in mine.values() for d in v)
            res = dict(cluster=shape, r=float(r), momenta_compared=len(diffs), max_abs_diff_per_k=max(diffs),
                       E0=E0, E0_ref=gsref, E0_diff=abs(E0 - gsref),
                       max_residual=max(d.get("residual", 0) for v in mine.values() for d in v), seconds=time.time() - t)
            res["pass"] = res["max_abs_diff_per_k"] < 1e-9 and res["E0_diff"] < 1e-9 and res["max_residual"] < 1e-8
            print(json.dumps(res), flush=True)
            results.append(res)
    # 27 spins
    ref27 = json.load(open(REF / "reference27_results.json"))
    c = Cluster("N27")
    for rr in ref27[:2]:
        t = time.time()
        mine = sector_mins(c, 14, rr["r"], vec=False)
        E0 = min(d["E0"] for v in mine.values() for d in v)
        kref = kkey((rr["k1"], rr["k2"]))
        kmin = min(mine, key=lambda kk: min(d["E0"] for d in mine[kk]))
        res = dict(cluster="N27", r=rr["r"], E0=E0, E0_ref=rr["energy"], E0_diff=abs(E0 - rr["energy"]),
                   gs_k=list(kmin), ref_k=list(kref), seconds=time.time() - t)
        res["k_match"] = kmin in (kref, kkey((-kref[0], -kref[1]))) or \
            abs(min(d["E0"] for d in mine.get(kref, mine.get(kkey((-kref[0], -kref[1])), [{"E0": 1e9}]))) - E0) < 1e-9
        res["pass"] = res["E0_diff"] < 1e-9 and res["k_match"]
        print(json.dumps(res), flush=True)
        results.append(res)
    out = HERE / "results" / "validation_reference_18_27.json"
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(results, indent=1))
    sys.exit(0 if all(r["pass"] for r in results) else 1)


if __name__ == "__main__":
    main()
