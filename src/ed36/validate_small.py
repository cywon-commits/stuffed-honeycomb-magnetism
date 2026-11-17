"""Independent validation of ed_sym on small clusters.

1. Builds the full (unsymmetrized) fixed-Sz Hamiltonian directly from the
   reference geometry bonds with numpy and diagonalizes it densely.
2. Exports every symmetry sector, builds each sector matrix with ed_sym
   (--mode dense), checks Hermiticity and dimension against Burnside, and
   checks that the union of all sector spectra equals the full spectrum.
"""
import json, math, pathlib, subprocess, sys, tempfile
import numpy as np
from symmetry import Cluster

HERE = pathlib.Path(__file__).resolve().parent
EXE = HERE / "ed_sym"


def full_spectrum(c, nup, r):
    states = [s for s in range(1 << c.N) if bin(s).count("1") == nup]
    idx = {s: i for i, s in enumerate(states)}
    D = len(states)
    H = np.zeros((D, D))
    for i, s in enumerate(states):
        for a, b, t in c.bonds:
            J = 1.0 if t == 0 else r
            if ((s >> a) & 1) == ((s >> b) & 1):
                H[i, i] += 0.25 * J
            else:
                H[i, i] -= 0.25 * J
                H[idx[s ^ (1 << a) ^ (1 << b)], i] += 0.5 * J
    return np.linalg.eigvalsh(H)


def sector_matrix(spec, r):
    out = subprocess.run([str(EXE), str(spec), "--mode", "dense", "--r", repr(r)], capture_output=True, text=True, check=True).stdout
    lines = out.splitlines()
    k = next(i for i, l in enumerate(lines) if l.startswith("DENSE"))
    D = int(lines[k].split()[1])
    M = np.zeros((D, D), complex)
    for l in lines[k + 1:k + 1 + D * D]:
        i, j, re, im = l.split()
        M[int(i), int(j)] = float(re) + 1j * float(im)
    return M


def main(name, r):
    c = Cluster(name)
    nup = (c.N + 1) // 2
    ref = full_spectrum(c, nup, r)
    secs = c.sectors(nup)
    allev = []
    worst_herm = 0.0
    with tempfile.TemporaryDirectory() as td:
        for i, s in enumerate(secs):
            p = pathlib.Path(td) / f"s{i}.spec"
            c.write_spec(s, nup, p)
            M = sector_matrix(p, r)
            assert M.shape[0] == c.burnside_dim(s, nup), (i, M.shape, c.burnside_dim(s, nup))
            worst_herm = max(worst_herm, np.abs(M - M.conj().T).max() if M.size else 0.0)
            allev.extend(np.linalg.eigvalsh((M + M.conj().T) / 2))
    allev = np.sort(allev)
    res = dict(cluster=name, r=r, nup=nup, full_dim=len(ref), sum_sector_dims=len(allev), sectors=len(secs),
               max_hermiticity_error=float(worst_herm),
               max_abs_spectrum_difference=float(np.abs(allev - ref).max()) if len(allev) == len(ref) else None,
               E0_full=float(ref[0]), E0_sectors=float(allev[0]))
    res["pass"] = bool(res["max_abs_spectrum_difference"] is not None and res["max_abs_spectrum_difference"] < 1e-10 and worst_herm < 1e-12)
    print(json.dumps(res))
    return res


if __name__ == "__main__":
    out = [main("N12", r) for r in (0.3, 0.77, 1.0)]
    pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else HERE / "validation_N12.json").write_text(json.dumps(out, indent=1))
    sys.exit(0 if all(o["pass"] for o in out) else 1)
