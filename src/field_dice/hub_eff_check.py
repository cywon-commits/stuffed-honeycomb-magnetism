import numpy as np, sys
from hub_eff import build, bogoliubov, S
def run(L, h_st, cut):
    N, A, B, sig, hexes = build(L, h_st); om, T = bogoliubov(A, B); r = np.sqrt(S/2)
    keep = om > cut
    def ax(R):
        w = np.zeros(2*N, complex)
        for i in hexes[R]: w[i] += r; w[N+i] += r
        return (w @ T)[N:]
    def pair(R):
        P = np.zeros((N, N), complex)
        for i in hexes[R]: P += -sig[i]*np.outer(T[i,:N].conj(), T[i,N:])
        return P + P.T
    a0 = ax(0); P0 = pair(0); den = om[:,None]+om[None,:]
    mask = np.outer(keep, keep); iu = np.triu_indices(N, 1)
    Ksum = 0; res = {}
    for R in range(L*L):
        aR = ax(R); Kx = np.sum((a0.conj()*aR/om)[keep]); Ksum += Kx.real
        x, y = divmod(R, L); dx, dy = (x+L//2)%L-L//2, (y+L//2)%L-L//2; d2 = dx*dx+dy*dy+dx*dy
        if d2 in (1, 3, 4) and d2 not in res:
            PR = pair(R); M = (P0.conj()*PR/den)*mask
            Kz = np.sum(M[iu]) + np.sum(np.diag(M))/2
            res[d2] = (-2*Kx.real, -2*Kz.real)
    return Ksum, res, om.min()
for L in [int(a) for a in sys.argv[1:]]:
    for h, cut in [(1e-4, 0), (1e-4, 0.05), (1e-6, 0.05)]:
        Ks, res, m = run(L, h, cut)
        print(f"L={L} h_st={h:g} drop_omega<{cut}: sum_R K^xx={Ks:.4f} (sum rule 9*chi_perp=1.5) | " +
              " ".join(f"d2={d}: jxy={v[0]:+.4f} jz={v[1]:+.4f}" for d, v in sorted(res.items())), flush=True)
