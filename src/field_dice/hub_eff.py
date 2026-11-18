"""Effective hub-hub couplings from second-order perturbation in J' around the Neel
honeycomb, with the rim treated in linear spin-wave theory (LSWT), on an L x L torus.

H_eff = sum_{R<R'} [ jxy (T^x T^x + T^y T^y) + jz T^z T^z ] (J'^2/J units, J = J' = 1)
jxy(R,R') = -2 Re K^xx_{RR'},  K^ab_{RR'} = sum_n <0|O^a_R|n><n|O^b_R'|0> / E_n,
O_R = sum of the six rim spins of hexagon R.  Transverse part: one-magnon states;
longitudinal part (O^z = sum sigma_i (S - n_i)): two-magnon states.
A tiny staggered field h_st regularizes the Goldstone zero mode.
"""
import numpy as np, sys

S = 0.5


def build(L, h_st):
    N = 2 * L * L
    idx = lambda x, y, s: 2 * ((x % L) * L + (y % L)) + s
    sig = np.array([1 if i % 2 == 0 else -1 for i in range(N)])
    bonds = []
    for x in range(L):
        for y in range(L):
            a = idx(x, y, 0)
            bonds += [(a, idx(x, y, 1)), (a, idx(x - 1, y, 1)), (a, idx(x, y - 1, 1))]
    A = np.zeros((N, N)); B = np.zeros((N, N))
    for i, j in bonds:                      # antiparallel AF bond: diagonal + pairing
        A[i, i] += S; A[j, j] += S; B[i, j] += S; B[j, i] += S
    A += h_st * np.eye(N)
    hexes = [[idx(x, y, 0), idx(x, y, 1), idx(x + 1, y, 0), idx(x + 1, y - 1, 1), idx(x + 1, y - 1, 0), idx(x, y - 1, 1)]
             for x in range(L) for y in range(L)]
    return N, A, B, sig, hexes


def bogoliubov(A, B):
    N = len(A); M = np.block([[A, B], [B, A]])
    K = np.linalg.cholesky(M).conj().T                  # M = K^+ K
    eta = np.diag(np.r_[np.ones(N), -np.ones(N)])
    w, U = np.linalg.eigh(K @ eta @ K.conj().T)
    o = np.argsort(-w); w = w[o]; U = U[:, o]
    om = w[:N]                                          # positive branch
    T1 = np.linalg.solve(K, U[:, :N] * np.sqrt(om))     # columns for gamma_n (2N x N)
    Uu, Vv = T1[:N], T1[N:]
    T = np.block([[Uu, Vv.conj()], [Vv, Uu.conj()]])    # Psi = T Phi, Phi = (gamma, gamma^+)
    return om, T


def couplings(L, h_st=1e-4):
    N, A, B, sig, hexes = build(L, h_st)
    om, T = bogoliubov(A, B)
    r = np.sqrt(S / 2)
    def amp(R, comp):                                    # <n|O^comp_R|0>, length N
        w = np.zeros(2 * N, complex)
        for i in hexes[R]:
            if comp == 'x':
                w[i] += r; w[N + i] += r
            else:                                        # lab S^y; B local frame flips y
                s = 1 if sig[i] > 0 else -1
                w[i] += -1j * r * s; w[N + i] += 1j * r * s
        return (w @ T)[N:]
    def pair(R):                                         # symmetrized pair-creation matrix of O^z
        P = np.zeros((N, N), complex)
        for i in hexes[R]:
            P += -sig[i] * np.outer(T[i, :N].conj(), T[i, N:])
        return P + P.T
    den = om[:, None] + om[None, :]
    a0x, a0y, P0 = amp(0, 'x'), amp(0, 'y'), pair(0)
    iu = np.triu_indices(N, 1)
    out = {}
    for R in range(1, L * L):
        x, y = divmod(R, L)
        dx, dy = (x + L // 2) % L - L // 2, (y + L // 2) % L - L // 2
        d2 = dx * dx + dy * dy + dx * dy                 # 60-degree metric
        Kx = np.sum(a0x.conj() * amp(R, 'x') / om)
        Ky = np.sum(a0y.conj() * amp(R, 'y') / om)
        PR = pair(R)
        Kz = np.sum((P0.conj() * PR)[iu] / den[iu]) + np.sum(np.diag(P0).conj() * np.diag(PR) / (2 * 2 * om))
        out.setdefault(d2, []).append((-2 * Kx.real, -2 * Ky.real, -2 * Kz.real))
    return om, {d2: np.mean(v, 0) for d2, v in sorted(out.items())}


if __name__ == "__main__":
    for L in [int(a) for a in sys.argv[1:]] or [9, 12]:
        om, c = couplings(L)
        print(f"L={L}: lowest magnon {om.min():.2e}")
        for k, (d2, v) in enumerate(list(c.items())[:5]):
            print(f"   shell d^2={d2:2d}: jxx={v[0]:+.6f} jyy={v[1]:+.6f} jz={v[2]:+.6f}")
