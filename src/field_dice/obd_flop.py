"""Order-by-disorder spin-flop picture of the 1/3 plateau on the honeycomb mirror (J_AB=1, J_BC=J_CA=r<1).
At h = 3 S r the classical ground states are {hub up, rim Neel vector n arbitrary}. LSWT zero-point energy
E_zp(theta) (theta = angle of n from the field) selects n || h if K = E_zp(pi/2) - E_zp(0) > 0.
The rim then stays collinear (m = S/3) until the spin-flop at |h - h0| = sqrt(2K/chi_perp)."""
import sys
import numpy as np
from plateau_1S import zero_point, classical

def neel_state(theta):
    n = np.array([np.sin(theta), 0.0, np.cos(theta)])
    return [-n, n, np.array([0, 0, 1.0])]

def Ezp(r, theta, S=0.5, n=48):
    return zero_point(neel_state(theta), (1.0, r, r), 3 * S * r, S, n)

def chi_perp(r, S=0.5, d=0.02):
    J = (1.0, r, r); h0 = 3 * S * r; ms = []
    for h in (h0 - d, h0 + d):
        sp, _ = classical(J, h, S); ms.append(S * sum(s[2] for s in sp) / 3)
    return (ms[1] - ms[0]) / (2 * d)

if __name__ == '__main__':
    S = 0.5
    for r in [float(x) for x in sys.argv[1:]]:
        th = np.linspace(0, np.pi / 2, 5)
        E = [Ezp(r, t, S)[0] for t in th]
        K = E[-1] - E[0]; chi = chi_perp(r, S)
        W = 2 * np.sqrt(2 * K / chi) if K > 0 else 0.0
        shape = [(e - E[0]) / K if K else 0 for e in E]
        print(f"r={r:.2f} K={K:.5f} chi={chi:.4f} W={W:.4f} E(theta)-E(0) / K at theta=0..90: {np.round(shape,3)}")
