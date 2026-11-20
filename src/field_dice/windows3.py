import numpy as np
from simplex import couplings, eps_max
from lswt3 import window
rows=[]
print("phi eps  J(AB,BC,CA)          down  window[lo,hi]   guesses")
for phd in [0,10,20,30,40,50,60]:
    ph=np.radians(phd)
    for e in [0.1,0.2,0.3,0.4]:
        if e>=eps_max(ph): continue
        J=couplings(e,ph); jab,jbc,jca=J
        # selected UUD: down = sublattice opposite the weakest bond
        weakest=int(np.argmin(J))           # 0:AB->C down, 1:BC->A down, 2:CA->B down
        sig=[(1,1,-1),(-1,1,1),(1,-1,1)][weakest]
        w=window(sig,J)
        print(f"{phd:3d} {e:.1f} ({jab:.3f},{jbc:.3f},{jca:.3f}) {'CAB'[weakest]}  {w if w is None else (round(w[0],4),round(w[1],4))}")
