"""Full S^z scan for K18 and C18: lowest energy per sector, plus plateau-sector
(S^z = N/6) levels with C2 parity, hub magnetization and (at J=J') t-eigenvalues."""
import csv, json, time, numpy as np
from ed_sz import geometry, Sector, lattice_maps
OUT = '../../data/field_dice/'
pts = sorted({round(x,4) for x in list(np.arange(0.1,2.0001,0.05))+list(np.arange(0.9,1.1001,0.01))})
grid = [(1.0, r) for r in pts] + [(round(x,4), 1.0) for x in np.arange(0.0,0.5001,0.05)]
shapes = geometry.catalog()
rows, lev = [], []
for name in ['K18','C18']:
    g = shapes[name]; N=g['N']; nplat = N//2 + N//6
    secs = {n: Sector(g, n) for n in range(N//2, N+1)}
    inv, sh = lattice_maps(g)
    Pinv = secs[nplat].perm_matrix(inv); Pt = secs[nplat].perm_matrix(sh)
    hub = [i for i in range(N) if i%3==2]
    s = secs[nplat].s
    mH = sum(((s>>i)&1) - 0.5 for i in hub)/len(hub)
    t0=time.time()
    for J,Jp in grid:
        E = {n-N//2: secs[n].lowest(J,Jp)[0][0] for n in secs if n!=nplat}
        w,v = secs[nplat].lowest(J,Jp,k=6); E[nplat-N//2]=w[0]
        rows.append(dict(shape=name,J=J,Jp=Jp,**{f'E{m}':E[m] for m in sorted(E)}))
        # quantum numbers within (near-)degenerate groups
        i=0
        while i<len(w):
            j=i
            while j+1<len(w) and w[j+1]-w[i]<1e-7: j+=1
            V=v[:,i:j+1]
            pe=np.linalg.eigvals(V.conj().T@(Pinv@V))
            te=np.linalg.eigvals(V.conj().T@(Pt@V)) if abs(J-Jp)<1e-12 else []
            mh=np.real(np.einsum('ij,i,ij->',V.conj(),mH,V))/(j-i+1)
            for a in range(i,j+1):
                lev.append(dict(shape=name,J=J,Jp=Jp,level=a,E=w[a],deg=j-i+1,
                    parities=json.dumps(sorted(np.round(pe.real,6).tolist())),
                    t_eigs=json.dumps([[round(z.real,5),round(z.imag,5)] for z in te]),
                    mH_avg=mh))
            i=j+1
    print(name,'done',round(time.time()-t0,1),'s',flush=True)
for fn,data in [('sz_scan_18.csv',rows),('plateau_levels_18.csv',lev)]:
    with open(OUT+fn,'w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=data[0].keys()); wr.writeheader(); wr.writerows(data)
