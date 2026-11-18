import csv, numpy as np, json
rows=list(csv.DictReader(open('../../data/field_dice/sz_scan_18.csv')))
out=[]
for r in rows:
    J,Jp=float(r['J']),float(r['Jp']); E=np.array([float(r[f'E{m}']) for m in range(10)])
    e0=E.min(); S=max(m for m in range(10) if E[m]<e0+1e-8)
    lo=max((E[3]-E[m])/(3-m) for m in range(3)); hi=min((E[m]-E[3])/(m-3) for m in range(4,10))
    hsat=E[9]-E[8]
    sc = J if J>0 else 1.0
    out.append(dict(shape=r['shape'],J=J,Jp=Jp,S=S,h_lo=lo,h_hi=hi,width=hi-lo,h_sat=hsat,
        h_sat_pred=0.5*max(6*J+3*Jp,9*Jp)))
with open('../../data/field_dice/plateau_summary_18.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=out[0].keys()); w.writeheader(); w.writerows(out)
print("shape  J    J'   S | plateau [lo, hi] width | cl. lower,upper | h_sat ED vs pred")
for o in out:
    J,Jp=o['J'],o['Jp']
    if J==1 and round(Jp,2) not in [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,1.0,1.05,1.1,1.2,1.4,1.6,1.8,2.0]: continue
    cl = f"{1.5*Jp:.3f},{1.5*Jp:.3f}" if Jp<=J else f"{max(0,1.5*(2*J-Jp)):.3f},{1.5*Jp:.3f}"
    print(f"{o['shape']} {J:4.2f} {Jp:4.2f} {o['S']} | [{o['h_lo']:.4f},{o['h_hi']:.4f}] {o['width']:+.4f} | {cl} | {o['h_sat']:.4f} {o['h_sat_pred']:.4f}")
