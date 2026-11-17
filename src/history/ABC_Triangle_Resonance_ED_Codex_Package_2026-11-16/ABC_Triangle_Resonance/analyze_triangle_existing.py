#!/usr/bin/env python3
"""Verified two-point-derived triangle quantities, N=27 triangular-stuffed-honeycomb ED.

No assumptions about higher-point correlations are made. Requires previously
supplied stage1/2 and level-crossing ZIPs. All quantities computed from their
pointgroup-averaged 27x27 scalar spin correlation matrices.
"""
import argparse,csv,io,re,zipfile,sys
from pathlib import Path
import numpy as np

L=3

def sid(x,y,s):return 3*((x%L)*L+y%L)+s

def triangle_list():
    # Geometry as in original ed.geometry(): 27 AB bonds and 54 hub-rim bonds.
    bonds={}
    for x in range(L):
      for y in range(L):
        a=sid(x,y,0)
        for xx,yy in [(x,y),(x-1,y),(x,y-1)]:
            b=sid(xx,yy,1);bonds[tuple(sorted((a,b)))]=1.0
        h=sid(x,y,2)
        for xx,yy,s in [(x,y,0),(x,y,1),(x+1,y,0),
                        (x+1,y-1,1),(x+1,y-1,0),(x,y-1,1)]:
            i=sid(xx,yy,s);bonds[tuple(sorted((h,i)))]=2.0
    triangles=[]
    for x in range(L):
      for y in range(L):
        h=sid(x,y,2)
        A=[sid(xx,yy,0) for xx,yy in [(x,y),(x+1,y),(x+1,y-1)]]
        B=[sid(xx,yy,1) for xx,yy in [(x,y),(x+1,y-1),(x,y-1)]]
        for a in A:
          for b in B:
            if bonds.get(tuple(sorted((a,b))))==1.0:
              assert bonds.get(tuple(sorted((a,h))))==2.0
              assert bonds.get(tuple(sorted((b,h))))==2.0
              triangles.append((a,b,h))
    assert len(bonds)==81, len(bonds)
    assert len(triangles)==54,len(triangles)
    assert len(set(tuple(sorted(t)) for t in triangles))==54
    return triangles,bonds

TRIANGLES,BONDS=triangle_list()

def read_npz_from_zip(zpath,member):
    if zpath.is_file():
        with zipfile.ZipFile(zpath) as z:
          with z.open(member) as f:
            raw=f.read()
          with np.load(io.BytesIO(raw)) as data:
            G=data['G_pointgroup'].copy()
    else:
        local=Path(__file__).resolve().parent/'reference_realspace'/Path(member).name
        with np.load(local) as data:
            G=data['G_pointgroup'].copy()
    assert G.shape==(27,27)
    assert np.max(np.abs(G-G.T))<1e-8
    assert np.max(abs(np.diag(G)-.75))<1e-7
    return G

def source_list(root):
    s_zip=root/'hub_static_stages12_verified.zip'
    c_zip=root/'hub_ground_state_crossings_verified.zip'
    out=[]
    for j,label in [(0.2,'j0.2'),(0.5,'j0.5'),(1.0,'j1')]:
      out.append((j,'ground',s_zip,'hub_static_stages12/static_'+label+'_realspace.npz'))
    jobs=[('0.237500000','n12_r0'),('0.256250000','n12_r0'),
          ('0.260937500','n00_r0'),('0.275000000','n00_r0'),
          ('0.350000000','n00_r0'),('0.387500000','n00_r0'),
          ('0.392187500','n00_r3'),('0.425000000','n00_r3')]
    for s,r in jobs:
      out.append((float(s),'ground',c_zip,f'hub_ground_state_crossings/obs_p{s}_{r}_e0_realspace.npz'))
    # At each crossing include the competitor branch for comparing discontinuities.
    for s,r in [('0.256250000','n00_r0'),('0.387500000','n00_r3')]:
      out.append((float(s),'competitor',c_zip,f'hub_ground_state_crossings/obs_p{s}_{r}_e0_realspace.npz'))
    return sorted(out,key=lambda v:(v[0],v[1]))

def summarize(G,j,role):
    ab=np.array([G[a,b] for a,b,h in TRIANGLES]);ah=np.array([G[a,h] for a,b,h in TRIANGLES]);bh=np.array([G[b,h] for a,b,h in TRIANGLES])
    ba,bh,bc=float(ab.mean()),float(ah.mean()),float(bh.mean())
    sum_b=ba+bh+bc
    p12=0.5-2*sum_b/3 # exact spin-1/2 identity
    S2=2.25+2*sum_b
    p32=1-p12
    dimer=ba+np.exp(2j*np.pi/3)*bh+np.exp(4j*np.pi/3)*bc
    scalar_chi2=3*p12/16  # exact chi^2=3 P12/16
    sym_dimer_amp=9*p12/8 # exact 1/2 {Ddag,D}=9P12/8
    decomp_ab=len([b for b in BONDS.values() if b==1.0]);decomp_ch=len([b for b in BONDS.values() if b==2.0]);assert decomp_ab==27 and decomp_ch==54
    E=27*ba+54*j*(bh+bc)/2
    # Angle averaged result for the chirality-odd part (not available from 2pt matrices).
    return {'Jp_over_J':j,'branch':role,'AB_bond':ba,'AC_bond':bh,'BC_bond':bc,
            'P_singlet_AB':0.25-ba,'P_singlet_AC':0.25-bh,'P_singlet_BC':0.25-bc,
            'P_triangle_doublet':p12,'P_triangle_quartet':p32,
            'S_triangle_squared':S2,'D_complex_real':float(dimer.real),
            'D_complex_imag':float(dimer.imag),'D_abs':float(abs(dimer)),
            'D_fluct_sym':sym_dimer_amp - abs(dimer)**2,
            'DdagD_plus_DDdag_over2':sym_dimer_amp,
            'chi_squared':scalar_chi2,'energy_reconstructed':E,
            'max_triangle_G_spread':float(max(np.ptp(ab),np.ptp(ah),np.ptp(bh)))}

def main():
    pa=argparse.ArgumentParser();pa.add_argument('--input-dir',type=Path,default=Path('/mnt/data'))
    pa.add_argument('--out',type=Path,default=Path('/mnt/data/triangle_resonance_work/triangle_27spin_existing.csv'))
    a=pa.parse_args(); rows=[]
    for j,role,zpath,member in source_list(a.input_dir):
       g=read_npz_from_zip(zpath,member); rows.append(summarize(g,j,role))
    a.out.parent.mkdir(parents=True,exist_ok=True)
    with a.out.open('w',newline='') as f:
       w=csv.DictWriter(f,fieldnames=rows[0].keys());w.writeheader();w.writerows(rows)
    for r in rows:
       print(f"{r['Jp_over_J']:>9.6f} {r['branch']:>10}  AB={r['AB_bond']:+.7f}  AC={r['AC_bond']:+.7f}  P12={r['P_triangle_doublet']:.8f}  |<D>|={r['D_abs']:.8f}  <chi²>={r['chi_squared']:.8f}  E={r['energy_reconstructed']:.8f}  spread={r['max_triangle_G_spread']:.3g}")

if __name__=='__main__':main()
