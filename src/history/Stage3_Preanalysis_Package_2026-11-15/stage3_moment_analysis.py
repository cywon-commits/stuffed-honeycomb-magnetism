#!/usr/bin/env python3
"""Exact SU(2)-averaged spectral first moments from Stage 2 static ED data.

Uses 27x27 real-space SU(2) dot correlations from hub_static_stages12_verified.zip.
No numerical dynamical Lanczos or magnon linewidths are calculated.
"""
from __future__ import annotations
import csv, io, json, math, zipfile
import numpy as np

SOURCE='/mnt/data/hub_static_stages12_verified.zip'
OUT='/mnt/data/stage3_first_moments.csv'
L=3; NC=L*L; N=3*NC

def sid(x,y,s): return 3*((x%L)*L+y%L)+s

def bonds(jp):
    result=[]
    for x in range(L):
        for y in range(L):
            a=sid(x,y,0)
            result.extend((a,sid(xx,yy,1),1.) for xx,yy in [(x,y),(x-1,y),(x,y-1)])
            result.extend((sid(x,y,2),sid(xx,yy,s),jp) for xx,yy,s in [(x,y,0),(x,y,1),(x+1,y,0),(x+1,y-1,1),(x+1,y-1,0),(x,y-1,1)])
    assert len(result)==81 and len({tuple(sorted((i,j))) for i,j,w in result})==81
    return result

def op(channel,k=(0,0)):
    c=np.zeros(N,complex)
    for x in range(L):
        for y in range(L):
            phase=np.exp(-2j*np.pi*(k[0]*x+k[1]*y))/np.sqrt(NC)
            for s in range(3):
                v={'h':(0,0,1),'m':(1,1,0),'n':(1,-1,0)}[channel][s]
                c[sid(x,y,s)]=v*phase
    return c

def m0_m1(G,bs,c):
    m0=np.vdot(c,G@c).real
    # SU(2) scalar spin-dot spectral function summed over x,y,z.
    # Exact symmetrized f-sum = - sum_{<ij>} J_ij <Si.Sj> |c_i-c_j|^2.
    m1=-sum(w*G[i,j]*abs(c[i]-c[j])**2 for i,j,w in bs)
    return m0,float(m1)

rows=[]
with zipfile.ZipFile(SOURCE) as z:
    for x,S in [(0.2,0.5),(0.5,2.5),(1.0,0.5)]:
        j='j0.2' if x==0.2 else 'j0.5' if x==0.5 else 'j1'
        folder='hub_static_stages12/'
        G=np.load(io.BytesIO(z.read(folder+f'static_{j}_realspace.npz')))['G_pointgroup']
        raw=json.loads(z.read(folder+f'static_{j}.json'))['summary']
        bs=bonds(x)
        er=sum(w*G[i,jj] for i,jj,w in bs)
        assert abs(er-raw['energy'])<2e-9, (x,er,raw['energy'])
        hr=np.mean([G[i,jj] for i,jj,w in bs if w!=1.]) if x!=1. else raw['hub_rim_bond_mean']
        mu_expected=-6*x*raw['hub_rim_bond_mean']
        for ch,k,kind in [('n',(0,0),'staggered_Gamma'),('m',(0,0),'uniform_Gamma'),('h',(0,0),'hub_Gamma'),('h',(0,1/3),'hub_generic_peak_orbit'),('h',(1/3,2/3),'hub_corner')]:
            m0,m1=m0_m1(G,bs,op(ch,k))
            if ch=='h': assert abs(m1-mu_expected)<2e-12, (x,k,m1,mu_expected)
            if ch=='m': assert abs(m1-mu_expected)<2e-12
            e0=None
            if k==(0,0) and ch in ('m','h') and x!=1.0:
                # Assume only SU(2) ground multiplet in this cell sector.
                hh,_=m0_m1(G,bs,op('h',(0,0)))
                mh=np.vdot(op('m'),G@op('h')).real
                S2=S*(S+1)
                alpha_h=NC*(hh+mh)/S2
                alpha=alpha_h if ch=='h' else 1-alpha_h
                e0=alpha**2*S2/NC
            elif k==(0,0) and ch=='n':
                e0=0.0 # inversion odd; in symmetric, nondegenerate-parity multiplet
            elif ch=='h' and k==(0,1/3):
                e0=0.0 # difference does not connect ground orbit
            elif ch=='h' and k==(1/3,2/3) and x==0.5:
                e0=0.0 # unique Gamma ground sector
            # x=1 at gamma has K/K' orbital degeneracy, requires full P0
            inel=None if e0 is None else m0-e0
            centroid_all=m1/m0 if m0>0 else None
            centroid_inel=None if inel is None or inel<1e-12 else m1/inel
            rows.append(dict(Jp_J=x, channel=ch, k1=k[0],k2=k[1], M0_total=m0,M1_exact=m1,
                             M0_ground_proj=e0, M0_inelastic=inel, mean_energy_all=m1/m0,
                             mean_energy_inelastic=centroid_inel, label=kind, S=S))
            print(x,ch,k,'M0',round(m0,9),'M1',round(m1,9),'elastic',None if e0 is None else round(e0,9),'mean_inel',None if centroid_inel is None else round(centroid_inel,9))
        # exact Gamma f-sum matrix for (m,h): (mu,-mu;-mu,mu)
        m=op('m'); h=op('h')
        mu=m0_m1(G,bs,h)[1]
        assert abs(m0_m1(G,bs,m+h)[1])<1e-12
        assert abs(m0_m1(G,bs,m-h)[1]-4*mu)<1e-12
        print('Gamma conserved f-sum identity verified', x,mu)

with open(OUT,'w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print('Wrote',OUT)

# Longitudinal (fixed Sz=1/2) response: accurate Stage 3 validation moments.
zzrows=[]
with zipfile.ZipFile(SOURCE) as z:
    for x in (0.2,0.5,1.):
        j='j0.2' if x==0.2 else 'j0.5' if x==0.5 else 'j1'
        d=np.load(io.BytesIO(z.read(f'hub_static_stages12/static_{j}_realspace.npz')))
        G=d['G_selected']; Gzz=d['Gzz_selected']; mags=d['magnetization']
        b=bonds(x)
        for ch,k in [('n',(0,0)),('m',(0,0)),('h',(0,0)),('h',(0,1/3)),('h',(1/3,2/3))]:
            c=op(ch,k)
            z0=np.vdot(c,Gzz@c).real
            mean=sum(c[sid(xx,yy,s)]*mags[s] for xx in range(3) for yy in range(3) for s in range(3))
            z0_connected=z0-abs(mean)**2
            z1=-.5*sum(w*(G[i,jj]-Gzz[i,jj])*abs(c[i]-c[jj])**2 for i,jj,w in b)
            zzrows.append(dict(Jp_J=x,channel=ch,k1=k[0],k2=k[1], M0zz_raw=z0, M0zz_connected=z0_connected, M1zz_sym=z1,mean_energy_sym_no_extra_zero_proj=z1/z0_connected if z0_connected>1e-10 else None))
            print('ZZ',x,ch,k, 'M0',round(z0,9),'conn',round(z0_connected,9),'M1_sym',round(z1,9),'centroid',round(z1/z0_connected,9))
with open('/mnt/data/stage3_first_moments_longitudinal.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(zzrows[0])); w.writeheader();w.writerows(zzrows)
