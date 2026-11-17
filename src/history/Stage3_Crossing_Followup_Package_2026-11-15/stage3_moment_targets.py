#!/usr/bin/env python3
"""Stage 3 preanalysis: 27-spin exact equal-time & symmetrized f-sum matrices.
Source is Codex crossing archive; no dynamical spectral poles are claimed.
"""
import zipfile,io,json,csv
import numpy as np
from scipy.linalg import eigh
L=3; NC=9; N=27
ROOT='hub_ground_state_crossings/'
Z='/mnt/data/hub_ground_state_crossings_verified.zip'
OUT='/mnt/data/stage3_targets_moments.csv'
POINTS=[('0.2375','p0.237500000_n12_r0_e0'),('0.25625','p0.256250000_n12_r0_e0'),('0.2609375','p0.260937500_n00_r0_e0'),('0.275','p0.275000000_n00_r0_e0'),('0.35','p0.350000000_n00_r0_e0'),('0.3875','p0.387500000_n00_r0_e0'),('0.3921875','p0.392187500_n00_r3_e0'),('0.425','p0.425000000_n00_r3_e0')]
K_MOMENTA=[(0.,0.),(0.,1/3),(1/3,2/3)]

def sid(x,y,s): return 3*((x%L)*L+(y%L))+s

def bonds(jp):
 out=[]
 for x in range(L):
  for y in range(L):
   a=sid(x,y,0)
   for xx,yy in [(x,y),(x-1,y),(x,y-1)]:out.append((a,sid(xx,yy,1),1.))
   for xx,yy,s in [(x,y,0),(x,y,1),(x+1,y,0),(x+1,y-1,1),(x+1,y-1,0),(x,y-1,1)]:out.append((sid(x,y,2),sid(xx,yy,s),jp))
 assert len(out)==81
 return out

def op(k,kind):
 vec={'n':[1,-1,0],'m':[1,1,0],'h':[0,0,1]}[kind]
 return np.array([np.exp(-2j*np.pi*(k[0]*x+k[1]*y))*vec[s]/np.sqrt(NC) for x in range(L) for y in range(L) for s in range(3)])

def first_moment(G,b,c,d):
 # Sum of spin components, symmetrized f-sum; <O_dagger(H-E0)O> when symmetry-averaged source.
 return -sum(j*G[u,v]*np.conj(c[u]-c[v])*(d[u]-d[v]) for u,v,j in b)

def scalar(s):return float(np.real(s)) if abs(np.imag(s))<1e-11 else complex(s)

rows=[];summary=[]
with zipfile.ZipFile(Z) as z:
 for jstr,fn in POINTS:
  jp=float(jstr)
  obj=json.loads(z.read(ROOT+f'obs_{fn}.json'))['summary'];S=obj['spin'];E=obj['energy']
  dd=np.load(io.BytesIO(z.read(ROOT+f'obs_{fn}_realspace.npz')))
  Gdot=dd['G_pointgroup'];Gzz=dd['Gzz']; bs=bonds(jp)
  ene=sum(j*Gdot[a,b] for a,b,j in bs)
  assert abs(ene-E)<5e-9,(ene,E)
  results={}
  for k in K_MOMENTA:
   names=['n','m','h']; v=np.array([op(k,name) for name in names]);
   M0=np.array([[np.vdot(v[a],Gdot@v[b]) for b in range(3)] for a in range(3)])
   M1=np.array([[first_moment(Gdot,bs,v[a],v[b]) for b in range(3)] for a in range(3)])
   M1zz=np.array([[-.5*sum(j*(Gdot[u,w]-Gzz[u,w])*np.conj(v[a,u]-v[a,w])*(v[b,u]-v[b,w]) for u,w,j in bs) for b in range(3)] for a in range(3)])
   M0zz=np.array([[np.vdot(v[a],Gzz@v[b]) for b in range(3)] for a in range(3)])
   assert np.max(abs(M0-M0.conj().T))<1e-9
   assert np.max(abs(M1-M1.conj().T))<1e-9
   assert abs(M1[2,2].imag)<1e-11
   assert np.linalg.eigvalsh(M0)[0]>-1e-8
   assert np.linalg.eigvalsh(M1)[0]>-1e-8
   assert abs(M1[2,2]-M1[1,1])<1e-9 if k==(0.,0.) else True
   elastic=np.zeros((3,3))
   if k==(0.,0.):
    assert abs(M0[0,1])<1e-10 and abs(M0[0,2])<1e-10
    # For a unique orbital ground multiplet, S-weight of h, m is given by projection onto S_tot.
    ah=float(np.real(NC*(M0[2,2]+M0[1,2])/(S*(S+1))))
    am=1-ah
    av=np.array([0.,am,ah]);elastic=np.outer(av,av)*S*(S+1)/NC
    assert abs(M1[1,1]+M1[1,2])<2e-10
    assert abs(M1[2,2]+M1[2,1])<2e-10
    assert abs(M1[0,1])<2e-9 and abs(M1[0,2])<2e-9
   M0inel=M0-elastic
   for a,an in enumerate(names):
    for b,bn in enumerate(names):
     r=dict(Jprime_over_J=jp,source=fn,total_spin=S,source_cell_k=str(tuple(np.asarray(obj['cell_n'])/3)),probe_k=str(k),channel=an+bn,M0_re=M0[a,b].real,M0_im=M0[a,b].imag,M1_dot_re=M1[a,b].real,M1_dot_im=M1[a,b].imag,M0_elastic_re=elastic[a,b].real,M0_inelastic_re=M0inel[a,b].real,M0_inelastic_im=M0inel[a,b].imag,M0zz_re=M0zz[a,b].real,M0zz_im=M0zz[a,b].imag,M1zz_sym_re=M1zz[a,b].real,M1zz_sym_im=M1zz[a,b].imag)
     rows.append(r)
   print('\nJp=',jp,'k=',k,'S=',S,'E=',E,'M0_nn,m,h',np.real(np.diag(M0)),'M1_nn,m,h',np.real(np.diag(M1)))
   print('M0 nh',M0[0,2],' M1 nh',M1[0,2]);print('M0 mh',M0[1,2],' M1 mh',M1[1,2]);print(' M0_zz diag',np.real(np.diag(M0zz)),' M1_zz diag',np.real(np.diag(M1zz)))
   print(' centroids n,m,h',np.real(np.diag(M1))/np.real(np.diag(M0inel)))
   if k==(0.,0.): print('Gamma elastic',np.real(np.diag(elastic)),'inel',np.real(np.diag(M0inel)))
   else:
    # Conservative moment-based multi-channel bound on some accessible mode
    A=np.linalg.eigvalsh(M0);  B=np.linalg.eigvalsh(M1)
    vals=eigh(M1,M0,eigvals_only=True)
    print('gen eig M1/M0',vals)
    if k==(0.,1/3):
     n0=float(np.real(M0[0,0]));h0=float(np.real(M0[2,2]));nh=abs(M0[0,2]);print(' integrated overlap / sqrt(M0nn*M0hh)',nh/(n0*h0)**.5)
   results[k]={'M0':M0,'M1':M1,'elastic':elastic}
  # full scalar elastic fraction at gamma h
  gamma=results[(0.,0.)];h0=gamma['M0'][2,2].real; hel=gamma['elastic'][2,2]
  summary.append({'Jp':jp,'S':S,'E':E,'rim_mN2':obj['rim_neel_m2'],'hub_S_Gamma':obj['hub_S_Gamma'],'Gamma_h_elastic_fraction':hel/h0,'Gamma_n_M0':float(gamma['M0'][0,0].real),'Gamma_n_M1':float(gamma['M1'][0,0].real),'Gamma_n_centroid':float((gamma['M1'][0,0]/gamma['M0'][0,0]).real), 'q01_h_M0':float(results[(0.,1/3)]['M0'][2,2].real),'q01_h_M1':float(results[(0.,1/3)]['M1'][2,2].real),'q01_h_centroid':float((results[(0.,1/3)]['M1'][2,2]/results[(0.,1/3)]['M0'][2,2]).real)})
with open(OUT,'w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
with open('/mnt/data/stage3_targets_summary.csv','w',newline='') as f:
 w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
print('\nSUMMARY:')
for x in summary:print(x)
print('\nWrote',OUT)
