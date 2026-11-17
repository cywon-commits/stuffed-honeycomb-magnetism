import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,itertools,json,numpy as np,scipy.linalg as la
from ed import Model
from observables import full_samples,cell_q
from local_geometry import faces,matrices,verify
from native import seeds
ROOT=pathlib.Path(__file__).resolve().parent;verify();L=3;up=2;src=Model(L=L,up=up,primitive=False,jp=.7);src.sector([0,0]);rng=np.random.default_rng(2026);v=rng.normal(size=src.nb)+1j*rng.normal(size=src.nb);v/=la.norm(v);allm=np.array([sum(1<<a for a in c) for c in itertools.combinations(range(27),up)],np.uint32);lookup={int(s):i for i,s in enumerate(allm)};psi=full_samples(src,v,allm);f=faces(L)[0];ops=matrices();errors=[]
for k in [[0,0],[0,1/3]]:
 target=Model(L=L,up=up,primitive=False,jp=.7);target.sector(cell_q(k));out,rho=seeds(src,target,v,k,rdm=k==[0,0]);direct=np.zeros((18,len(allm)),complex)
 for R in range(9):
  x,y=divmod(R,L);weight=np.exp(-2j*np.pi*(k[0]*x+k[1]*y))/L
  for mu in range(6):
   sites=f[R*6+mu];bits=sum(1<<int(a) for a in sites)
   for i,state in enumerate(allm):
    col=sum(((int(state)>>int(a))&1)<<b for b,a in enumerate(sites))
    for row in range(8):
     if row.bit_count()!=col.bit_count():continue
     t=(int(state)&~bits)|sum(((row>>b)&1)<<int(a) for b,a in enumerate(sites));ii=lookup[t]
     for kind in range(3):direct[kind*6+mu,ii]+=weight*ops[mu,kind,row,col]*psi[i]
 expected=direct[:,[lookup[int(t)] for t in target.masks]]*np.sqrt(target.lengths.astype(float));err=float(la.norm(out-expected));fullerr=float(la.norm(np.array([full_samples(target,row,allm) for row in out])-direct));assert fullerr<1e-9;gramerr=float(la.norm(out.conj()@out.T-direct.conj()@direct.T));errors.append({'k':k,'seed_error':err,'full_basis_reconstruction_error':fullerr,'Gram_error':gramerr});assert max(err,gramerr)<1e-9
 if k==[0,0]:
  rhos=[]
  for mu in range(6):
   sites=f[mu];bits=sum(1<<int(a) for a in sites);amps={}
   for i,t in enumerate(allm):
    env=int(t)&~bits;local=sum(((int(t)>>int(a))&1)<<b for b,a in enumerate(sites));amps.setdefault(env,np.zeros(8,complex))[local]=psi[i]
   exact=sum((np.outer(a,a.conj()) for a in amps.values()),np.zeros((8,8),complex));rhos.append(float(la.norm(rho[mu]-exact)));assert rhos[-1]<1e-10
(ROOT/'fast_native_validation.json').write_text(json.dumps({'status':'PASS','independent_full_27spin_two_up_seed_and_covariance':errors,'RDM_errors':rhos},indent=2));print(errors,rhos)
