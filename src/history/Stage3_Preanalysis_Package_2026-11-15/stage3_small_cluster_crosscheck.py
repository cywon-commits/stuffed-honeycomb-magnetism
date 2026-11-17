import math
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

L=2;N=3*L*L;UP=6
site=lambda x,y,s:3*((x%L)*L+y%L)+s
edges=[]
for x in range(L):
 for y in range(L):
  a=site(x,y,0)
  edges += [(a,site(xx,yy,1),'J') for xx,yy in [(x,y),(x-1,y),(x,y-1)]]
  edges += [(site(x,y,2),site(xx,yy,s),'P') for xx,yy,s in [(x,y,0),(x,y,1),(x+1,y,0),(x+1,y-1,1),(x+1,y-1,0),(x,y-1,1)]]
assert len(set(tuple(sorted((a,b))) for a,b,_ in edges))==len(edges)==36
states=np.array([s for s in range(1<<N) if s.bit_count()==UP],dtype=np.int32)
idx={int(s):i for i,s in enumerate(states)}

for jp in (.2,.5):
 I=[];J=[];V=[];D=np.zeros(len(states))
 for a,b,typ in edges:
  coupling=1. if typ=='J' else jp
  mask=(1<<a)|(1<<b)
  for si,s in enumerate(states):
   ai=(int(s)>>a)&1;bi=(int(s)>>b)&1
   D[si]+=coupling*(.25 if ai==bi else -.25)
   if ai!=bi:
    I.append(si);J.append(idx[int(s)^mask]);V.append(coupling*.5)
 I+=list(range(len(states)));J+=list(range(len(states)));V+=list(D)
 H=coo_matrix((V,(I,J)),shape=(len(states),len(states))).tocsr()
 E,vec=eigsh(H,k=1,which='SA',tol=1e-12)
 psi=vec[:,0];e=E[0]
 def dot_corr(a,b):
  ai=(states>>a)&1;bi=(states>>b)&1
  diag=(psi.conj()*psi*(ai-.5)*(bi-.5)).sum()
  flipped=(states^((1<<a)|(1<<b)))
  mix=np.where(ai!=bi)[0]
  off=.5*np.vdot(psi[mix],np.array([psi[idx[int(t)]] for t in flipped[mix]]))
  return (diag+off).real
 c=np.zeros(N)
 for x in range(L):
  for y in range(L):
   c[site(x,y,2)]=(-1)**x/np.sqrt(L*L)
 opvalues=np.sum(((states[:,None]>>np.arange(N))&1-.5)*c[None,:],axis=1) if False else np.array([sum(c[ii]*(((int(s)>>ii)&1)-.5) for ii in range(N)) for s in states])
 w=psi*opvalues
 direct=np.vdot(w,(H@w)-e*w).real
 formula=-.5*sum((1. if typ=='J' else jp)*(dot_corr(a,b)-float(np.vdot(psi,psi*((((states>>a)&1)-.5)*((((states>>b)&1)-.5)))).real))*(c[a]-c[b])**2 for a,b,typ in edges)
 # Check total spin S2 by S_tot^2=N*3/4+2 sum Si.Sj
 s2=N*.75+2*sum(dot_corr(i,j) for i in range(N) for j in range(i+1,N))
 print(f'L=2 12-spin jp={jp}: E={e:.10f}, S2={s2:.8f}, longitudinal moment direct={direct:.12f}, from static f-sum={formula:.12f}, |diff|={abs(direct-formula):.2e}')
