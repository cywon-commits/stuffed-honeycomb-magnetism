"""Independent full-351-state operator/H/moment test of physical unfolding path."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools,json,numpy as np,scipy.linalg as la
from experiment import geometry_weights,ROOT
from ed import Model
from native import seeds
from local_geometry import faces,matrices
from observables import full_samples
L=3;up=2;s=Model(up=up);s.sector([1/3,2/3]);rng=np.random.default_rng(212);v=rng.normal(size=s.nb)+1j*rng.normal(size=s.nb);v/=la.norm(v)
states=np.array([sum(1<<i for i in x) for x in itertools.combinations(range(27),up)],np.uint32);ix={int(x):i for i,x in enumerate(states)};psi=full_samples(s,v,states);H=np.zeros((len(states),len(states)),complex)
for col,state in enumerate(states):
 for a,b,j in s.bonds:
  opposite=((int(state)>>a)^(int(state)>>b))&1;H[col,col]+=j*(.25-.5*opposite)
  if opposite:H[ix[int(state)^((1<<a)|(1<<b))],col]+=.5*j
f=faces()[0];ops=matrices();report=[]
for Q in [[0,0],[8/9,7/9],[1/3,2/3]]:
 Q=np.array(Q);k=np.mod([Q[0]+Q[1],2*Q[0]-Q[1]],1);t=Model(up=up);t.sector((s.q+Q)%1);out,_=seeds(s,t,v,k);B=out.T@geometry_weights(Q);direct=np.zeros((len(states),18),complex)
 for R in range(9):
  x,y=divmod(R,3);weight=np.exp(-2j*np.pi*(k@[x,y]))/3
  for mu in range(6):
   sites=f[6*R+mu];bits=sum(1<<int(a) for a in sites)
   for col,state in enumerate(states):
    pattern=sum(((int(state)>>int(a))&1)<<i for i,a in enumerate(sites))
    for row in range(8):
     if row.bit_count()!=pattern.bit_count():continue
     new=(int(state)&~bits)|sum(((row>>i)&1)<<int(a) for i,a in enumerate(sites))
     for kind in range(3):direct[ix[new],kind*6+mu]+=weight*ops[mu,kind,row,pattern]*psi[col]
 direct=direct@geometry_weights(Q);reconstructed=np.column_stack([full_samples(t,b,states) for b in B.T]);HB=np.column_stack([t.mv(b) for b in B.T]);e,U=la.eigh(H);a=U.conj().T@direct
 errors={'Q':Q.tolist(),'seed':float(la.norm(reconstructed-direct)),'H_seed':float(la.norm(np.column_stack([full_samples(t,b,states) for b in HB.T])-H@direct)),'moments':[]}
 for p in range(3):
  exact=direct.conj().T@np.linalg.matrix_power(H,p)@direct;spectral=a.conj().T@(e[:,None]**p*a);errors['moments'].append(float(la.norm(exact-spectral)))
 assert max(errors['seed'],errors['H_seed'],*errors['moments'])<1e-10;report.append(errors)
(ROOT/'primitive_dynamic_validation.json').write_text(json.dumps({'status':'PASS','tests':report},indent=2));print(report)
