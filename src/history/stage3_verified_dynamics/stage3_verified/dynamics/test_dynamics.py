import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,pathlib,numpy as np,scipy.linalg as la
from ed import Model,geometry
from observables import full_samples,spin_seed,cell_q,site_weights,coefficients
ROOT=pathlib.Path(__file__).resolve().parent
L=2;up=6;j=.2375;s=Model(L=L,up=up,primitive=False,jp=j);s.sector([0,0]);H=np.column_stack([s.mv(np.eye(s.nb,dtype=complex)[:,i]) for i in range(s.nb)]);e,V=la.eigh(H);v=V[:,0];k=np.array([0,.5]);t=Model(L=L,up=up,primitive=False,jp=j);t.sector(cell_q(k));Ht=np.column_stack([t.mv(np.eye(t.nb,dtype=complex)[:,i]) for i in range(t.nb)]);et,U=la.eigh(Ht);B=np.column_stack([spin_seed(s,t,v,k,c) for c in ['n','m','H']]);errors=[]
# Independent bit-diagonal longitudinal operator evaluated directly at target orbit representatives.
psi=full_samples(s,v,t.masks)
for i,c in enumerate(['n','m','H']):
 weights=site_weights(L,k,coefficients(k,c));diag=sum(weights[a]*(((t.masks>>a)&1)-.5) for a in range(t.n));exact=psi*diag*np.sqrt(t.lengths.astype(float));errors.append(float(la.norm(exact-B[:,i])))
A=U.conj().T@B;D=Ht-e[0]*np.eye(t.nb);direct=[B.conj().T@B,B.conj().T@D@B,(D@B).conj().T@(D@B)];momenterrors=[float(la.norm(A.conj().T@((et-e[0])[:,None]**p*A)-direct[p])) for p in range(3)]
assert max(errors)<1e-10 and max(momenterrors)<1e-9
result={'independent_longitudinal_seed_errors':errors,'complete_spectral_moment_errors_M0_M1_M2':momenterrors,'status':'PASS','cluster':'12 spins, L=2, Nup=6; full target-block spectral decomposition; independent bit-diagonal Fourier seeds'}
(ROOT/'small_dynamics_validation.json').write_text(json.dumps(result,indent=2));print(result)
