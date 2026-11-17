"""Exact spatial covariance of physical-momentum operator residue matrices."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import itertools,json,numpy as np,scipy.linalg as la
from experiment import ROOT,geometry_weights,pack,dump
from ensemble import maps,unfold,fourier
from sectors import point_permutations
from ed import geometry,translations

def V(Q):
 Q=np.asarray(Q);k=np.mod([Q[0]+Q[1],2*Q[0]-Q[1]],1);W=geometry_weights(Q);v=np.zeros((162,6),complex)
 for R,(x,y) in enumerate(itertools.product(range(3),repeat=2)):
  for a in range(3):v[a*54+R*6:a*54+R*6+6,2*a:2*a+2]=np.exp(-2j*np.pi*(k@[x,y]))/3*W[a*6:a*6+6,2*a:2*a+2]
 return v

def build():
 q0=np.array([8/9,7/9]);v0=V(q0);raw=np.load(ROOT/'inputs/raw_observables_j1.npz')['selected_full_raw_Gram'];S=unfold(fourier(raw[:54,:54]),'D')[0];Qs=[r['Q'] for r in S];coord,_=geometry(3);gs,ts=translations(3,coord,True);lookup={tuple(p):i for i,p in enumerate(ts)};qK=np.array([1/3,2/3]);qsources=[qK,(-qK)%1];pg=point_permutations();swap=np.r_[np.arange(54,108),np.arange(54),np.arange(108,162)];sign=np.r_[np.ones(108),-np.ones(54)];rawprime=raw[np.ix_(swap,swap)].conj()*sign[:,None]*sign[None,:];records=[]
 for src,sourceQ in enumerate(qsources):
  C=raw if src==0 else rawprime;G0=v0.conj().T@C@v0;chars=np.exp(-2j*np.pi*(gs@sourceQ))
  for i,p in enumerate(pg):
   ix,c,_=maps(p);vp=np.zeros_like(v0);vp[ix]=c[:,None]*v0
   possibilities=[]
   for Q in Qs:
    v=V(Q);T=v.conj().T@vp;err=la.norm(v@T-vp);possibilities.append((err,Q,T))
   err,Q,T=min(possibilities,key=lambda r:r[0]);assert err<1e-10 and la.norm(T.conj().T@T-np.eye(6))<1e-10
   inv=np.argsort(p);mapped=np.array([chars[lookup[tuple(inv[t[p]])]] for t in ts]);srcnew=int(np.argmin([la.norm(mapped-np.exp(-2j*np.pi*(gs@q))) for q in qsources]));Cnew=raw if srcnew==0 else rawprime;exact=V(Q).conj().T@Cnew@V(Q);pred=T@G0@T.conj().T;ge=float(la.norm(exact-pred));assert ge<1e-7
   records.append({'initial_source':'K' if src==0 else 'Kprime','source_after_map':'K' if srcnew==0 else 'Kprime','point_group_index':i,'physical_Q':Q.tolist(),'operator_transform':pack(T),'operator_covariance_error':float(err),'static_covariance_error':ge,'rule':'R_mapped = T R_original T^dagger; frequency and direct residual unchanged by exact unitary map.'})
 dump(ROOT/'symmetry_partner_maps.json',{'records':records,'maximum_static_covariance_error':max(r['static_covariance_error'] for r in records),'scope':'Exact symmetry-derived partners, not separately run eigensolves. Includes transformation of all complex cross residues.'});print('PASS',len(records),'maps; max static error',max(r['static_covariance_error'] for r in records))
if __name__=='__main__':build()
