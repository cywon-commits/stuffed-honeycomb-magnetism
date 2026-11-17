"""Static SU(2)-scalar correlations and independently validated spin-operator seeds."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import ctypes as ct, pathlib,subprocess,time,json
import numpy as np,scipy.linalg as la
from ed import Model,ptr,P
ROOT=pathlib.Path(__file__).resolve().parent

def loadops():
 so=ROOT/'observables_kernel.so';src=ROOT/'observables.cpp'
 if not so.exists() or so.stat().st_mtime<src.stat().st_mtime:subprocess.run(['g++','-O3','-march=native','-fopenmp','-shared','-fPIC',str(src),'-o',str(so)],check=True)
 lib=ct.CDLL(str(so));lib.pair_stats.argtypes=[P]*17;lib.spin_seed.argtypes=[P]*14;return lib

def cell_q(k):return np.array([(k[0]+k[1])/3,(2*k[0]-k[1])/3])
def rank_array(m,states):
 lo=states&((1<<m.split)-1);return m.rl[lo]+m.rh[m.pop[lo],states>>m.split]
def full_samples(m,v,states):
 rr=rank_array(m,states);j=m.rawbasis[m.rid[rr]];out=np.zeros(len(states),complex);ok=j>=0
 out[ok]=v[j[ok]]*m.phase[m.shift[rr[ok]]]/np.sqrt(m.lengths[j[ok]].astype(float));return out

def convert_basis(source,target,v):
 return full_samples(source,v,target.masks)*np.sqrt(target.lengths.astype(float))

def pair_classes(m):
 pairs=[(a,b) for a in range(m.n) for b in range(a+1,m.n)]
 keys=[min(tuple(sorted((int(p[a]),int(p[b])))) for p in m.perms) for a,b in pairs]
 unique=sorted(set(keys));lookup={k:i for i,k in enumerate(unique)}
 aa=np.array([p[0] for p in pairs],np.int32);bb=np.array([p[1] for p in pairs],np.int32);cls=np.array([lookup[k] for k in keys],np.int32)
 return aa,bb,cls,np.bincount(cls),unique

def pair_expectations(m,v):
 assert len(m.gs)==m.L*m.L,'Use cell translations for sublattice-resolved observables.'
 aa,bb,cls,counts,keys=pair_classes(m);nc=len(counts)
 config=np.array([m.nb,len(aa),nc,m.n,m.L*m.L,m.split,m.rh.shape[1],m.threads],np.int32)
 dot=np.zeros(nc);zz=np.zeros(nc);mag=np.zeros(3);t=time.perf_counter()
 loadops().pair_stats(*[ptr(a) for a in [config,m.masks,m.lengths,m.rid,m.shift,m.rawbasis,m.pop,m.rl,m.rh,aa,bb,cls,m.phase,np.ascontiguousarray(v),dot,zz,mag]])
 G=np.eye(m.n)*.75;Gz=np.eye(m.n)*.25
 for a,b,c in zip(aa,bb,cls):G[a,b]=G[b,a]=dot[c]/counts[c];Gz[a,b]=Gz[b,a]=zz[c]/counts[c]
 return {'Gdot':G,'Gzz':Gz,'magnetization':mag,'seconds':time.perf_counter()-t,'pair_classes':nc,'pair_class_counts':counts.tolist()}

def site_weights(L,k,coeff):
 return np.ascontiguousarray([coeff[s]*np.exp(-2j*np.pi*(k[0]*x+k[1]*y))/L for x in range(L) for y in range(L) for s in range(3)],complex)

def coefficients(k,channel):
 if channel in ['physical_K','physical_Kprime'] and np.max(abs(np.exp(2j*np.pi*np.asarray(k))-1))>1e-10:raise ValueError('The fixed physical K/Kprime probes fold to cell Gamma; use k=(0,0).')
 a,b=k;fa=1+np.exp(2j*np.pi*a)+np.exp(2j*np.pi*(a-b));fb=1+np.exp(2j*np.pi*(a-b))+np.exp(-2j*np.pi*b)
 if channel=='A':return np.array([1,0,0],complex)
 if channel=='B':return np.array([0,1,0],complex)
 if channel=='H':return np.array([0,0,1],complex)
 if channel=='total':return np.ones(3,complex)
 if channel=='n':return np.array([1,-1,0],complex)
 if channel=='m':return np.array([1,1,0],complex)
 if channel=='Phi':return np.array([fa,fb,0])
 if channel=='physical_K':return np.exp(-2j*np.pi*np.array([0,1/3,2/3]))/np.sqrt(3)
 if channel=='physical_Kprime':return np.exp(-2j*np.pi*np.array([0,2/3,1/3]))/np.sqrt(3)
 raise ValueError(channel)

def spin_seed(source,target,v,k,channel,kind='z'):
 assert len(source.gs)==source.L**2 and len(target.gs)==target.L**2
 expected=source.up+{'z':0,'+':1,'-':-1}[kind];assert target.up==expected
 diff=target.gs@(target.q-source.q-cell_q(k));assert np.max(abs(np.exp(2j*np.pi*diff)-1))<1e-8
 w=site_weights(source.L,k,coefficients(k,channel));cfg=np.array([target.nb,source.n,source.split,source.rh.shape[1],source.threads,{'z':0,'+':1,'-':-1}[kind]],np.int32);out=np.zeros(target.nb,complex)
 loadops().spin_seed(*[ptr(a) for a in [cfg,target.masks,target.lengths,source.rid,source.shift,source.rawbasis,source.lengths,source.pop,source.rl,source.rh,source.phase,w,np.ascontiguousarray(v),out]])
 return out

def correlation_matrices(G,L):
 # Translation average explicitly; this also handles a symmetry-averaged ensemble.
 out={}
 for n1 in range(L):
  for n2 in range(L):
   k=np.array([n1,n2])/L;C=np.zeros((3,3),complex)
   for x in range(L):
    for y in range(L):
     for xp in range(L):
      for yp in range(L):
       ph=np.exp(-2j*np.pi*(k[0]*(xp-x)+k[1]*(yp-y)))/(L*L)
       C+=ph*G[3*(x*L+y):3*(x*L+y)+3,3*(xp*L+yp):3*(xp*L+yp)+3]
   out[(n1,n2)]=C
 return out

def inversion_vector(m,v,perm):
 # Inversion is its own inverse, and maps k to -k.
 lo=np.zeros(1<<m.split,np.uint32);hi=np.zeros(1<<(m.n-m.split),np.uint32)
 for table,start in [(lo,0),(hi,m.split)]:
  for mask in range(1,len(table)):
   bit=mask&-mask;table[mask]=table[mask^bit]|(1<<int(perm[start+bit.bit_length()-1]))
 pre=lo[m.masks&((1<<m.split)-1)]|hi[m.masks>>m.split]
 return full_samples(m,v,pre)*np.sqrt(m.lengths.astype(float))

def total_spin_check(m,v):
 old=(m.flips,m.js,m.diag)
 try:
  pairs=[(a,b,2.) for a in range(m.n) for b in range(a+1,m.n)]
  m.flips=np.array([(1<<a)|(1<<b) for a,b,j in pairs],np.uint32);m.js=np.full(len(pairs),2.)
  m.diag=np.full(m.nb,.75*m.n)
  for a,b,j in pairs:m.diag+=2*(.25-.5*(((m.masks>>a)^(m.masks>>b))&1))
  # mv uses len(m.bonds), so temporarily replace that list as well.
  bonds=m.bonds;m.bonds=pairs
  w=m.mv(v);ev=np.vdot(v,w).real;S=(np.sqrt(1+4*ev)-1)/2;near=round(S-.5)+.5
  return {'S2':float(ev),'spin':float(S),'nearest_half_integer_spin':near,'S2_eigen_residual':float(la.norm(w-near*(near+1)*v))}
 finally:
  m.flips,m.js,m.diag=old
  if 'bonds' in locals():m.bonds=bonds
