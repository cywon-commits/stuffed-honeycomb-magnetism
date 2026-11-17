import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,scipy.linalg as la,itertools,pathlib,ctypes as ct,subprocess,time
from scipy.sparse.linalg import LinearOperator,eigsh
from geometry import ROOT,make
LIB=None
def lib():
 global LIB
 if LIB is None:
  if not (ROOT/'kernel18.so').exists() or (ROOT/'kernel18.cpp').stat().st_mtime>(ROOT/'kernel18.so').stat().st_mtime:subprocess.run(['g++','-O3','-march=native','-shared','-fPIC','-fopenmp',str(ROOT/'kernel18.cpp'),'-o',str(ROOT/'kernel18.so')],check=True)
  LIB=ct.CDLL(str(ROOT/'kernel18.so'));LIB.action.argtypes=[ct.c_int,ct.c_int]+[ct.c_void_p]*11+[ct.c_int]
 return LIB
def ptr(x):return x.ctypes.data_as(ct.c_void_p)
class Basis:
 def __init__(self,g,up=None,full=False):
  self.g=g;self.n=g['N'];self.up=self.n//2 if up is None else up;assert self.n<=22;self.states=np.sort(np.array([sum(1<<a for a in c) for c in itertools.combinations(range(self.n),self.up)],np.uint32));self.lookup=np.full(1<<self.n,-1,np.int32);self.lookup[self.states]=np.arange(len(self.states));self.full=full
  perms=np.arange(self.n)[None,:] if full else g['perms'];self.trans=np.zeros((1,2),int) if full else g['cells'];self.perms=perms;self.split=self.n//2
  def lut(offset,size):
   t=np.zeros((len(perms),size),np.uint32)
   for k,p in enumerate(perms):
    for m in range(1,size):b=m&-m;t[k,m]=t[k,m^b]|(1<<int(p[offset+b.bit_length()-1]))
   return t
  self.lo=lut(0,1<<self.split);self.hi=lut(self.split,1<<(self.n-self.split));ts=self.lo[:,self.states&((1<<self.split)-1)]|self.hi[:,self.states>>self.split];mins=ts.min(axis=0);self.reps,inv,count=np.unique(mins,return_inverse=True,return_counts=True);self.lengths_all=count.astype(np.uint8);self.repidx=np.full(1<<self.n,-1,np.int32);self.repidx[self.states]=inv;self.shift=np.zeros(1<<self.n,np.uint8)
  for k in reversed(range(len(perms))):
   tt=self.lo[k,self.reps&((1<<self.split)-1)]|self.hi[k,self.reps>>self.split];self.shift[tt]=k
  self.rawbasis=np.arange(len(self.reps),dtype=np.int32);self.q=np.zeros(2);self.sector([0,0])
 def sector(self,q):
  self.q=np.array(q,float);self.phase=np.ascontiguousarray(np.exp(-2j*np.pi*(self.trans@q)));valid=np.ones(len(self.reps),bool)
  for k,ph in enumerate(self.phase):
   if abs(ph-1)>1e-9:valid&=(self.lo[k,self.reps&((1<<self.split)-1)]|self.hi[k,self.reps>>self.split])!=self.reps
  self.rawbasis=np.full(len(self.reps),-1,np.int32);self.rawbasis[valid]=np.arange(sum(valid),dtype=np.int32);self.masks=self.reps[valid].copy();self.lens=self.lengths_all[valid].copy();self.nb=len(self.masks);self.calls=0;self.mv_seconds=0.;self.set_r(1.)
 def set_r(self,r):
  self.r=r;self.edges=self.g['bonds'];self.flips=np.array([(1<<a)|(1<<b) for a,b,t in self.edges],np.uint32);self.js=np.array([1. if t==0 else r for a,b,t in self.edges]);self.diag=sum(j*(.25-.5*(((self.masks>>a)^(self.masks>>b))&1)) for (a,b,t),j in zip(self.edges,self.js))
 def action(self,v,flips=None,js=None,diag=None):
  v=np.ascontiguousarray(v,complex);out=np.empty_like(v);args=[self.masks,self.lens,self.flips if flips is None else flips,self.js if js is None else js,self.diag if diag is None else diag,self.repidx,self.shift,self.rawbasis,self.phase,v,out]
  # Eleven pointers after integer sizes.
  lib().action(self.nb,len(self.js if js is None else js),*[ptr(a) for a in args],2);return out
 def mv(self,v):
  t=time.perf_counter();out=self.action(v);self.calls+=1;self.mv_seconds+=time.perf_counter()-t;return out
 def spin_action(self,v):
  pairs=list(itertools.combinations(range(self.n),2));fl=np.array([(1<<a)|(1<<b) for a,b in pairs],np.uint32);js=np.full(len(pairs),2.);diag=np.full(self.nb,.75*self.n)
  for a,b in pairs:diag+=2*(.25-.5*(((self.masks>>a)^(self.masks>>b))&1))
  return self.action(v,fl,js,diag)
 def expand(self,v):
  ids=self.rawbasis[self.repidx[self.states]];out=np.zeros(len(self.states),complex);ok=ids>=0;out[ok]=v[ids[ok]]*self.phase[self.shift[self.states[ok]]]/np.sqrt(self.lens[ids[ok]].astype(float));return out
 def spin(self,v):
  s=self.spin_action(v);S2=np.vdot(v,s).real;S=round((np.sqrt(max(0,1+4*S2))-1)/2);return S,float(S2),float(la.norm(s-S*(S+1)*v))
