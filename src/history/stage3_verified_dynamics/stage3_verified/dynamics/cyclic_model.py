"""Additional compatible cyclic little-group reduction; no spin sector assumed."""
import numpy as np,time
from ed import Model,ptr
from sectors import point_permutations
class CyclicModel(Model):
 def __init__(self,L=3,up=14,jp=.2,cell_n=(0,0),angular=0,threads=2):
  started=time.perf_counter();super().__init__(L,up,False,jp,threads);self.cell_n=tuple(cell_n);self.angular=angular
  k=np.array(cell_n)/L;self.baseq=np.array([(k[0]+k[1])/3,(2*k[0]-k[1])/3]);oldperm=self.perms.copy();oldgs=self.gs.copy();pg=point_permutations(L,jp)
  if self.cell_n==(0,0):self.order=6;gen=pg[1]
  elif L==3 and self.cell_n in [(1,2),(2,1)]:self.order=3;gen=pg[2]
  else:
   self.order=2;gen=None;chars=np.exp(-2j*np.pi*(oldgs@self.baseq));lookup={tuple(p):i for i,p in enumerate(oldperm)}
   for p in pg[6:]:
    inv=np.argsort(p)
    if all(abs(chars[lookup[tuple(p[t[inv]])]]-chars[j])<1e-9 for j,t in enumerate(oldperm)):gen=p;break
   if gen is None:raise ValueError('No compatible reflection')
  assert 0<=angular<self.order
  powers=[];p=np.arange(self.n)
  for r in range(self.order):powers.append(p.copy());p=gen[p]
  assert np.array_equal(p,np.arange(self.n))
  self.perms=np.array([t[p] for p in powers for t in oldperm]);self.gs=np.tile(oldgs,(self.order,1));self.rotation_powers=np.repeat(np.arange(self.order),len(oldperm))
  split=self.split;nl=1<<split;nh=1<<(self.n-split)
  def lut(start,size):
   out=np.zeros((len(self.gs),size),np.uint32)
   for g,p in enumerate(self.perms):
    for mask in range(1,size):
     bit=mask&-mask;out[g,mask]=out[g,mask^bit]|(1<<int(p[start+bit.bit_length()-1]))
   return out
  self.lo=lut(0,nl);self.hi=lut(split,nh);self.rid=np.full(self.dim,-1,np.int32);self.shift=np.zeros(self.dim,np.uint8)
  reps=np.empty(self.dim,np.uint32);lens=np.empty(self.dim,np.uint8)
  nr=self.lib.orbits(self.n,up,self.dim,len(self.gs),ptr(self.lo),ptr(self.hi),split,ptr(self.rid),ptr(self.shift),ptr(reps),ptr(lens))
  self.reps=reps[:nr].copy();self.lens=lens[:nr].copy();self.setup_seconds=time.perf_counter()-started;self.sector(self.baseq)
 def sector(self,q=None):
  q=self.baseq if q is None else q;assert np.max(abs(np.asarray(q)-self.baseq))<1e-10
  self.q=np.array(q);self.phase=np.ascontiguousarray(np.exp(-2j*np.pi*(self.gs@self.q+self.angular*self.rotation_powers/self.order)))
  valid=np.ones(len(self.reps),bool)
  for g,ph in enumerate(self.phase):
   if abs(ph-1)>1e-9:
    t=self.lo[g,self.reps&((1<<self.split)-1)]|self.hi[g,self.reps>>self.split];valid&=t!=self.reps
  self.rawbasis=np.full(len(valid),-1,np.int32);self.rawbasis[valid]=np.arange(valid.sum(),dtype=np.int32)
  self.masks=self.reps[valid];self.lengths=self.lens[valid];self.nb=len(self.masks);self.diag=np.zeros(self.nb)
  for a,b,j in self.bonds:self.diag+=j*(.25-.5*(((self.masks>>a)^(self.masks>>b))&1))
  self.calls=0;self.mv_seconds=0.
  from scipy.sparse.linalg import LinearOperator
  return LinearOperator((self.nb,self.nb),matvec=self.mv,dtype=np.complex128)
