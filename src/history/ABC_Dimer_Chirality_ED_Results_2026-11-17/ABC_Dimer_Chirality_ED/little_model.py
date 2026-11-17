"""Primitive translations plus an abelian little-group character; optional C3v A1/A2."""
import time,numpy as np
from ed import Model,ptr
from sectors import point_permutations
from observables import full_samples
class LittleModel(Model):
 def __init__(self,q,angular=0,reflection=None,threads=2,primitive=True,jp=1.):
  super().__init__(primitive=primitive,jp=jp,threads=threads);self.baseq=np.asarray(q);pg=point_permutations();oldperm=self.perms.copy();oldgs=self.gs.copy();lookup={tuple(p):i for i,p in enumerate(oldperm)};chars=np.exp(-2j*np.pi*(oldgs@q))
  def fixes(p):
   inv=np.argsort(p)
   return all(abs(chars[lookup[tuple(inv[t[p]])]]-chars[j])<1e-9 for j,t in enumerate(oldperm))
  rotations=[p for p in pg[:6] if fixes(p)];self.rot_order=len(rotations);assert self.rot_order in [1,3,6];assert 0<=angular<self.rot_order
  if reflection is not None:assert abs(np.exp(-2j*np.pi*angular/self.rot_order).imag)<1e-10
  group=rotations.copy();char=[np.exp(-2j*np.pi*angular*r/self.rot_order) for r in range(self.rot_order)]
  self.mirror_index=None
  if reflection is not None:
   self.mirror_index=next(i for i in range(6,12) if fixes(pg[i]));mirror=pg[self.mirror_index];group+=[mirror[p] for p in rotations];char+=list(reflection*np.asarray(char));
  self.little_group=np.array(group);self.little_char=np.array(char);self.angular=angular;self.reflection=reflection;self.perms=np.array([t[p] for p in group for t in oldperm]);self.gs=np.tile(oldgs,(len(group),1));assert len(self.perms)<256
  self.phase=np.ascontiguousarray((self.little_char[:,None]*np.exp(-2j*np.pi*(oldgs@q))[None,:]).ravel());self.q=np.asarray(q)
  def lut(start,size):
   out=np.zeros((len(self.gs),size),np.uint32)
   for g,p in enumerate(self.perms):
    for mask in range(1,size):
     bit=mask&-mask;out[g,mask]=out[g,mask^bit]|(1<<int(p[start+bit.bit_length()-1]))
   return out
  self.lo=lut(0,1<<self.split);self.hi=lut(self.split,1<<(self.n-self.split));self.rid=np.full(self.dim,-1,np.int32);self.shift=np.zeros(self.dim,np.uint8);reps=np.empty(self.dim,np.uint32);lens=np.empty(self.dim,np.uint8);nr=self.lib.orbits(self.n,self.up,self.dim,len(self.gs),ptr(self.lo),ptr(self.hi),self.split,ptr(self.rid),ptr(self.shift),ptr(reps),ptr(lens));self.reps=reps[:nr].copy();self.lens=lens[:nr].copy();valid=np.ones(nr,bool)
  for g,ph in enumerate(self.phase):
   if abs(ph-1)>1e-9:valid&=(self.lo[g,self.reps&((1<<self.split)-1)]|self.hi[g,self.reps>>self.split])!=self.reps
  self.rawbasis=np.full(nr,-1,np.int32);self.rawbasis[valid]=np.arange(valid.sum(),dtype=np.int32);self.masks=self.reps[valid];self.lengths=self.lens[valid];self.nb=len(self.masks);self.diag=np.zeros(self.nb)
  for a,b,j in self.bonds:self.diag+=j*(.25-.5*(((self.masks>>a)^(self.masks>>b))&1))
  self.calls=0;self.mv_seconds=0.
 def project(self,source,B):
  out=np.zeros((self.nb,B.shape[1]),complex)
  for p,ch in zip(self.little_group,self.little_char):
   tables=[]
   for start,size in [(0,1<<self.split),(self.split,1<<(self.n-self.split))]:
    tab=np.zeros(size,np.uint32)
    for mask in range(1,size):
     bit=mask&-mask;tab[mask]=tab[mask^bit]|(1<<int(p[start+bit.bit_length()-1]))
    tables.append(tab)
   pre=tables[0][self.masks&((1<<self.split)-1)]|tables[1][self.masks>>self.split]
   for a,b in enumerate(B.T):out[:,a]+=ch.conjugate()*full_samples(source,b,pre)*np.sqrt(self.lengths.astype(float))/len(self.little_group)
  return out
