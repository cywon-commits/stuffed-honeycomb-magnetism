"""Translation-reduced spin-1/2 triangular Heisenberg ED. See report.md."""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import numpy as np, scipy.linalg as la
from scipy.sparse.linalg import LinearOperator
from math import comb
import ctypes as ct, subprocess, pathlib, time, json, resource, argparse
ROOT=pathlib.Path(__file__).resolve().parent
P=ct.c_void_p; I=ct.c_int

def libload():
    target=ROOT/'kernel.so'
    if not target.exists() or target.stat().st_mtime<(ROOT/'kernel.cpp').stat().st_mtime:
        subprocess.run(['g++','-O3','-march=native','-fopenmp','-shared','-fPIC',str(ROOT/'kernel.cpp'),'-o',str(target)],check=True)
    lib=ct.CDLL(str(target));lib.init()
    lib.orbits.argtypes=[I,I,I,I,P,P,I,P,P,P,P];lib.orbits.restype=I
    lib.mv.argtypes=[I,I,P,P,P,P,P,P,P,P,P,P,P,I,I,P,P,P,I]
    return lib

def ptr(a): return a.ctypes.data_as(P)
def geometry(L,jp=1.):
    def sid(x,y,s):return 3*((x%L)*L+y%L)+s
    coords=[];b=[]
    for x in range(L):
      for y in range(L):
        coords.extend([(x+2*y,x-y),(x+2*y+1,x-y),(x+2*y,x-y+1)])
        a=sid(x,y,0)
        b.extend((a,sid(xx,yy,1),1.) for xx,yy in [(x,y),(x-1,y),(x,y-1)])
        b.extend((sid(x,y,2),sid(xx,yy,s),jp) for xx,yy,s in [(x,y,0),(x,y,1),(x+1,y,0),(x+1,y-1,1),(x+1,y-1,0),(x,y-1,1)])
    return np.array(coords),b

def translations(L,coords,primitive=True):
    def key(u,v):return ((u+2*v)%(3*L),(u-v)%(3*L))
    lookup={key(*c):i for i,c in enumerate(coords)}
    gs=coords.copy() if primitive else np.array([(x+2*y,x-y) for x in range(L) for y in range(L)])
    perms=np.array([[lookup[key(*(c+g))] for c in coords] for g in gs])
    return gs,perms

class Model:
 def __init__(self,L=3,up=14,primitive=True,jp=1.,threads=2):
    t=time.perf_counter();self.L=L;self.up=up;self.n=3*L*L;self.dim=comb(self.n,up);self.threads=threads
    assert self.n<=30
    if primitive and jp!=1:raise ValueError('Primitive translations require Jp=J=1')
    self.lib=libload();self.coords,self.bonds=geometry(L,jp);self.gs,self.perms=translations(L,self.coords,primitive)
    self.split=self.n//2;nl=1<<self.split;nh=1<<(self.n-self.split)
    def lut(start,size):
      out=np.zeros((len(self.gs),size),np.uint32)
      for g,perm in enumerate(self.perms):
       for m in range(1,size):
        bit=m&-m;out[g,m]=out[g,m^bit]|(1<<int(perm[start+bit.bit_length()-1]))
      return out
    self.lo=lut(0,nl);self.hi=lut(self.split,nh)
    self.rid=np.full(self.dim,-1,np.int32);self.shift=np.zeros(self.dim,np.uint8)
    reps=np.empty(self.dim,np.uint32);lens=np.empty(self.dim,np.uint8)
    nr=self.lib.orbits(self.n,up,self.dim,len(self.gs),ptr(self.lo),ptr(self.hi),self.split,ptr(self.rid),ptr(self.shift),ptr(reps),ptr(lens))
    self.reps=reps[:nr].copy();self.lens=lens[:nr].copy()
    self.pop=np.array([m.bit_count() for m in range(nl)],np.uint8)
    self.rl=np.zeros(nl,np.int32);self.rh=np.zeros((self.split+1,nh),np.int32)
    for m in range(nl):self.rl[m]=sum(comb(p,k) for k,p in enumerate([p for p in range(self.split) if m>>p&1],1))
    for m in range(nh):
      ps=[p+self.split for p in range(self.n-self.split) if m>>p&1]
      for kl in range(self.split+1):self.rh[kl,m]=sum(comb(p,kl+k) if kl+k<=p else 0 for k,p in enumerate(ps,1))
    self.flips=np.array([(1<<a)|(1<<b) for a,b,j in self.bonds],np.uint32)
    self.js=np.array([j for a,b,j in self.bonds],float)
    self.setup_seconds=time.perf_counter()-t
 def sector(self,q):
    self.q=np.array(q);self.phase=np.ascontiguousarray(np.exp(-2j*np.pi*(self.gs@self.q)))
    # Character must be single valued on the periodic torus.
    for a in [(self.L,self.L),(2*self.L,-self.L)]:
     if abs(np.exp(2j*np.pi*np.dot(q,a))-1)>1e-9:raise ValueError('Momentum incompatible with periodic boundary conditions')
    valid=np.ones(len(self.reps),bool)
    for g,ph in enumerate(self.phase):
     if abs(ph-1)>1e-9:
      translated=self.lo[g,self.reps&((1<<self.split)-1)]|self.hi[g,self.reps>>self.split]
      valid &= translated!=self.reps
    self.rawbasis=np.full(len(valid),-1,np.int32);self.rawbasis[valid]=np.arange(valid.sum(),dtype=np.int32)
    self.masks=self.reps[valid];self.lengths=self.lens[valid];self.nb=len(self.masks)
    self.diag=np.zeros(self.nb)
    for a,b,j in self.bonds:self.diag+=j*(.25-.5*(((self.masks>>a)^(self.masks>>b))&1))
    self.calls=0;self.mv_seconds=0.
    return LinearOperator((self.nb,self.nb),matvec=self.mv,dtype=np.complex128)
 def mv(self,x):
    x=np.ascontiguousarray(x,dtype=np.complex128);y=np.empty_like(x);t=time.perf_counter()
    self.lib.mv(self.nb,len(self.bonds),*[ptr(a) for a in [self.masks,self.lengths,self.diag,self.flips,self.js,self.rid,self.shift,self.rawbasis,self.pop,self.rl,self.rh]],self.split,self.rh.shape[1],ptr(self.phase),ptr(x),ptr(y),self.threads)
    self.calls+=1;self.mv_seconds+=time.perf_counter()-t
    return y

# Explicitly reorthogonalized Hermitian Lanczos with exact Ritz residual verification.
# Bounded basis; on reaching maxiter, saves a partial Ritz vector and labels it unconverged.
def lanczos(model,maxiter=180,tol=1e-8,check_every=10,seed=2026,checkpoint=None,copy_reorth=False):
    rng=np.random.default_rng(seed);v=rng.standard_normal(model.nb)+1j*rng.standard_normal(model.nb);v/=la.norm(v)
    # Row-major basis permits BLAS block reorthogonalization. Allocated storage is bounded.
    Q=np.empty((maxiter,model.nb),complex);alpha=[];beta=[];history=[];t=time.perf_counter();converged=False
    for it in range(maxiter):
      Q[it]=v;w=model.mv(v);a=np.vdot(v,w)
      if abs(a.imag)>1e-10:raise RuntimeError('Non-Hermitian Lanczos diagonal')
      alpha.append(a.real);w-=a.real*v
      if it:w-=beta[-1]*Q[it-1]
      for repeat in range(2):
        if copy_reorth:
          coeff=Q[:it+1].conj()@w;w-=coeff@Q[:it+1]
        else:
          coeff=la.blas.zgemv(1.,Q[:it+1].T,w,trans=2)
          w-=la.blas.zgemv(1.,Q[:it+1].T,coeff)
      b=la.norm(w)
      if (it+1)%check_every==0 or it+1==maxiter or b<1e-13:
        d,z=la.eigh_tridiagonal(np.array(alpha),np.array(beta),select='i',select_range=(0,0)) if it else (np.array(alpha),np.ones((1,1)))
        est=b*abs(z[-1,0]);record={'iteration':it+1,'ritz':float(d[0]),'recurrence_residual':float(est),'elapsed_seconds':time.perf_counter()-t}
        if est<tol or it+1==maxiter or b<1e-13:
          psi=z[:,0]@Q[:it+1];psi/=la.norm(psi);hpsi=model.mv(psi);energy=np.vdot(psi,hpsi).real;res=la.norm(hpsi-energy*psi)
          record.update(energy=float(energy),residual=float(res));converged=res<tol
        history.append(record);print(json.dumps(record),flush=True)
        if checkpoint:pathlib.Path(str(checkpoint)+'.history.json').write_text(json.dumps(history,indent=2))
        if converged or b<1e-13:break
      beta.append(b);v=w/b
    if checkpoint:np.savez_compressed(str(checkpoint)+'.npz',vector=psi,energy=energy,q=model.q)
    return {'energy':float(energy),'residual':float(res),'converged':bool(converged),'iterations':it+1,'matvec_calls':model.calls,'matvec_seconds':model.mv_seconds,'solve_seconds':time.perf_counter()-t,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'history':history}

def main():
 p=argparse.ArgumentParser();p.add_argument('--copy-reorth',action='store_true');p.add_argument('--maxiter',type=int,default=180);p.add_argument('--tol',type=float,default=1e-8);p.add_argument('--q',nargs=2,type=float,default=[0,0]);p.add_argument('--name',default='run');p.add_argument('--threads',type=int,default=2);a=p.parse_args()
 m=Model(threads=a.threads);m.sector(a.q);print(json.dumps({'setup_seconds':m.setup_seconds,'dim':m.dim,'sector_dim':m.nb,'q':a.q}),flush=True)
 rng=np.random.default_rng(0);x=rng.normal(size=m.nb)+1j*rng.normal(size=m.nb);x/=la.norm(x);u=rng.normal(size=m.nb)+1j*rng.normal(size=m.nb);u/=la.norm(u)
 t=time.perf_counter();hx=m.mv(x);benchmark=time.perf_counter()-t;hu=m.mv(u);herm=abs(np.vdot(u,hx)-np.vdot(hu,x))
 m.calls=0;m.mv_seconds=0
 result=lanczos(m,a.maxiter,a.tol,checkpoint=ROOT/a.name,copy_reorth=a.copy_reorth)
 result.update(maxiter=a.maxiter,tolerance=a.tol,seed=2026,threads=a.threads,copy_reorth=a.copy_reorth,q=a.q,full_dim=m.dim,sector_dim=m.nb,setup_seconds=m.setup_seconds,benchmark_matvec_seconds=benchmark,hermiticity_absolute=float(herm))
 (ROOT/(a.name+'.json')).write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=='__main__':main()
