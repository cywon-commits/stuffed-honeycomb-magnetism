import pathlib,ctypes as ct,subprocess,numpy as np,time
from ed import ptr,P
from local_geometry import faces,matrices
ROOT=pathlib.Path(__file__).resolve().parent
_lib=None
def load():
 global _lib
 if _lib is None:
  src=ROOT/'triangle_ops.cpp';so=ROOT/'triangle_ops.so'
  if not so.exists() or so.stat().st_mtime<src.stat().st_mtime:subprocess.run(['g++','-O3','-march=native','-shared','-fPIC','-fopenmp',str(src),'-o',str(so)],check=True)
  _lib=ct.CDLL(str(so));_lib.triangle_seeds.argtypes=[P]*17
 return _lib

def seeds(source,target,v,k,rdm=False):
 f=faces(source.L)[0];ops=matrices();fast=bool(np.all(source.lengths==source.L**2) and np.all(target.lengths==source.L**2) and np.array_equal(source.masks,target.masks));cfg=np.array([target.nb,source.L**2,source.split,source.rh.shape[1],source.threads,int(fast),int(rdm)],np.int32);weights=np.ascontiguousarray([np.exp(-2j*np.pi*(k[0]*x+k[1]*y))/source.L for x in range(source.L) for y in range(source.L)],complex);out=np.zeros((18,target.nb),complex);rho=np.zeros((6,8,8),complex)
 if rdm:assert np.array_equal(source.masks,target.masks)
 load().triangle_seeds(*[ptr(a) for a in [cfg,target.masks,target.lengths,source.rid,source.shift,source.rawbasis,source.lengths,source.pop,source.rl,source.rh,source.phase,np.ascontiguousarray(v),f,weights,ops,out,rho]])
 return out,rho
