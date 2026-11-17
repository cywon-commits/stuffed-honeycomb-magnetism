import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import ctypes as ct,subprocess,math,itertools,json,numpy as np,scipy.linalg as la
from geometry import ROOT,make
from pilot36 import ptr,P
if not (ROOT/'kernel36pilot.so').exists() or (ROOT/'kernel36pilot.cpp').stat().st_mtime>(ROOT/'kernel36pilot.so').stat().st_mtime:subprocess.run(['g++','-O3','-march=native','-shared','-fPIC','-fopenmp',str(ROOT/'kernel36pilot.cpp'),'-o',str(ROOT/'kernel36pilot.so')],check=True)
G=make([[2,2],[-2,4]]);lib=ct.CDLL(str(ROOT/'kernel36pilot.so'));lib.init36();lib.orbits36.argtypes=[ct.c_int,ct.c_int,ct.c_uint64,ct.c_int]+[P]*5;lib.orbits36.restype=ct.c_int;lib.mv36.argtypes=[ct.c_int,ct.c_int]+[P]*10+[ct.c_int]
n=36;up=2;dim=math.comb(n,up);rid=np.full(dim,-1,np.int32);shift=np.zeros(dim,np.uint8);reps=np.empty(dim,np.uint64);lens=np.empty(dim,np.uint8);perms=np.ascontiguousarray(G['perms'],np.int32);nr=lib.orbits36(n,up,dim,12,*[ptr(x) for x in [perms,rid,shift,reps,lens]]);masks=reps[:nr].copy();lengths=lens[:nr].copy();flips=np.array([(1<<a)|(1<<b) for a,b,t in G['bonds']],np.uint64);js=np.array([1. if t==0 else .3 for a,b,t in G['bonds']]);diag=sum(j*(.25-.5*(((masks>>a)^(masks>>b))&1)) for (a,b,t),j in zip(G['bonds'],js));phase=np.ones(12,complex);rng=np.random.default_rng(36);v=rng.normal(size=nr)+1j*rng.normal(size=nr);v/=la.norm(v);hv=np.empty_like(v);lib.mv36(nr,len(flips),*[ptr(a) for a in [masks,lengths,flips,js,diag,rid,shift,phase,v,hv]],2)
states=np.sort(np.array([sum(1<<a for a in c) for c in itertools.combinations(range(36),2)],np.uint64));lookup={int(s):i for i,s in enumerate(states)};psi=v[rid]/np.sqrt(lengths[rid].astype(float));expected=np.zeros(dim,complex)
for i,s0 in enumerate(states):
 s=int(s0)
 for a,b,t in G['bonds']:
  j=1 if t==0 else .3;diff=((s>>a)^(s>>b))&1;expected[i]+=j*(.25-.5*diff)*psi[i]
  if diff:expected[lookup[s^((1<<a)|(1<<b))]]+=.5*j*psi[i]
expanded=hv[rid]/np.sqrt(lengths[rid].astype(float));err=float(la.norm(expected-expanded));assert err<1e-11
(ROOT/'pilot36_independent_validation.json').write_text(json.dumps({'status':'PASS','N':36,'number_of_down_spins':2,'full_basis_dim':dim,'translation_Gamma_dim':nr,'independent_Python_full_vs_Cpp_reduced_action_error':err,'scope':'36-site high-Sz kernel validation, not a low-spin ground-state calculation'},indent=2));print('PASS',err)
