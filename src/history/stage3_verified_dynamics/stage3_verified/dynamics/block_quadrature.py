"""Finite block-Krylov spectral quadrature. Nodes are NOT claimed exact eigenvalues."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,time,resource,numpy as np,scipy.linalg as la
from ed import Model
from crossing_transitions import pack
ROOT=pathlib.Path(__file__).resolve().parent
def run(blocks=24):
 start=time.perf_counter();z=np.load(ROOT/'weak_seeds.npz');B=z['B'];HB=z['HB'];E=float(z['E0']);m=Model(primitive=False,jp=.2375);m.sector(z['q']);del z
 d=3*blocks;Q=np.empty((d,m.nb),complex);q,r0=la.qr(B,mode='economic');Q[:3]=q.T;del q;T=np.zeros((d,d),complex);orth=0.
 for step in range(blocks):
  a=3*step;b=a+3;W=np.column_stack([m.mv(Q[i])-E*Q[i] for i in range(a,b)]);C=(Q[:b]@W.conj()).conj();T[:b,a:b]=C;T[a:b,:b]=C.conj().T
  for repeat in range(2):
   c=(Q[:b]@W.conj()).conj();W-=Q[:b].T@c
  q,r=la.qr(W,mode='economic');orth=max(orth,float(la.norm((Q[:b]@q.conj()).conj())))
  if b<d:Q[b:b+3]=q.T
  print(json.dumps({'block':step+1,'dimension':b,'seconds':time.perf_counter()-start,'orth_error':orth}),flush=True)
 T=(T+T.conj().T)/2;energies,U=la.eigh(T);A=U[:3].conj().T@r0;res=la.norm(r@U[-3:],axis=0);mom=[A.conj().T@(energies[:,None]**p*A) for p in range(3)];exact=[B.conj().T@B,B.conj().T@HB,HB.conj().T@HB];errors=[float(la.norm(mom[p]-exact[p])) for p in range(3)];assert max(errors)<1e-7
 nodes=[{'omega':float(e),'estimated_eigen_residual':float(res[i]),'residue_matrix':pack(np.outer(A[i].conj(),A[i])),'status':'KRYLOV_QUADRATURE_NODE_NOT_VERIFIED_EIGENVALUE'} for i,e in enumerate(energies)]
 result={'kind':'Finite block-Krylov quadrature, separate from verified physical poles','basis_dimension':d,'nodes':nodes,'moments':{f'M{p}':pack(mom[p]) for p in range(3)},'moment_errors':errors,'orthogonality_error':orth,'H_calls':m.calls,'wall_seconds':time.perf_counter()-start,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'broadening':None}
 (ROOT/'weak_quadrature.json').write_text(json.dumps(result,indent=2));print(json.dumps({k:v for k,v in result.items() if k not in ['nodes','moments']}),flush=True)
if __name__=='__main__':run()
