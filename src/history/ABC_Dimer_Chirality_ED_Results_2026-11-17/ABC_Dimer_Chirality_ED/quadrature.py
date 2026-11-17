"""Fully reorthogonalized block Krylov quadrature, explicitly not verified poles."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,time,resource,json,numpy as np,scipy.linalg as la
from experiment import ROOT,Model,pack,dump,gram

def run(name,dimension=96):
 if name=='j09_q20':name='j09_q'
 start=time.perf_counter();cpu_start=time.process_time();z=np.load(ROOT/(name+'_seeds.npz'));j=float(z['j']);requested_dimension=dimension;dimension=min(dimension,60) if j!=1 else dimension;m=Model(primitive=j==1,jp=j);m.sector(z['q']);B=z['B'];E0=float(z['E0']);M=gram(B);e,U=la.eigh(M);keep=e>max(e.max()*1e-12,1e-13);r=int(sum(keep));d=(dimension//r)*r;Q=np.empty((d,m.nb),complex);Q[:r]=(B@U[:,keep]/np.sqrt(e[keep])).T;R=gram(Q[:r].T,B);del M,e,U;T=np.zeros((d,d),complex);worstorth=0.
 for a in range(0,d,r):
  b=a+r;W=np.asfortranarray(np.column_stack([m.mv(Q[i])-E0*Q[i] for i in range(a,b)]));C=la.blas.zgemm(1.,Q[:b].T,W,trans_a=2);T[:b,a:b]=C;T[a:b,:b]=C.conj().T
  for _ in range(2):
   C=la.blas.zgemm(1.,Q[:b].T,W,trans_a=2);W-=la.blas.zgemm(1.,Q[:b].T,C)
  q,rr=la.qr(W,mode='economic');worstorth=max(worstorth,float(la.norm(la.blas.zgemm(1.,Q[:b].T,q,trans_a=2))))
  if b<d:Q[b:b+r]=q.T
  print(json.dumps({'name':name,'basis_dimension':b,'block_rank':r,'seconds':time.perf_counter()-start,'orth_error':worstorth}),flush=True)
 T=(T+T.conj().T)/2;omega,V=la.eigh(T);A=V[:r].conj().T@R;res=la.norm(rr@V[-r:],axis=0);moments=[A.conj().T@(omega[:,None]**p*A) for p in range(3)];err=[float(la.norm(moments[p]-z[f'M{p}'])) for p in range(3)];assert max(err)<1e-6
 np.savez_compressed(ROOT/(name+'_quadrature.npz'),omega=omega,amplitudes=A,recurrence_residual=res)
 result={'name':name,'status':'KRYLOV_QUADRATURE_NOT_INDEPENDENTLY_CONVERGED_EIGENVALUES','dimension':d,'requested_dimension':requested_dimension,'block_rank':r,'moment_errors':err,'orthogonality_error':worstorth,'H_calls':m.calls,'CPU_seconds':time.process_time()-cpu_start,'wall_seconds':time.perf_counter()-start,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'nodes':[{'omega':float(x),'recurrence_residual':float(res[i]),'residue':pack(np.outer(A[i].conj(),A[i]))} for i,x in enumerate(omega)]};dump(ROOT/(name+'_quadrature.json'),result);print(json.dumps({k:v for k,v in result.items() if k!='nodes'}),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('name');p.add_argument('--dimension',type=int,default=96);a=p.parse_args();run(a.name,a.dimension)
