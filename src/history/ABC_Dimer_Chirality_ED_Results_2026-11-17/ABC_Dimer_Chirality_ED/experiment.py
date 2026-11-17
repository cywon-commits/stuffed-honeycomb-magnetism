"""Matrix-free scalar spectroscopy; unbroadened residues and independently checked poles."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='2'
import pathlib,json,time,resource,argparse,gc,csv
import numpy as np,scipy.linalg as la
from scipy.sparse.linalg import eigsh,LinearOperator,ArpackNoConvergence
from ed import Model
from cyclic_model import CyclicModel
from native import seeds
from observables import convert_basis,cell_q,total_spin_check
from local_geometry import faces,w
ROOT=pathlib.Path(__file__).resolve().parent

def pack(x):return {'real':np.asarray(x).real.tolist(),'imag':np.asarray(x).imag.tolist()}
def dump(p,x):pathlib.Path(p).write_text(json.dumps(x,indent=2))
def gram(x,y=None):return la.blas.zgemm(1.,np.asfortranarray(x),np.asfortranarray(x if y is None else y),trans_a=2)
def geometry_weights(Q):
 _,_,c,verts=faces();g=np.ones(6,complex)
 for mu,pts in enumerate(verts[:6]):
  order=sorted(range(3),key=lambda i:(round(pts[i,1],10),round(pts[i,0]+pts[i,1]/2,10)));g[mu]=w**(-order[0])
 W=np.zeros((18,6),complex)
 for kind,geo in enumerate([g,g.conj(),np.ones(6)]):
  for mu in range(6):W[6*kind+mu,2*kind+mu%2]=geo[mu]*np.exp(-2j*np.pi*(c[mu]@Q))/np.sqrt(3)
 return W

def monitored(m,label):
 start=time.perf_counter();hist=[];base=m.calls
 def mv(v):
  y=m.mv(v)
  if (m.calls-base)%25==0:
   r={'label':label,'H_calls':m.calls-base,'seconds':time.perf_counter()-start};hist.append(r);print(json.dumps(r),flush=True);dump(ROOT/(label+'_progress.json'),hist)
  return y
 return LinearOperator((m.nb,m.nb),matvec=mv,dtype=complex)

def prepare_source():
 dest=ROOT/'source_j1.npz'
 if dest.exists():return
 old=np.load(ROOT/'inputs/K_primitive.npz');m=Model();m.sector(old['q']);v=old['vector'];E=float(old['energy']);original={'energy':E,'H_residual':float(la.norm(m.mv(v)-E*v)),'spin':total_spin_check(m,v)};start=time.perf_counter();vals,V=eigsh(monitored(m,'refine_source'),k=1,which='SA',v0=v,ncv=24,tol=2e-14,maxiter=150);v=V[:,0];E=float(vals[0]);checks={'energy':E,'H_residual':float(la.norm(m.mv(v)-E*v)),'spin':total_spin_check(m,v),'wall_seconds':time.perf_counter()-start,'H_calls':m.calls,'original':original};assert checks['H_residual']<1e-10 and checks['spin']['S2_eigen_residual']<1e-8;np.savez_compressed(dest,vector=v,energy=E,q=m.q);dump(ROOT/'source_j1.json',checks);print(json.dumps(checks),flush=True)

def construct(name,j,Q,partner=False,cell=False):
 start=time.perf_counter();cpu_start=time.process_time();Q=np.asarray(Q,float)
 if j==1:
  z=np.load(ROOT/'source_j1.npz');s=Model();q=z['q'];q=(-q)%1 if partner else q;s.sector(q);v=z['vector'].conj() if partner else z['vector'];v=v/la.norm(v);E=float(z['energy']);t=Model();t.sector((q+Q)%1);k=np.mod([Q[0]+Q[1],2*Q[0]-Q[1]],1);out,_=seeds(s,t,v,k);B=np.asfortranarray(out.T@geometry_weights(Q));del out
  names=['D_up','D_down','Ddag_up','Ddag_down','chi_up','chi_down']
  # Ground orbital K/Kprime pair; only one lies in any target primitive sector.
  G=[]
  for qp,vg in [(z['q'],z['vector']),((-z['q'])%1,z['vector'].conj())]:
   if np.max(abs(np.exp(2j*np.pi*(t.gs@(t.q-qp)))-1))<1e-8:G.append(vg/la.norm(vg))
 else:
  z=np.load(ROOT/f'inputs/p{j:.9f}_n00_r3_e0.npz');small=CyclicModel(jp=j,angular=3);s=Model(primitive=False,jp=j);s.sector([0,0]);v=convert_basis(small,s,z['vector']);del small;v/=la.norm(v);E=float(z['energy']);t=Model(primitive=False,jp=j);t.sector(cell_q(Q));B0,_=seeds(s,t,v,Q);B=np.asfortranarray(B0.T);del B0;k=Q;names=[f'{a}_mu{mu}' for a in ['D','Ddag','chi'] for mu in range(6)];G=[v] if la.norm(Q)<1e-10 else []
 sourcecheck={'H_residual':float(la.norm(s.mv(v)-E*v)),'spin':total_spin_check(s,v)};Mraw=gram(B);A=np.vstack([g.conj()@B for g in G]) if G else np.empty((0,B.shape[1]),complex);elastic=gram(A)
 for g,a in zip(G,A):B-=g[:,None]*a[None,:]
 del s,v,G;gc.collect();HB=np.column_stack([t.mv(b)-E*b for b in B.T]);M0=gram(B);M1=gram(B,HB);M2=gram(HB);err=float(la.norm(Mraw-elastic-M0));assert err<1e-8
 checks={'source':sourcecheck,'projection_decomposition_error':err,'moment_Hermiticity':[float(la.norm(M-M.conj().T)) for M in [M0,M1,M2]],'moment_minimum_eigenvalues':[float(la.eigvalsh((M+M.conj().T)/2)[0]) for M in [M0,M1,M2]],'seed_spin':[]}
 for a in ([4,5] if j==1 else [12,13]):
  if la.norm(B[:,a])>1e-10:checks['seed_spin'].append(total_spin_check(t,B[:,a]/la.norm(B[:,a])))
 r={'name':name,'j':j,'Q':Q.tolist(),'cell_transfer':np.asarray(k).tolist(),'target_q':t.q.tolist(),'source_partner':partner,'dimension':t.nb,'E0':E,'channel_names':names,'Mraw':pack(Mraw),'M0_elastic':pack(elastic),'M0':pack(M0),'M1':pack(M1),'M2':pack(M2),'validation':checks,'H_calls':t.calls,'wall_seconds':time.perf_counter()-start,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
 np.savez(ROOT/(name+'_seeds.npz'),B=B,E0=E,q=t.q,j=j,Q=Q,elastic_amplitudes=A,Mraw=Mraw,M0=M0,M1=M1,M2=M2);r.update(wall_seconds=time.perf_counter()-start,CPU_seconds=time.process_time()-cpu_start,wall_scope='source setup, seeds, direct moments, validation, and seed checkpoint write');dump(ROOT/(name+'_moments.json'),r);print(json.dumps({k:r[k] for k in ['name','dimension','validation','wall_seconds','peak_RSS_MiB']}),flush=True)

def solve(name,count=10,ncv=64,maxiter=1500,tol=2e-12):
 if name=='j09_q20':
  # Compatibility with the in-flight queue: retain the measured seeds and split targets.
  import subprocess,sys
  old=ROOT/(name+'_seeds.npz');new=ROOT/'j09_q_seeds.npz'
  if old.exists():old.replace(new)
  oldmeta=ROOT/(name+'_moments.json')
  if oldmeta.exists():
   rec=json.loads(oldmeta.read_text());rec['name']='j09_q';dump(ROOT/'j09_q_moments.json',rec);oldmeta.unlink()
  for parity in ['1','-1']:subprocess.run([sys.executable,str(ROOT/'reduced_experiment.py'),'j09_q','--reflection',parity,'--count','6'],check=True)
  return
 start=time.perf_counter();z=np.load(ROOT/(name+'_seeds.npz'));j=float(z['j']);m=Model(primitive=j==1,jp=j);m.sector(z['q']);B=z['B'];E0=float(z['E0']);previous=ROOT/(name+'_poles.npz');G=[]
 if j==1:
  g=np.load(ROOT/'source_j1.npz')
  for q,v in [(g['q'],g['vector']),((-g['q'])%1,g['vector'].conj())]:
   if np.max(abs(np.exp(2j*np.pi*(m.gs@(m.q-q)))-1))<1e-8:G.append(v/la.norm(v))
 elif la.norm(z['Q'])<1e-10:
  c=CyclicModel(jp=j,angular=3);g=np.load(ROOT/f'inputs/p{j:.9f}_n00_r3_e0.npz');v=convert_basis(c,m,g['vector']);G=[v/la.norm(v)];del c
 op=monitored(m,name+'_solve')
 def action(v):
  y=op@v
  for g in G:y+=25*g*np.vdot(g,v)
  return y
 Aop=LinearOperator((m.nb,m.nb),matvec=action,dtype=complex)
 # Pure scalar source mixture, reproducible and nonzero in all accessible face channels.
 rng=np.random.default_rng(20261116);v0=B@(rng.normal(size=B.shape[1])+1j*rng.normal(size=B.shape[1]));v0/=la.norm(v0)
 complete=True
 try:vals,V=eigsh(Aop,k=count,which='SA',v0=v0,ncv=ncv,tol=tol,maxiter=maxiter)
 except ArpackNoConvergence as exc:vals,V=exc.eigenvalues,exc.eigenvectors;complete=False
 order=np.argsort(vals);vals=vals[order];V=V[:,order];poles=[];Rs=[];amplitudes=[]
 for i,e in enumerate(vals):
  v=V[:,i];v/=la.norm(v);h=m.mv(v);energy=float(np.vdot(v,h).real);res=float(la.norm(h-energy*v));spin=total_spin_check(m,v);a=v.conj()@B;R=np.outer(a.conj(),a);verified=res<1e-8 and spin['S2_eigen_residual']<1e-7;scalar=abs(spin['nearest_half_integer_spin']-.5)<1e-8;d={'index':i,'energy':energy,'omega':energy-E0,'H_residual':res,'spin':spin,'verified_eigenpair':verified,'scalar_accessible':scalar,'total_seed_weight':float(np.trace(R).real),'residue':pack(R),'amplitudes':pack(a),'status':'VERIFIED_SCALAR_POLE' if verified and scalar else ('VERIFIED_SPIN_FORBIDDEN_STATE' if verified else 'UNCONVERGED_NOT_PHYSICAL_POLE')};poles.append(d);Rs.append(R);amplitudes.append(a);print(json.dumps({k:v for k,v in d.items() if k not in ['residue','amplitudes']}),flush=True)
 np.savez_compressed(previous,vectors=V,energy=vals,q=m.q,residues=np.array(Rs),amplitudes=np.array(amplitudes));sums=[]
 for p in range(3):sums.append(sum((r['omega']**p*R for r,R in zip(poles,Rs) if r['verified_eigenpair'] and r['scalar_accessible']),np.zeros((B.shape[1],B.shape[1]),complex)))
 result={'name':name,'j':j,'dimension':m.nb,'E0':E0,'requested_eigenpairs':count,'ncv':ncv,'arpack_returned_all':complete,'poles':poles,'partial_moments':[pack(x) for x in sums],'remainder_moments':[pack(z[f'M{p}']-sums[p]) for p in range(3)],'H_calls_including_spin_checks':m.calls,'H_seconds_including_spin_checks':m.mv_seconds,'wall_seconds':time.perf_counter()-start,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024};dump(ROOT/(name+'_spectrum.json'),result);print(json.dumps({k:v for k,v in result.items() if k not in ['poles','partial_moments','remainder_moments']}),flush=True)

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('task',choices=['source','seeds','solve']);p.add_argument('--name',default='j1_G');p.add_argument('--j',type=float,default=1);p.add_argument('--Q',nargs=2,type=float,default=[0,0]);p.add_argument('--partner',action='store_true');p.add_argument('--count',type=int,default=10);p.add_argument('--ncv',type=int,default=64);p.add_argument('--maxiter',type=int,default=1500);a=p.parse_args()
 if a.task=='source':prepare_source()
 elif a.task=='seeds':construct(a.name,a.j,a.Q,a.partner)
 else:solve(a.name,a.count,a.ncv,a.maxiter)
