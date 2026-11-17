import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,time,json,resource,numpy as np,scipy.linalg as la
from scipy.sparse.linalg import eigsh,LinearOperator,ArpackNoConvergence
from experiment import ROOT,pack,dump,gram,monitored,Model
from little_model import LittleModel
from observables import total_spin_check

def run(name,angular=0,reflection=-1,count=8):
 start=time.perf_counter();cpu_start=time.process_time();z=np.load(ROOT/(name+'_seeds.npz'));q=z['q'];j=float(z['j']);full=Model(primitive=j==1,jp=j);full.sector(q);m=LittleModel(q,angular,reflection,primitive=j==1,jp=j);B=m.project(full,z['B']);del full;E0=float(z['E0']);M0=gram(B);tag=name+f'_a{angular}_p{reflection}';print(json.dumps({'tag':tag,'dimension':m.nb,'projected_seed_M0_trace':float(np.trace(M0).real),'rot_order':m.rot_order,'mirror_index':m.mirror_index}),flush=True)
 if np.trace(M0).real<1e-12:dump(ROOT/(tag+'_empty.json'),{'M0':pack(M0),'dimension':m.nb});return
 G=[];g=np.load(ROOT/'source_j1.npz') if j==1 else None
 # Project the matching primitive ground state, if present.
 if j==1:
  prim=Model();prim.sector(q)
  for qg,vg in [(g['q'],g['vector']),((-g['q'])%1,g['vector'].conj())]:
   if np.max(abs(np.exp(2j*np.pi*(prim.gs@(q-qg)))-1))<1e-8:
    vv=m.project(prim,vg[:,None])[:,0]
    if la.norm(vv)>1e-6:G.append(vv/la.norm(vv))
  del prim
 elif la.norm(z['Q'])<1e-10:
  from cyclic_model import CyclicModel
  from observables import convert_basis
  small=CyclicModel(jp=j,angular=3);old=np.load(ROOT/f'inputs/p{j:.9f}_n00_r3_e0.npz');base=Model(primitive=False,jp=j);base.sector([0,0]);vg=convert_basis(small,base,old['vector']);vv=m.project(base,vg[:,None])[:,0]
  if la.norm(vv)>1e-6:G.append(vv/la.norm(vv))
  del small,base,vg

 op=monitored(m,tag)
 def action(v):
  y=op@v
  for g in G:y+=25*g*np.vdot(g,v)
  return y
 A=LinearOperator((m.nb,m.nb),matvec=action,dtype=complex);rng=np.random.default_rng(17);v0=B@(rng.normal(size=B.shape[1])+1j*rng.normal(size=B.shape[1]));v0/=la.norm(v0);complete=True;method='INDEPENDENT_ARPACK_WITH_DIRECT_RESIDUAL'
 partner_file=ROOT/(name+f'_a{6-angular}_pNone_poles.npz')
 if j!=1 and la.norm(z['Q'])<1e-10 and angular in [4,5] and reflection is None and partner_file.exists():
  old=np.load(partner_file);E=old['energy'];V=old['vectors'].conj();method='EXACT_CONJUGATE_PARTNER_FRESHLY_REVERIFIED'
 else:
  try:E,V=eigsh(A,k=count,which='SA',v0=v0,ncv=max(48,4*count),tol=(1e-12 if name=='j1_G' else 1e-11),maxiter=800)
  except ArpackNoConvergence as exc:E,V=exc.eigenvalues,exc.eigenvectors;complete=False
 order=np.argsort(E);E=E[order];V=V[:,order];poles=[];Rs=[];aa=[]
 for i,e in enumerate(E):
  v=V[:,i]/la.norm(V[:,i]);h=m.mv(v);energy=float(np.vdot(v,h).real);res=float(la.norm(h-energy*v));spin=total_spin_check(m,v);a=v.conj()@B;R=np.outer(a.conj(),a);ok=res<1e-8 and spin['S2_eigen_residual']<1e-7;scalar=spin['nearest_half_integer_spin']==.5;p={'index':i,'energy':energy,'omega':energy-E0,'H_residual':res,'spin':spin,'verified_eigenpair':ok,'scalar_accessible':scalar,'angular':angular,'reflection':reflection,'residue':pack(R),'total_seed_weight':float(np.trace(R).real),'status':'VERIFIED_SCALAR_POLE' if ok and scalar else ('VERIFIED_SPIN_FORBIDDEN_STATE' if ok else 'UNCONVERGED_NOT_PHYSICAL_POLE')};poles.append(p);Rs.append(R);aa.append(a);print(json.dumps({k:v for k,v in p.items() if k!='residue'}),flush=True)
 np.savez_compressed(ROOT/(tag+'_poles.npz'),vectors=V,energy=E,residues=np.array(Rs),amplitudes=np.array(aa));result={'name':name,'tag':tag,'dimension':m.nb,'angular':angular,'reflection':reflection,'rotation_order':m.rot_order,'mirror_generator_index':m.mirror_index,'projected_M0':pack(M0),'arpack_returned_all':complete,'method':method,'poles':poles,'H_calls_including_spin_checks':m.calls,'iterative_H_calls':m.calls-2*len(poles),'CPU_seconds':time.process_time()-cpu_start,'wall_seconds':time.perf_counter()-start,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024};dump(ROOT/(tag+'_spectrum.json'),result)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('name');p.add_argument('--angular',type=int,default=0);p.add_argument('--reflection',default='-1');p.add_argument('--count',type=int,default=8);a=p.parse_args();run(a.name,a.angular,None if a.reflection=='none' else int(a.reflection),a.count)
