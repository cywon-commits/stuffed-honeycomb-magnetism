"""Use exact conjugation of verified primitive-K target states for K-transfer poles.
This changes the operator overlaps; they are recomputed against K-transfer seeds.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time,resource,numpy as np,scipy.linalg as la
from experiment import ROOT,Model,pack,dump
from observables import total_spin_check

def run():
 for name,checkpoint,conjugate in [('j1_K','j1_G_poles.npz',True),('j1_KprimeSource','inputs/Gamma_primitive.npz',False)]:
  p=ROOT/(name+'_seeds.npz')
  if not p.exists():continue
  start=time.perf_counter();z=np.load(p);m=Model();m.sector(z['q']);B=z['B'];E0=float(z['E0']);old=np.load(ROOT/checkpoint)
  if 'vectors' in old:V=old['vectors'].conj() if conjugate else old['vectors'];Es=old['energy']
  else:V=old['vector'][:,None];Es=np.array([float(old['energy'])])
  poles=[];amps=[];Rs=[]
  for i,v in enumerate(V.T):
   v/=la.norm(v);h=m.mv(v);E=float(np.vdot(v,h).real);res=float(la.norm(h-E*v));spin=total_spin_check(m,v);a=v.conj()@B;R=np.outer(a.conj(),a);ok=res<1e-8 and spin['S2_eigen_residual']<1e-7;scalar=spin['nearest_half_integer_spin']==.5;poles.append({'index':i,'energy':E,'omega':E-E0,'H_residual':res,'spin':spin,'verified_eigenpair':ok,'scalar_accessible':scalar,'total_seed_weight':float(np.trace(R).real),'residue':pack(R),'amplitudes':pack(a),'status':'VERIFIED_SCALAR_POLE' if ok and scalar else ('VERIFIED_SPIN_FORBIDDEN_STATE' if ok else 'UNCONVERGED_NOT_PHYSICAL_POLE')});amps.append(a);Rs.append(R)
  np.savez_compressed(ROOT/(name+'_poles.npz'),vectors=V,energy=Es,q=m.q,amplitudes=amps,residues=Rs);r={'name':name,'dimension':m.nb,'method':'REUSED_EXACT_CONJUGATE_OR_PRIOR_PRIMITIVE_TARGET_CHECKPOINT_FRESHLY_REVERIFIED','poles':poles,'H_calls_including_spin_checks':m.calls,'iterative_H_calls':0,'wall_seconds':time.perf_counter()-start,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024};dump(ROOT/(name+'_spectrum.json'),r);print(name,[(p['omega'],p['spin']['nearest_half_integer_spin'],p['total_seed_weight']) for p in poles],flush=True)
if __name__=='__main__':run()
