import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import time,json,resource,numpy as np,scipy.linalg as la
from experiment import ROOT,pack,dump,Model
from observables import convert_basis
from native import seeds
from local_geometry import matrices,P
start=time.perf_counter();z=np.load(ROOT/'source_j1.npz');p=Model();p.sector(z['q']);s=Model(primitive=False);s.sector([0,0]);v=convert_basis(p,s,z['vector']);del p;v/=la.norm(v);vprime=v.conj();assert abs(np.vdot(v,vprime))<1e-10;out,rho=seeds(s,s,v,[0,0],rdm=True);A=np.vstack([out@v.conj(),out@vprime.conj()]);del out;ops=matrices();rows=[];allG=[]
for mu in range(6):
 G=np.array([[np.trace(rho[mu]@ops[mu,a].conj().T@ops[mu,b]) for b in range(3)] for a in range(3)]);cols=[mu,6+mu,12+mu];G0=A[:,cols].conj().T@A[:,cols]/9;Gp=G-G0;ch=float(np.trace(rho[mu]@ops[mu,2]).real);pD=float(np.trace(rho[mu]@P).real);rhs=3*np.sqrt(3)*ch;errors={'total_circular_asymmetry':float(abs(G[0,0]-G[1,1]-rhs)),'positive_frequency_asymmetry_with_zero_weight_removed':float(abs(Gp[0,0]-Gp[1,1]-(rhs-G0[0,0]+G0[1,1]))),'chi_squared':float(abs(G[2,2]-3*pD/16)),'dimer_symmetric':float(abs((G[0,0]+G[1,1])/2-9*pD/8))};assert max(errors.values())<1e-10 and la.eigvalsh(Gp)[0]>-1e-9;rows.append({'orientation':mu,'P_doublet':pD,'chi_CCW':ch,'M0_total':pack(G),'M0_ground':pack(G0),'M0_positive':pack(Gp),'errors':errors});allG.append(G)
np.savez_compressed(ROOT/'fresh_local_sum_rule_data.npz',rho_ABC=rho,ground_cell_amplitudes=A,local_total_Gram=np.array(allG));dump(ROOT/'fresh_local_sum_rules.json',{'scope':'Fresh refined source, direct local RDM and all verified ground-pair amplitudes; exact integrated sum rules, not a resolved all-Q excitation spectrum. Local circular identity includes omega=0. Kprime reverses chirality, orbital-pair average cancels its mean.','rows':rows,'wall_seconds':time.perf_counter()-start,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024});print('PASS fresh local identities',max(max(r['errors'].values()) for r in rows),'seconds',time.perf_counter()-start,flush=True)
