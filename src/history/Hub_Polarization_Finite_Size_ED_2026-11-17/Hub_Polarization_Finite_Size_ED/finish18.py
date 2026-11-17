import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,time,itertools,numpy as np,scipy.linalg as la
from scipy.sparse.linalg import eigsh,LinearOperator
from geometry import ROOT,catalog
from ed18 import Basis
from scan18 import Scan,write_csv
from observables import full_pair_action,correlations,values

def run():
 shapes=catalog();crossrows=[];fresh=[];sanity=[];variances=[]
 for name in ['K18','C18']:
  g=shapes[name];scanner=Scan(name,g)
  for br in json.loads((ROOT/(name+'_crossings.json')).read_text()):
   a,b=br['r_lo'],br['r_hi'];A=scanner.point(a);B=scanner.point(b);wa,wb=A['winner'],B['winner']
   def branch(P,label):
    candidates=[l for l in P['levels'] if l['S']==label['S'] and abs(l['k1']-label['k1'])+abs(l['k2']-label['k2'])<1e-8];assert candidates;return min(candidates,key=lambda l:l['energy'])
   aa,ab=branch(A,wa),branch(A,wb);ba,bb=branch(B,wa),branch(B,wb);da=ab['energy']-aa['energy'];db=bb['energy']-ba['energy'];err=max(x['H_residual'] for x in [aa,ab,ba,bb]);assert da>5*err and db<-5*err,(name,br,da,db,err)
   estimate=round(a+(b-a)*da/(da-db),10);P=scanner.point(estimate);pa,pb=branch(P,wa),branch(P,wb);delta=pb['energy']-pa['energy'];assert P['winner']['S'] in [wa['S'],wb['S']];crossrows.append({**br,'delta_E_hi_minus_lo_branch_at_rlo':da,'delta_E_hi_minus_lo_branch_at_rhi':db,'maximum_endpoint_H_residual':err,'secant_intersection_estimate':estimate,'actual_branch_gap_at_estimate':delta,'interpretation':'Finite-size branch crossing; endpoints verified with full momentum/SU2 coverage. Secant estimate is not an exact critical coupling.'})
  # Whole-spin-space control (all M=0 states), with an independent random start.
  for r in [.4,1.]:
   P=scanner.point(r);m=Basis(g,full=True);m.set_r(r);rng=np.random.default_rng(701);v0=rng.normal(size=m.nb)+1j*rng.normal(size=m.nb);E,V=eigsh(LinearOperator((m.nb,m.nb),matvec=m.mv,dtype=complex),k=2,which='SA',ncv=40,tol=1e-12,v0=v0);err=abs(min(E)-P['winner']['energy']);assert err<1e-8;sanity.append({'shape':name,'r':r,'independent_full_Sz0_dimension':m.nb,'full_vs_reduced_ground_energy_error':float(err),'full_Sz0_H_calls':m.calls})
  # Degenerate r=0: independent rim singlet times chosen maximal-hub-spin Dicke state.
  m=Basis(g,full=True);m.set_r(0);rim=[i for i in range(18) if i%3!=2];hub=list(range(2,18,3));rs=np.array([sum(1<<a for a in c) for c in itertools.combinations(range(12),6)],np.uint32);lookup={int(s):i for i,s in enumerate(rs)};H=np.zeros((len(rs),len(rs)))
  for col,s0 in enumerate(rs):
   s=int(s0)
   for a,c,t in g['bonds']:
    if t:continue
    u=ri=rim.index(a);v=rim.index(c);diff=((s>>u)^(s>>v))&1;H[col,col]+=.25-.5*diff
    if diff:H[lookup[s^((1<<u)|(1<<v))],col]+=.5
  es,us=la.eigh(H,subset_by_index=[0,0]);psi=np.zeros(m.nb,complex)
  for i,s in enumerate(rs):
   mask=sum(((int(s)>>u)&1)<<a for u,a in enumerate(rim))
   for hh in itertools.combinations(hub,3):psi[m.lookup[mask+sum(1<<a for a in hh)]]=us[i,0]/np.sqrt(20)
  res=float(la.norm(m.mv(psi)-es[0]*psi));S,S2,sres=m.spin(psi);assert res<1e-9 and S==3 and sres<1e-9;sanity.append({'shape':name,'r':0,'energy':float(es[0]),'H_residual':res,'S':S,'S2_residual':sres,'note':'Rim singlet times hub S=3,M=0 test vector; free hubs have a degenerate manifold. This does not select a polarized physical ground state at r=0.'})
  required=[.2,.3,.4,.5,.6,.68]+[b[x] for b in json.loads((ROOT/(name+'_crossings.json')).read_text()) for x in ['r_lo','r_hi']]
  for r in sorted(set(required)):
   P=scanner.point(r);z=np.load(ROOT/'checkpoints'/name/f'r{r:.9f}.npz');m=Basis(g);m.sector(z['q']);m.set_r(r);v=z['vector'];hv=m.mv(v);res=float(la.norm(hv-float(z['energy'])*v));S,S2,sres=m.spin(v);G,psi,im=correlations(m,v);gerr=float(la.norm(G-z['Gdot']));assert res<1e-8 and sres<1e-7 and gerr<1e-9;fresh.append({'shape':name,'r':r,'checkpoint_H_residual':res,'checkpoint_S2_residual':sres,'full_reconstruction_norm_error':float(abs(la.norm(psi)-1)),'recomputed_G_error':gerr})
   row={'shape':name,'r':r}
   for label,ix in [('hub',hub),('rim',rim)]:
    pairs=[(a,b,2.) for i,a in enumerate(ix) for b in ix[i+1:]];action=full_pair_action(m.states,m.lookup,psi,pairs,.75*len(ix));mean=float(np.vdot(psi,action).real);var=float(max(0,la.norm(action)**2-mean**2));row[label+'_S2']=mean;row[label+'_S2_variance']=var
   variances.append(row)
 write_csv(ROOT/'crossing_brackets_18.csv',crossrows);write_csv(ROOT/'checkpoint_reconstruction_validation.csv',fresh);write_csv(ROOT/'subsystem_spin_variances_18.csv',variances);(ROOT/'sanity_checks.json').write_text(json.dumps({'status':'PASS','rows':sanity},indent=2));print('PASS 18-spin endpoints/intersections, full-space controls, checkpoint reconstruction and r=0/r=1 sanity checks.',flush=True)
if __name__=='__main__':run()
