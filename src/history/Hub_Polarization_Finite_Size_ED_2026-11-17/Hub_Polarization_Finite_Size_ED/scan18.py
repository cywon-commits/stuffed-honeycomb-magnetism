import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,time,resource,argparse,csv,numpy as np,scipy.linalg as la
from scipy.sparse.linalg import eigsh,LinearOperator
from geometry import ROOT,catalog
from ed18 import Basis
from observables import correlations,values,checks
SWEEP=[.1,.15,.2,.225,.25,.275,.3,.35,.4,.45,.5,.55,.6,.65,.675,.68,.7]
def dump(path,x):path.write_text(json.dumps(x,indent=2))
def write_csv(path,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with path.open('w') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
class Scan:
 def __init__(self,name,g):self.name=name;self.g=g;self.b=Basis(g);self.rows=[];self.start=time.perf_counter()
 def point(self,r,force=False):
  folder=ROOT/'checkpoints'/self.name;folder.mkdir(parents=True,exist_ok=True);label=f'r{r:.9f}';out=folder/(label+'.json')
  if out.exists() and not force:return json.loads(out.read_text())
  start=time.perf_counter();cpu=time.process_time();b=self.b;levels=[];saved={};history=[]
  for iq,q in enumerate(self.g['momenta']):
   b.sector(q);b.set_r(r);rng=np.random.default_rng(20261117+iq);v0=rng.normal(size=b.nb)+1j*rng.normal(size=b.nb);op=LinearOperator((b.nb,b.nb),matvec=b.mv,dtype=complex);tt=time.perf_counter();E,V=eigsh(op,k=5,which='SA',ncv=48,tol=3e-13,v0=v0,maxiter=1500);ix=np.argsort(E);E=E[ix];V=V[:,ix]
   # Joint resolution of H and S² near crossings, then H diagonalization within each S.
   done=set();resolved=[]
   for i,e in enumerate(E):
    if i in done:continue
    group=[k for k,en in enumerate(E) if abs(en-e)<1e-4 and k not in done];done.update(group);Q,_=la.qr(V[:,group],mode='economic');SQ=np.column_stack([b.spin_action(v) for v in Q.T]);sv,U=la.eigh(Q.conj().T@SQ);sp=np.rint((np.sqrt(np.maximum(0,1+4*sv))-1)/2).astype(int)
    for spin in sorted(set(sp)):
     Z=Q@U[:,sp==spin];HZ=np.column_stack([b.mv(v) for v in Z.T]);en,rot=la.eigh((Z.conj().T@HZ+HZ.conj().T@Z)/2);W=Z@rot
     for en0,vec in zip(en,W.T):resolved.append((en0,vec))
   resolved.sort(key=lambda x:x[0]);E=np.array([x[0] for x in resolved]);V=np.column_stack([x[1] for x in resolved])
   for i,v in enumerate(V.T):
    v=v/la.norm(v);h=b.mv(v);energy=float(np.vdot(v,h).real);res=float(la.norm(h-energy*v));S,S2,sres=b.spin(v);assert res<1e-8 and sres<1e-6,(self.name,r,q,energy,res,S2,sres)
    record={'shape':self.name,'r':r,'iq':iq,'k1':float(q[0]),'k2':float(q[1]),'index':i,'energy':energy,'S':S,'S2':S2,'H_residual':res,'S2_residual':sres,'dimension':b.nb};levels.append(record);saved[(iq,i)]=v.copy()
   history.append({'iq':iq,'H_calls':b.calls,'wall_seconds':time.perf_counter()-tt,'matvec_seconds':b.mv_seconds,'dimension':b.nb})
  ordered=sorted(levels,key=lambda x:x['energy']);winner=ordered[0];near=[x for x in ordered if x['energy']-winner['energy']<1e-8];v=saved[(winner['iq'],winner['index'])];b.sector([winner['k1'],winner['k2']]);b.set_r(r);G,psi,imag=correlations(b,v);obs=values(self.g,G,winner['S']);cc=checks(b,psi,G);energy_corr=sum((1 if t==0 else r)*G[a,c] for a,c,t in self.g['bonds']);cc.update(energy_from_G_error=abs(energy_corr-winner['energy']),max_pair_imaginary=imag);assert max(cc.values())<1e-7 and obs['spin_sum_rule_error']<1e-7
  distinct=[x for x in ordered if (x['S']!=winner['S'] or abs(x['k1']-winner['k1'])+abs(x['k2']-winner['k2'])>1e-8) and x['energy']-winner['energy']>1e-8];gap=distinct[0]['energy']-winner['energy'] if distinct else None
  result={'shape':self.name,'r':r,'winner':winner,'degenerate_ground_states':near,'global_ground_coverage':'All allowed cell momenta in Sz=0. Every integer-S multiplet contains M=0; translations commute with SU(2). No total-spin multiplet omitted. Cold independent starts in each momentum.','levels':levels,'observables':obs,'validation':cc,'gap_to_lowest_distinct_spin_or_momentum':gap,'cost':{'wall_seconds':time.perf_counter()-start,'CPU_seconds':time.process_time()-cpu,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'sectors':history}}
  np.savez_compressed(folder/(label+'.npz'),vector=v,q=np.array([winner['k1'],winner['k2']]),energy=winner['energy'],Gdot=G,T=self.g['T'],r=r,S=winner['S']);dump(out,result);print(json.dumps({'shape':self.name,'r':r,'E':winner['energy'],'S':winner['S'],'q':[winner['k1'],winner['k2']],'pH':obs['p_H'],'wall_seconds':result['cost']['wall_seconds']}),flush=True);return result
 def run(self):
  points={r:self.point(r) for r in SWEEP};brackets=[]
  def signature(p):return {(x['S'],round(x['k1'],7),round(x['k2'],7)) for x in p['degenerate_ground_states']}
  # Refine every observed change; recursive midpoint screens for additional intervening sectors.
  queue=[(a,b) for a,b in zip(SWEEP,SWEEP[1:]) if signature(points[a])!=signature(points[b])]
  while queue:
   a,b=queue.pop(0);A=points[a];B=points[b]
   if b-a<=.0025:
    brackets.append({'shape':self.name,'r_lo':a,'r_hi':b,'width':b-a,'S_lo':A['winner']['S'],'S_hi':B['winner']['S'],'k_lo':json.dumps([[x['k1'],x['k2']] for x in A['degenerate_ground_states']]),'k_hi':json.dumps([[x['k1'],x['k2']] for x in B['degenerate_ground_states']]),'status':'VERIFIED_ENDPOINT_GROUND_SECTORS; scan-dependent crossing bracket'});continue
   m=round((a+b)/2,10);points[m]=self.point(m)
   if signature(A)!=signature(points[m]):queue.append((a,m))
   if signature(points[m])!=signature(B):queue.append((m,b))
  # Additional grid spacing <=0.0125 screens sampled intervals for missed plateaus.
  for r in np.arange(.1,.7000001,.0125):
   r=round(float(r),8)
   if r not in points:points[r]=self.point(r)
  # Refine changes discovered in the denser grid too.
  xs=sorted(points);new=[(a,b) for a,b in zip(xs,xs[1:]) if signature(points[a])!=signature(points[b]) and b-a>.0025]
  while new:
   a,b=new.pop();m=round((a+b)/2,10);points[m]=self.point(m)
   for u,v in [(a,m),(m,b)]:
    if signature(points[u])!=signature(points[v]):
     if v-u>.0025:new.append((u,v))
  brackets=[];xs=sorted(points)
  for a,b in zip(xs,xs[1:]):
   A,B=points[a],points[b]
   if signature(A)!=signature(B):brackets.append({'shape':self.name,'r_lo':a,'r_hi':b,'width':b-a,'S_lo':A['winner']['S'],'S_hi':B['winner']['S'],'k_lo':json.dumps(sorted(signature(A))),'k_hi':json.dumps(sorted(signature(B))),'status':'VERIFIED_GROUND_ENDPOINTS; numerical energy-intersection refinement follows'})
  dump(ROOT/(self.name+'_crossings.json'),brackets);return points
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--shape',choices=['K18','C18'],required=True);p.add_argument('--single',type=float);a=p.parse_args();g=catalog()[a.shape];s=Scan(a.shape,g)
 if a.single is None:s.run()
 else:s.point(a.single)
