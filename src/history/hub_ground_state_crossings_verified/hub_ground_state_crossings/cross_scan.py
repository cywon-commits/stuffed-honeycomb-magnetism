"""Cold lowest-root scans, physical eigenvector deflation, and verified affine-energy slopes."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,time,resource,numpy as np,scipy.linalg as la
from scipy.sparse.linalg import LinearOperator,eigsh,ArpackNoConvergence
from cyclic_model import CyclicModel
from observables import total_spin_check
ROOT=pathlib.Path(__file__).resolve().parent;OLD=(ROOT/'endpoint_reference') if (ROOT/'endpoint_reference').exists() else ROOT.parent/'static';CACHE={}
PLAN=[((0,0),r) for r in [0,1,2,3]]+[((0,1),r) for r in [0,1]]+[((1,2),r) for r in [0,1]]
class TimeLimit(Exception):pass

def tag(j,n,r,index=0):return f'p{j:.9f}_n{n[0]}{n[1]}_r{r}_e{index}'
def get_model(j,n,r,L=3,up=14):
 key=(L,up,tuple(n),r);start=time.perf_counter();hit=key in CACHE
 if not hit:
  m=CyclicModel(L=L,up=up,jp=j,cell_n=n,angular=r);m.hub_flags=np.array([a%3==2 or b%3==2 for a,b,_ in m.bonds]);m.diag_rim=np.zeros(m.nb);m.diag_hub=np.zeros(m.nb)
  for a,b,_ in m.bonds:
   part=.25-.5*(((m.masks>>a)^(m.masks>>b))&1)
   if a%3==2 or b%3==2:m.diag_hub+=part
   else:m.diag_rim+=part
  CACHE[key]=m
 m=CACHE[key];m.js=np.where(m.hub_flags,j,1.).astype(float);m.diag=m.diag_rim+j*m.diag_hub;m.bonds=[(a,b,j if a%3==2 or b%3==2 else 1.) for a,b,_ in m.bonds];m.calls=0;m.mv_seconds=0.
 return m,hit,time.perf_counter()-start

def slope(m,v):
 oldjs,olddiag=m.js,m.diag
 try:
  m.js=m.hub_flags.astype(float);m.diag=m.diag_hub;w=m.mv(v);return float(np.vdot(v,w).real)
 finally:m.js,m.diag=oldjs,olddiag

def solve_one(j,n,r,index=0,ncv=64,seconds=900,initial=None,force=False,L=3,up=14):
 name=tag(j,n,r,index);path=ROOT/(name+'.json')
 if path.exists() and not force:
  prev=json.loads(path.read_text())
  if prev.get('verified') and (ROOT/(name+'.npz')).exists():return prev
 m,hit,setup=get_model(j,n,r,L,up);phase_real=np.max(abs(m.phase.imag))<1e-13;start=time.perf_counter();calls=0;history=[];Q=[]
 for idx in range(index):
  p=ROOT/(tag(j,n,r,idx)+'.json');assert p.exists() and json.loads(p.read_text()).get('verified')
  if not (ROOT/(tag(j,n,r,idx)+'.npz')).exists():solve_one(j,n,r,index=idx,ncv=ncv,seconds=seconds,L=L,up=up)
  q=np.load(ROOT/(tag(j,n,r,idx)+'.npz'))['vector'];assert q.size==m.nb
  for old in Q:assert abs(np.vdot(old,q))<1e-7
  Q.append(q)
 mu=1+3*L*L+6*L*L*j # exceeds full Hamiltonian spectral width for positive couplings.
 def matvec(x):
  nonlocal calls
  if time.perf_counter()-start>seconds:raise TimeLimit()
  z=x.astype(complex) if phase_real else x[:m.nb]+1j*x[m.nb:];y=m.mv(z)
  for q in Q:y+=mu*q*np.vdot(q,z)
  calls+=1
  if calls%50==0:
   record={'name':name,'H_calls':calls,'elapsed_seconds':time.perf_counter()-start};history.append(record);print(json.dumps(record),flush=True)
  return y.real.copy() if phase_real else np.concatenate([y.real,y.imag])
 dim=m.nb if phase_real else 2*m.nb;A=LinearOperator((dim,dim),matvec=matvec,dtype=float);rng=np.random.default_rng(2026)
 if initial is None:v0=rng.normal(size=dim)
 else:
  initial=np.asarray(initial);v0=initial.real.copy() if phase_real else np.concatenate([initial.real,initial.imag])
 for q in Q:
  z=v0.astype(complex) if phase_real else v0[:m.nb]+1j*v0[m.nb:];z-=q*np.vdot(q,z);v0=z.real.copy() if phase_real else np.concatenate([z.real,z.imag])
 v0/=la.norm(v0);record={'name':name,'Jp':j,'cell_n':list(n),'cell_k':(np.array(n)/L).tolist(),'L':L,'up':up,'character':r,'cyclic_order':m.order,'physical_root_index':index,'dimension':m.nb,'ncv':ncv,'seed':2026,'start_type':'independent_random' if initial is None else 'continuation','model_cache_hit':hit,'setup_seconds':setup,'arpack_tol':1e-11,'verified':False,'deflation_shift':mu if Q else 0,'history':history}
 v=None
 try:
  e,z=eigsh(A,k=1,which='SA',ncv=ncv,tol=1e-11,maxiter=3000,v0=v0);v=z[:,0].astype(complex) if phase_real else z[:m.nb,0]+1j*z[m.nb:,0]
 except ArpackNoConvergence:record['status']='NO_CONVERGED_RITZ_PAIR'
 except TimeLimit:record['status']='WALL_LIMIT_NO_VERIFIED_EIGENPAIR'
 hseconds=m.mv_seconds
 if v is not None:
  v/=la.norm(v);hv=m.mv(v);energy=float(np.vdot(v,hv).real);res=float(la.norm(hv-energy*v));orth=max([abs(np.vdot(q,v)) for q in Q],default=0.);checkstart=time.perf_counter();spin=total_spin_check(m,v);derivative=slope(m,v)
  record.update(energy=energy,residual=res,orthogonality_to_lower_roots=float(orth),spin_check=spin,spin=spin['nearest_half_integer_spin'],slope_dE_dJp=derivative,verification_seconds=time.perf_counter()-checkstart,verified=bool(res<1e-8 and spin['S2_eigen_residual']<1e-5 and orth<1e-7))
  record['status']='VERIFIED' if record['verified'] else 'FAILED_DIRECT_CHECK'
  np.savez_compressed(ROOT/(name+'.npz'),vector=v,energy=energy,q=m.q,jp=j,cell_n=n,angular=r,cyclic=True,spin=record['spin'])
 record.update(H_calls=calls,H_kernel_seconds=hseconds,wall_seconds=time.perf_counter()-start,peak_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
 temporary=path.with_name('.'+path.name+'.tmp');temporary.write_text(json.dumps(record,indent=2));os.replace(temporary,path);print(json.dumps({k:v for k,v in record.items() if k!='history'}),flush=True)
 return record

def old_record(j,n,r):
 raw=json.loads((OLD/f'j{j:g}_n{n[0]}{n[1]}_r{r}.json').read_text());return {'Jp':j,'energy':raw['energy'],'residual':raw['residual'],'verified':raw['converged']}
def points(n,r):
 out=[old_record(j,n,r) for j in [.2,.5]]
 for p in ROOT.glob(f'p*_n{n[0]}{n[1]}_r{r}_e0.json'):
  d=json.loads(p.read_text())
  if d.get('verified'):out.append(d)
 return sorted(out,key=lambda x:x['Jp'])
def lower_chord(j,n,r):
 pts=points(n,r)
 for a,b in zip(pts[:-1],pts[1:]):
  if a['Jp']<=j<=b['Jp']:
   w=(j-a['Jp'])/(b['Jp']-a['Jp']);value=(1-w)*(a['energy']-max(1e-9,5*a['residual']))+w*(b['energy']-max(1e-9,5*b['residual']));return float(value),[a['Jp'],b['Jp']]
 raise ValueError(j)

def certify_reflection():
 from sectors import point_permutations
 from ed import geometry,translations
 c,b=geometry(3,.35);gs,ts=translations(3,c,False);pg=point_permutations(3,.35);F=pg[6];R=pg[2];invF=np.argsort(F);invR=np.argsort(R);q=np.array([1/3,0]);lookup={tuple(p):i for i,p in enumerate(ts)};chars=np.exp(-2j*np.pi*(gs@q))
 preserves=all(abs(chars[lookup[tuple(F[t[invF]])]]-chars[i])<1e-9 for i,t in enumerate(ts));relation=np.array_equal(F[R[invF]],invR);assert preserves and relation
 (ROOT/'corner_character_equivalence.json').write_text(json.dumps({'weighted_graph_reflection_verified':True,'corner_translation_character_preserved':preserves,'F_R_F_inverse_equals_R_inverse':relation,'conclusion':'At cell (1/3,2/3), C3 characters 1 and 2 have identical spectra by an exact weighted-graph reflection; character 2 is symmetry-derived, not separately timed.'},indent=2))

def initial_scan():
 certify_reflection()
 for j in [.275,.35,.425]:
  results={}
  # Independent random starts for every explicitly solved lowest root.
  for n,r in [((0,0),3),((0,0),0),((1,2),0),((0,1),0)]:results[(n,r)]=solve_one(j,n,r)
  for n,r in [((0,0),1),((0,0),2),((0,1),1),((1,2),1)]:
   current=min(x['energy'] for (nn,rr),x in results.items() if nn==n and x.get('verified'));bound,ends=lower_chord(j,n,r)
   if bound>current+1e-8:
    entry={'Jp':j,'cell_n':n,'character':r,'kind':'EXCLUDED_BY_SECTOR_CONCAVITY_BOUND','energy_lower_bound':bound,'endpoint_couplings':ends,'computed_sector_upper_bound':current,'note':'No eigenvalue or spin is assigned to this uncomputed character.'};(ROOT/(tag(j,n,r)+'_bound.json')).write_text(json.dumps(entry,indent=2));print(json.dumps(entry),flush=True)
   else:results[(n,r)]=solve_one(j,n,r)
  # A second physical root in Gamma character 3 distinguishes competing spin branches.
  solve_one(j,(0,0),3,index=1)
  roots=[x for x in results.values() if x.get('verified')];winner=min(roots,key=lambda x:x['energy']);(ROOT/f'initial_j{j:g}.json').write_text(json.dumps({'Jp':j,'lowest_checked':winner,'sector_minima':{str(n):min([x for (nn,rr),x in results.items() if nn==n and x.get('verified')],key=lambda x:x['energy']) for n in [(0,0),(0,1),(1,2)]},'uncomputed_characters_checked_by_concavity_bounds':True},indent=2))
if __name__=='__main__':initial_scan()
