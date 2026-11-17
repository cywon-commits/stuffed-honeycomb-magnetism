import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import ctypes as ct,subprocess,time,resource,json,math,numpy as np,scipy.linalg as la
from geometry import make,ROOT
P=ct.c_void_p
def ptr(a):return a.ctypes.data_as(P)
def dims(g,particles):
 traces=[]
 for p in g['perms']:
  seen=set();cycles=[]
  for a in range(g['N']):
   if a in seen:continue
   k=a;l=0
   while k not in seen:seen.add(k);l+=1;k=int(p[k])
   cycles.append(l)
  poly=[0]*(particles+1);poly[0]=1
  for l in cycles:
   for k in reversed(range(l,particles+1)):poly[k]+=poly[k-l]
  traces.append(poly[particles])
 out=[]
 for q in g['momenta']:
  value=sum(t*np.exp(2j*np.pi*(R@q)) for t,R in zip(traces,g['cells']))/g['Ncell'];assert abs(value.imag)<1e-3 and abs(value.real-round(value.real))<1e-3;out.append({'k':q.tolist(),'dimension':int(round(value.real))})
 assert sum(x['dimension'] for x in out)==math.comb(g['N'],particles);return out

def run():
 g=make([[2,2],[-2,4]]);n=36;exact=dims(g,18);ram=int(open('/sys/fs/cgroup/memory.max').read());subprocess.run(['g++','-O3','-march=native','-shared','-fPIC','-fopenmp',str(ROOT/'kernel36pilot.cpp'),'-o',str(ROOT/'kernel36pilot.so')],check=True);lib=ct.CDLL(str(ROOT/'kernel36pilot.so'));lib.init36();lib.orbits36.argtypes=[ct.c_int,ct.c_int,ct.c_uint64,ct.c_int]+[P]*5;lib.orbits36.restype=ct.c_int;lib.mv36.argtypes=[ct.c_int,ct.c_int]+[P]*10+[ct.c_int];results=[];perms=np.ascontiguousarray(g['perms'],np.int32)
 for k in [4,6]:
  full=math.comb(n,k);assert full*16<ram/4;start=time.perf_counter();rid=np.full(full,-1,np.int32);shift=np.zeros(full,np.uint8);reps=np.empty(full,np.uint64);lens=np.empty(full,np.uint8);nb=lib.orbits36(n,k,full,12,*[ptr(a) for a in [perms,rid,shift,reps,lens]]);masks=reps[:nb].copy();lengths=lens[:nb].copy();del reps,lens;setup=time.perf_counter()-start;flips=np.array([(1<<a)|(1<<b) for a,b,t in g['bonds']],np.uint64);js=np.array([1 if t==0 else .3 for a,b,t in g['bonds']],float);diag=sum(j*(.25-.5*(((masks>>a)^(masks>>b))&1)) for (a,b,t),j in zip(g['bonds'],js));phase=np.ones(12,complex);rng=np.random.default_rng(36+k);x=rng.normal(size=nb)+1j*rng.normal(size=nb);x/=la.norm(x);y=np.empty_like(x)
  def action(v):
   out=np.empty_like(v);lib.mv36(nb,len(flips),*[ptr(a) for a in [masks,lengths,flips,js,diag,rid,shift,phase,v,out]],2);return out
  hx=action(x);u=rng.normal(size=nb)+1j*rng.normal(size=nb);u/=la.norm(u);herm=abs(np.vdot(u,hx)-np.vdot(action(u),x));assert herm<1e-10;ts=[]
  for _ in range(5):t=time.perf_counter();action(x);ts.append(time.perf_counter()-t)
  results.append({'N':36,'r':.3,'number_of_down_spins':k,'Sz_abs':18-k,'full_sector_dim':full,'translation_Gamma_dim':nb,'basis_generation_seconds':setup,'median_H_action_seconds':float(np.median(ts)),'H_action_samples_seconds':ts,'Hermiticity_error':float(herm),'metadata_bytes_measured_arrays':rid.nbytes+shift.nbytes+masks.nbytes+lengths.nbytes+diag.nbytes,'complex_vector_bytes':x.nbytes,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'status':'MEASURED_HIGH_SZ_PILOT; NOT_A_GROUND_STATE_CALCULATION'});print(results[-1],flush=True)
 low=min(x['dimension'] for x in exact);high=max(x['dimension'] for x in exact);big=results[-1];estimate=big['median_H_action_seconds']*high/big['translation_Gamma_dim'];out={'status':'PHASE_C_LOW_SPIN_ED_BLOCKED_BY_MEMORY; REAL_36SPIN_HIGH_SZ_PILOTS_COMPLETED','available_RAM_bytes':ram,'exact_Sz0_full_dimension':math.comb(36,18),'exact_translation_sector_dimensions':exact,'smallest_complex_vector_GiB':low*16/2**30,'largest_complex_vector_GiB':high*16/2**30,'Gamma_real_vector_GiB':exact[0]['dimension']*8/2**30,'dense_combinatorial_lookup_GiB':math.comb(36,18)*5/2**30,'measured_pilots':results,'translation_only_Sz0_H_seconds_estimate_from_high_Sz_pilot':estimate,'estimate_limitations':'Dimension-linear extrapolation from high magnetization; not a measured low-spin matvec. Half filling has more opposite-spin bonds and rank operations, so this can be optimistic. No converged low-spin 36-spin eigenvalue exists.','hypothetical_300_to_1000_H_actions_hours_at_two_CPU_estimate':[estimate*300/3600,estimate*1000/3600],'HPC_recommendation':'At least 128 GiB for selected Gamma one-dimensional point-group sectors with compatible indexing and restarted eigensolver; approximately 256 GiB and 64-128 CPU cores for meaningful Gamma/corner coverage. These are design estimates, not completed benchmarks.'};(ROOT/'feasibility36.json').write_text(json.dumps(out,indent=2))
 import csv
 with (ROOT/'results_36.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=results[0].keys());w.writeheader();w.writerows(results)
if __name__=='__main__':run()
