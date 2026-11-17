import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,json
from geometry import ROOT,make
from scan18 import write_csv
D=ROOT/'reference27';g=make([[3,0],[0,3]]);rows=[]
for r in json.loads((ROOT/'reference27_results.json').read_text()):
 n=r['source_file'];meta=n.replace('_realspace.npz','.json').replace('raw_observables_','observables_').replace('.npz','.json');data=json.loads((D/meta).read_text());old=data.get('summary',data);G=np.load(D/f"correlations_r{r['r']:.9f}.npz")['Gdot'];phase=np.exp(-2j*np.pi*(g['coords']@np.array([1/3,2/3])))/np.sqrt(27);mag=float(np.vdot(phase,G@phase).real)
 for key,new in [('rim_neel_m2',r['observables']['m_N2']),('hub_S_Gamma',r['observables']['hub_S_Gamma']),('physical_K_120_intensity',mag),('physical_K_intensity',mag)]:
  if key not in old:continue
  err=abs(float(old[key])-new);assert err<1e-7;rows.append({'r':r['r'],'benchmark':key,'archived_value':old[key],'recomputed_value':new,'absolute_error':err,'source_metadata':meta})
write_csv(ROOT/'reference_benchmark_comparison.csv',rows);(ROOT/'reference_benchmark_validation.json').write_text(json.dumps({'status':'PASS','number_of_benchmarks':len(rows),'maximum_error':max(r['absolute_error'] for r in rows)},indent=2));print('PASS',len(rows),'archived 27-spin correlation benchmarks.')
