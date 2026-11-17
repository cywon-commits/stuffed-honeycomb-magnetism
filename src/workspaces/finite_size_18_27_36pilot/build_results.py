import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,csv,numpy as np,pathlib,resource,platform,subprocess,hashlib
from geometry import ROOT,make
from observables import values
from scan18 import write_csv

def run():
 scan=[];sectors=[];levels=[];comparison=[];hub=[];cost=[]
 for shape in ['K18','C18']:
  for f in sorted((ROOT/'checkpoints'/shape).glob('*.json')):
   p=json.loads(f.read_text());w=p['winner'];o=p['observables'];row={'shape':shape,'N':18,'Ncell':6,'r':p['r'],'energy':w['energy'],'S':w['S'],'k1':w['k1'],'k2':w['k2'],'H_residual':w['H_residual'],'S2_residual':w['S2_residual'],'ground_sectors':json.dumps([[x['S'],x['k1'],x['k2']] for x in p['degenerate_ground_states']]),**{k:v for k,v in o.items() if not isinstance(v,(list,dict))},'legacy_normalized_Neel_mode':2*6*o['m_N2'],'gap_to_lowest_measured_distinct_spin_or_momentum':p['gap_to_lowest_distinct_spin_or_momentum'],'coverage':'All six allowed cell momenta; all integer S via SU2 M=0','source':'NEW_18_SPIN_ED'};scan.append(row);comparison.append(row)
   for iq,q in enumerate(make(np.load(f.with_suffix('.npz'))['T'])['momenta']):hub.append({'shape':shape,'N':18,'r':p['r'],'k1':q[0],'k2':q[1],'S_h_spin_summed':o['hub_structure_spin_summed'][iq],'S_h_zz_SU2_average':o['hub_structure_spin_summed'][iq]/3})
   for l in p['levels']:levels.append(l)
   for iq in range(6):
    loc=[l for l in p['levels'] if l['iq']==iq];first=min(loc,key=lambda l:l['energy']);sectors.append({**first,'number_of_returned_lowest_eigenpairs':len(loc),'M':0,'SU2_completeness':'Every integer S has M=0; translations commute with S².'})
   c=p['cost'];cost.append({'shape':shape,'r':p['r'],'wall_seconds':c['wall_seconds'],'CPU_seconds':c['CPU_seconds'],'peak_process_RSS_MiB':c['peak_RSS_MiB'],'H_calls':sum(x['H_calls'] for x in c['sectors']),'scope':'complete six-sector point and static correlators; excludes checkpoint compression'})
 for p in json.loads((ROOT/'reference27_results.json').read_text()):
  o=p['observables'];row={k:v for k,v in p.items() if k!='observables'};row.update({k:v for k,v in o.items() if not isinstance(v,(list,dict))});row['legacy_normalized_Neel_mode']=18*o['m_N2'];row['source']='VERIFIED_27_SPIN_DATA_REINGESTED';comparison.append(row)
  for q,s in zip(make([[3,0],[0,3]])['momenta'],o['hub_structure_spin_summed']):hub.append({'shape':'N27','N':27,'r':p['r'],'k1':q[0],'k2':q[1],'S_h_spin_summed':s,'S_h_zz_SU2_average':s/3})
 scan.sort(key=lambda x:(x['shape'],x['r']));comparison.sort(key=lambda x:(x['N'],x['shape'],x['r']));write_csv(ROOT/'gs_18_scan.csv',scan);write_csv(ROOT/'size_comparison.csv',comparison);write_csv(ROOT/'sector_completeness_18.csv',sectors);write_csv(ROOT/'all_low_levels_18.csv',levels);write_csv(ROOT/'hub_structure_factors.csv',hub);write_csv(ROOT/'computational_costs.csv',cost)
 validation={'status':'PASS','number_of_18spin_couplings_and_shapes':len(scan),'number_of_18spin_momentum_minima':len(sectors),'number_of_18spin_returned_eigenpairs':len(levels),'max_ground_H_residual':max(x['H_residual'] for x in scan),'max_all_eigenpair_H_residual':max(x['H_residual'] for x in levels),'max_ground_S2_residual':max(x['S2_residual'] for x in scan),'max_all_eigenpair_S2_residual':max(x['S2_residual'] for x in levels),'max_18spin_total_spin_sum_rule_error':max(x['spin_sum_rule_error'] for x in scan),'max_27spin_total_spin_sum_rule_error':max(x['observables']['spin_sum_rule_error'] for x in json.loads((ROOT/'reference27_results.json').read_text())),'max_18spin_energy_from_G_error':max(json.loads(f.read_text())['validation']['energy_from_G_error'] for shape in ['K18','C18'] for f in (ROOT/'checkpoints'/shape).glob('*.json')),'note':'Global ground-state coverage is at all sampled points. Numerical scans cannot rule out arbitrarily narrow re-entrant plateaus between unchanged endpoints.'}
 assert validation['max_ground_H_residual']<1e-8 and validation['max_ground_S2_residual']<1e-8 and validation['max_all_eigenpair_S2_residual']<1e-8 and validation['max_18spin_total_spin_sum_rule_error']<1e-7
 (ROOT/'validation.json').write_text(json.dumps(validation,indent=2));print(validation)
 # Remove only redundant expanded copies; vectors, quantum numbers and full G matrices remain.
 for f in (ROOT/'checkpoints').rglob('*.npz'):
  z=np.load(f)
  if 'full_vector' in z:
   data={k:z[k] for k in z.files if k!='full_vector'};z.close();np.savez_compressed(f,**data)
 env={'CPU_model':next(l.split(':',1)[1].strip() for l in pathlib.Path('/proc/cpuinfo').read_text().splitlines() if l.startswith('model name')),'CPU_quota':pathlib.Path('/sys/fs/cgroup/cpu.max').read_text().strip(),'RAM_limit_bytes':int(pathlib.Path('/sys/fs/cgroup/memory.max').read_text()),'Python':platform.python_version(),'NumPy':np.__version__,'OpenMP_threads':2,'BLAS_threads':1};(ROOT/('cached_validation_environment.json' if (ROOT/'environment.json').exists() else 'environment.json')).write_text(json.dumps(env,indent=2))
if __name__=='__main__':run()
