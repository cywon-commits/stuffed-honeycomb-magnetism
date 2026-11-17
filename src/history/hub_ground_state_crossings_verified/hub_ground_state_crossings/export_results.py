"""Export measured levels, certified sample classifications, branch gaps and static data."""
import pathlib,json,csv,numpy as np
import cross_scan as cs
from refine_crossings import summarize,all_j,root_record,FAMILIES
ROOT=cs.ROOT

def write_csv(name,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for row in rows for k in row))
 with (ROOT/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
  for row in rows:w.writerow({k:json.dumps(v) if isinstance(v,(list,dict,tuple)) else v for k,v in row.items()})
def build():
 levels=[];raw=[];parities=json.loads((ROOT/'reflection_labels.json').read_text()) if (ROOT/'reflection_labels.json').exists() else {}
 profile=ROOT/'branch_symmetry_and_overlaps.json'
 if profile.exists():parities.update({d['name']:d for d in json.loads(profile.read_text())})
 for p in ROOT.glob('p*_e*.json'):
  d=json.loads(p.read_text())
  if 'name' in d:raw.append(d)
 for j in [.2,.5]:
  for n,r in FAMILIES:raw.append(root_record(j,n,r))
 for d in raw:
  if d is None:continue
  phase=parities.get(d['name'],{})
  if d['cell_n']==[0,1]:phase={'reflection_parity':1 if d['character']==0 else -1}
  v=d.get('slope_dE_dJp');levels.append({'name':d['name'],'Jp_over_J':d['Jp'],'cell_n':d['cell_n'],'cell_k':[n/3 for n in d['cell_n']],'character':d['character'],'physical_root_index':d.get('physical_root_index',0),'spin':d.get('spin'),'reflection_parity':phase.get('reflection_parity'),'energy_J':d.get('energy'),'direct_residual_J':d.get('residual'),'S2_residual':d.get('spin_check',{}).get('S2_eigen_residual'),'verified':d.get('verified',False),'source':'REUSED_ENDPOINT' if d.get('reference_endpoint') else 'MEASURED_NEW','V_expectation':v,'hub_rim_bond_mean':v/54 if v is not None else None,'rim_bond_mean':(d['energy']-d['Jp']*v)/27 if v is not None else None,'dimension':d.get('dimension'),'H_calls':d.get('H_calls'),'wall_seconds':d.get('wall_seconds'),'setup_seconds':d.get('setup_seconds'),'peak_rss_MiB':d.get('peak_rss_MiB'),'ncv':d.get('ncv')})
 levels.sort(key=lambda d:(d['Jp_over_J'],d['cell_n'],d['character'],d['physical_root_index']));write_csv('all_levels.csv',levels)
 ground=[];sectors=[];points=[];bounds=[]
 for j in all_j():
  d=summarize(j);points.append(d)
  if not d or not d['global_ground_identified']:continue
  for item in d['family_evidence']:
   row={'Jp_over_J':j,'cell_n':item['family'][0],'character':item['family'][1],'kind':item['kind']}
   if item['kind']=='CONCAVITY_LOWER_BOUND':row.update(energy_lower_bound_J=item['energy_lower_bound'],bound_endpoint_couplings=item['bounds_endpoint_couplings'])
   else:row.update(energy_J=item['record']['energy'],spin=item['record']['spin'],residual_J=item['record']['residual'])
   bounds.append(row)
  w=d['winner'];higher=[x['energy'] for x in raw if x and x.get('verified') and x['Jp']==j and x['energy']>w['energy']+1e-8];v=w.get('slope_dE_dJp')
  if v is None and j in [.2,.5]:v=json.loads((cs.OLD/f'static_j{j:g}.json').read_text())['summary']['hub_rim_energy_from_Fourier']/j
  ground.append({'Jp_over_J':j,'name':w['name'],'cell_n':w['cell_n'],'cell_k':[n/3 for n in w['cell_n']],'character':w['character'],'spin':w['spin'],'reflection_parity':w.get('reflection_parity'),'energy_J':w['energy'],'direct_residual_J':w['residual'],'S2_residual':w['spin_check']['S2_eigen_residual'],'next_measured_distinct_energy_separation':min(higher)-w['energy'] if higher else None,'separation_note':'Not a complete excitation gap; symmetry/spin partners and unsolved same-block excited levels are separate.','V_expectation':v,'hub_rim_bond_mean':v/54 if v is not None else None,'rim_bond_mean':(w['energy']-j*v)/27 if v is not None else None})
  for key,s in d['sector_minima'].items():
   best=s['lowest_measured'];sectors.append({'Jp_over_J':j,'cell_n':s['cell_n'],'full_cell_sector_dimension':2228700,'minimum_identified':s['minimum_identified'],'energy_lower_bound_J':s['energy_lower_bound'],'energy_upper_from_verified_eigenpair_J':best['energy'] if best else None,'spin_of_lowest_measured_pair':best['spin'] if best else None,'character_of_lowest_measured_pair':best['character'] if best else None,'direct_residual_J':best['residual'] if best else None,'kind':'VERIFIED_CELL_MINIMUM' if s['minimum_identified'] else 'BOUNDED_CELL_MINIMUM_NOT_EXPLICITLY_SOLVED'})
 write_csv('ground_samples.csv',ground);write_csv('cell_sector_comparison.csv',sectors);write_csv('character_evidence.csv',bounds)
 obs=[]
 for p in ROOT.glob('obs_p*_e*.json'):obs.append(json.loads(p.read_text())['summary'])
 for j in [.2,.5]:
  a=json.loads((cs.OLD/f'static_j{j:g}.json').read_text())['summary'];obs.append({'name':f'reference_j{j:g}','Jp':j,'role':'GROUND_REUSED_ENDPOINT','spin':a['S_total'],'cell_n':[int(round(3*k))%3 for k in a['cell_ground_momentum']],'energy':a['energy'],'residual':a['energy_residual'],'rim_neel_m2':a['rim_neel_m2'],'hub_S_Gamma':a['hub_S_Gamma'],'hub_S_cell01':a['hub_S_cell01'],'hub_S_cell12':a['hub_S_cell12'],'rim_hub_uniform_Gamma':a['rim_hub_uniform_Gamma']})
 obs.sort(key=lambda x:(x['Jp'],x['role']));write_csv('static_observables.csv',obs)
 intervals=json.loads((ROOT/'crossing_intervals.json').read_text()) if (ROOT/'crossing_intervals.json').exists() else None
 if intervals:
  rows=[]
  for index,c in enumerate(intervals['crossings'],1):
   row={k:v for k,v in c.items() if k!='competing_branch_records'};row['index']=index
   if 'gap_right_minus_left_at_a' in c and 'gap_right_minus_left_at_b' in c:
    ga=c['gap_right_minus_left_at_a'];gb=c['gap_right_minus_left_at_b'];row['linear_interpolation_estimate']=(c['a']*gb-c['b']*ga)/(gb-ga) if gb!=ga else None;row['estimate_note']='Interpolation estimate only; the measured bracket is authoritative.'
   rows.append(row)
  write_csv('crossing_brackets.csv',rows)
 cost={'new_verified_eigenpairs':sum(bool(d.get('verified')) for d in raw if d and not d.get('reference_endpoint')),'new_incomplete_eigenpairs':sum(not d.get('verified') for d in raw if d and not d.get('reference_endpoint')),'measured_solve_wall_seconds_sum':sum(d.get('wall_seconds',0) for d in raw if d and not d.get('reference_endpoint')),'measured_setup_seconds_sum':sum(d.get('setup_seconds',0) for d in raw if d and not d.get('reference_endpoint')),'measured_H_calls_sum':sum(d.get('H_calls',0) for d in raw if d and not d.get('reference_endpoint')),'process_peak_rss_MiB':max([d.get('peak_rss_MiB',0) for d in raw if d]+[0]),'cost_scope':'Sum of instrumented new eigenpair solves, including verification/checkpoint work; excludes reused endpoint solves, endpoint profiling, plotting and static-observable runs. Peak RSS is process high water with cached bases.','CPU_quota':2,'OpenMP_threads':2,'BLAS_threads':1,'RAM_cap_GiB':8,'static_analysis_seconds_sum':sum(d.get('analysis_seconds',0) for d in obs),'static_pair_kernel_seconds_sum':sum(d.get('pair_kernel_seconds',0) for d in obs),'unrecorded_work_note':'Endpoint profiling and a discarded optional audit startup at j=0.209375 are excluded from these sums; no eigenvalue from that discarded startup is reported.'}
 (ROOT/'computational_cost.json').write_text(json.dumps(cost,indent=2));(ROOT/'exported_points.json').write_text(json.dumps(points,indent=2));print(json.dumps(cost,indent=2));return ground,sectors,obs,intervals
if __name__=='__main__':build()
