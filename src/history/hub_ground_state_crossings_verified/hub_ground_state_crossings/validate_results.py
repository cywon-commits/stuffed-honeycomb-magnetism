"""Validate completed result records, sample classifications and static identities."""
import pathlib,json,csv
ROOT=pathlib.Path(__file__).resolve().parent

def main():
 levels=[json.loads(p.read_text()) for p in ROOT.glob('p*_e*.json')];levels=[d for d in levels if 'name' in d]
 assert levels and all(d.get('verified') for d in levels)
 assert all(d.get('up',14)==14 and d.get('L',3)==3 for d in levels)
 # Early production records used the same hardcoded/default L=3, up=14 source.
 expected={(0,0,0):371808,(0,0,1):371130,(0,0,2):371700,(0,0,3):371232,(0,1,0):1119396,(0,1,1):1109304,(1,2,0):742900,(1,2,1):742900}
 assert all(d['dimension']==expected[(*d['cell_n'],d['character'])] for d in levels)
 assert max(d['residual'] for d in levels)<1e-8
 assert max(d['spin_check']['S2_eigen_residual'] for d in levels)<1e-5
 for j in [.275,.35,.425]:
  p=json.loads((ROOT/f'point_{j:.9f}.json').read_text());assert p['global_ground_identified'] and all(s['minimum_identified'] for s in p['sector_minima'].values())
 intervals=json.loads((ROOT/'crossing_intervals.json').read_text());assert len(intervals['crossings'])==2 and intervals['all_bracket_widths_met']
 for c in intervals['crossings']:assert c['width']<.005 and c['gap_right_minus_left_at_a']>0 and c['gap_right_minus_left_at_b']<0 and len(c['competing_branch_records'])==4
 static=[json.loads(p.read_text()) for p in ROOT.glob('obs_p*_e*.json')]
 assert len(static)>=10
 for d in static:
  v=d['validation'];assert v['energy_reconstruction_error']<1e-8 and v['minimum_C_eigenvalue']>-1e-8 and v['Hermiticity_max']<1e-8 and v['onsite_sum_rule_error']<1e-8 and v['slope_from_correlations_error']<1e-8
  assert d['summary']['neel_identity_error']<1e-8 and d['summary']['total_spin_partition_error']<1e-7
 summary={'status':'PASS','new_verified_eigenpairs':len(levels),'max_direct_H_residual_J':max(d['residual'] for d in levels),'max_S2_eigen_residual':max(d['spin_check']['S2_eigen_residual'] for d in levels),'static_states_validated':len(static),'maximum_static_energy_reconstruction_error_J':max(d['validation']['energy_reconstruction_error'] for d in static),'small_tests':['crossing_validation.json','complex_deflation_validation.json','solver_validation.json'],'no_unconverged_records_used':True,'legacy_metadata_note':'Early initial-scan records omit explicit L/up fields; L=3 and up=14 follow their production source defaults and the checked symmetry-block dimensions.','coupling_bracket_count':len(intervals['crossings'])}
 (ROOT/'final_validation.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
