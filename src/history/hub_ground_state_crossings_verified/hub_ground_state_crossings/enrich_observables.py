"""Sublattice moment identities and same-coupling branch contrasts."""
import pathlib,json,csv,numpy as np
from ed import geometry
ROOT=pathlib.Path(__file__).resolve().parent

def run():
 states={}
 for path in ROOT.glob('obs_p*_e*.json'):
  d=json.loads(path.read_text());s=d['summary'];z=d['C']['(0, 0)'];C=np.array(z['real'])+1j*np.array(z['imag']);u=np.array([1,1,0.]);h=np.array([0,0,1.]);s.update(SA2=float(9*C[0,0].real),SB2=float(9*C[1,1].real),SR2=float(9*np.vdot(u,C@u).real),SH2=float(9*C[2,2].real),SR_dot_SH=float(9*np.vdot(u,C@h).real));s['neel_identity_error']=float(abs(324*s['rim_neel_m2']-(2*(s['SA2']+s['SB2'])-s['SR2'])));s['total_spin_partition_error']=float(abs(s['spin']*(s['spin']+1)-(s['SR2']+s['SH2']+2*s['SR_dot_SH'])));assert s['neel_identity_error']<1e-8 and s['total_spin_partition_error']<1e-7
  G=np.load(path.with_name('obs_'+s['name']+'_realspace.npz'))['Gdot'];_,bonds=geometry(3,s['Jp']);rec=json.loads((ROOT/(s['name']+'.json')).read_text());d['validation']['slope_from_correlations_error']=float(abs(sum(G[a,b] for a,b,_ in bonds if a%3==2 or b%3==2)-rec['slope_dE_dJp']));assert d['validation']['slope_from_correlations_error']<1e-8
  path.write_text(json.dumps(d,indent=2));states[s['name']]=s
 file=ROOT/'crossing_intervals.json'
 if not file.exists():file=ROOT/'crossing_intervals_preview.json'
 if not file.exists():return
 intervals=json.loads(file.read_text());rows=[]
 for i,c in enumerate(intervals['crossings'],1):
  j=c['a'];A=c['left_quantum_numbers'];B=c['right_quantum_numbers'];ca=next((x for x in c['competing_branch_records'] if x['Jp']==j and x['cell_n']==A[0] and x['character']==A[1] and x['spin']==A[2]),None);cb=next((x for x in c['competing_branch_records'] if x['Jp']==j and x['cell_n']==B[0] and x['character']==B[1] and x['spin']==B[2]),None)
  if not ca or not cb or ca['name'] not in states or cb['name'] not in states:continue
  a,b=states[ca['name']],states[cb['name']];row={'crossing':i,'common_Jp':j,'left_branch_name':a['name'],'right_branch_name':b['name'],'energy_right_minus_left':b['energy']-a['energy'],'left_branch_role':a['role'],'right_branch_role':b['role']}
  for key in ['rim_neel_m2','hub_S_Gamma','hub_S_cell01','hub_S_cell12','rim_hub_uniform_Gamma','SA2','SB2','SR2','SH2','SR_dot_SH']:row[key+'_left']=a[key];row[key+'_right']=b[key];row['delta_'+key]=b[key]-a[key]
  row['delta_neel_from_sublattice_moment_change']=2*((b['SA2']+b['SB2'])-(a['SA2']+a['SB2']))/324;row['delta_neel_from_rim_total_spin_change']=-(b['SR2']-a['SR2'])/324;rows.append(row)
 if rows:
  with (ROOT/'same_coupling_branch_contrasts.csv').open('w') as f:
   writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
 print(json.dumps({'states_enriched':len(states),'same_coupling_contrasts':len(rows)}))
if __name__=='__main__':run()
