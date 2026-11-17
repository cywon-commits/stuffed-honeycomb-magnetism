"""Validate complete momentum/little-group coverage before declaring a ground state."""
import pathlib,json
ROOT=pathlib.Path(__file__).resolve().parent
PLAN=[((0,0),ell) for ell in [0,1,2,3]]+[((0,1),ell) for ell in [0,1]]+[((1,2),ell) for ell in [0,1,2]]
def verified_ground(jp,root=ROOT):
 records=[]
 for n,ell in PLAN:
  name=f'j{jp:g}_n{n[0]}{n[1]}_r{ell}';p=root/(name+'.json')
  if not p.exists():return None
  r=json.loads(p.read_text())
  if not r['converged'] or r.get('residual',1)>1e-8:return None
  r['name']=name;records.append(r)
 by={(tuple(r['cell_n']),r['angular']):r for r in records}
 gamma=sum(by[((0,0),ell)]['dimension']*(2 if ell in [1,2] else 1) for ell in [0,1,2,3]);generic=sum(by[((0,1),ell)]['dimension'] for ell in [0,1]);corner=sum(by[((1,2),ell)]['dimension'] for ell in [0,1,2])
 assert gamma==generic==corner==2228700
 assert gamma+6*generic+2*corner==20058300
 best=min(records,key=lambda r:r['energy']);report={'Jp':jp,'all_representative_blocks_converged':True,'verified_blocks':len(records),'cell_sector_dimension_checks':[gamma,generic,corner],'full_dimension_covered':gamma+6*generic+2*corner,'ground_name':best['name'],'ground_energy':best['energy'],'ground_cell_n':best['cell_n'],'ground_cyclic_order':best['cyclic_order'],'ground_angular':best['angular'],'ground_residual':best['residual'],'lowest_energy_separation_between_checked_blocks':sorted(r['energy'] for r in records)[1]-best['energy'],'separation_note':'This is not the full first excitation gap; same-block excited roots were not calculated.','measured_solve_seconds_sum':sum(r['solve_seconds'] for r in records),'maximum_measured_peak_rss_MiB':max(r['peak_rss_MiB'] for r in records),'complex_conjugate_Gamma_partners':'ell=4 matches ell=2; ell=5 matches ell=1 for the real Hamiltonian at Gamma.'}
 (root/f'ground_j{jp:g}_coverage.json').write_text(json.dumps(report,indent=2));return report
if __name__=='__main__':
 for jp in [.2,.5]:print(jp,verified_ground(jp))
