"""Prevent accidental observable/ground labels that contradict verified conventions."""
import pathlib,json
from analyze_state import analyze
from observables import coefficients
from coverage import verified_ground
ROOT=pathlib.Path(__file__).resolve().parent
result={}
for channel in ['physical_K','physical_Kprime']:
 try:coefficients([1/3,0],channel)
 except ValueError:result[channel+'_rejects_nonzero_cell_momentum']=True
 else:raise AssertionError(channel)
try:analyze(.7,ROOT/'isotropic_primitive_ground.npz')
except ValueError:result['incomplete_sector_coverage_rejected']=True
else:raise AssertionError('Incomplete coverage')
if verified_ground(.2) is not None:
 try:analyze(.2,ROOT/'j0.2_n01_r0.npz')
 except ValueError as e:
  assert 'does not match' in str(e);result['non_ground_checkpoint_rejected']=True
 else:raise AssertionError('Non-ground checkpoint')
else:result['non_ground_checkpoint_check']='SKIPPED: requires completed coupling 0.2 coverage.'
(ROOT/'driver_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2));print('DRIVER GUARD TESTS PASSED')
