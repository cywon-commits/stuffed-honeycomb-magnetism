"""Independent static SU(2) raising check of the unexpected S=5/2 ground multiplet."""
import pathlib,json,numpy as np,scipy.linalg as la,time,resource
from ed import Model
from cyclic_model import CyclicModel
from observables import spin_seed,convert_basis
from coverage import verified_ground
ROOT=pathlib.Path(__file__).resolve().parent
jp=.5;cov=verified_ground(jp);assert cov is not None
z=np.load(ROOT/(cov['ground_name']+'.npz'));src=CyclicModel(jp=jp,cell_n=tuple(z['cell_n']),angular=int(z['angular']));source=Model(up=14,primitive=False,jp=jp);source.sector(z['q']);v=convert_basis(src,source,z['vector']);del src
t=time.perf_counter();out=[];E=cov['ground_energy']
for targetup,expected in [(15,8.),(16,5.),(17,0.)]:
 target=Model(up=targetup,primitive=False,jp=jp);target.sector(z['q']);w=3*spin_seed(source,target,v,[0,0],'total','+');norm2=float(np.vdot(w,w).real)
 row={'source_up':source.up,'target_up':targetup,'target_M':targetup-13.5,'target_cell_dimension':target.nb,'raising_norm_squared':norm2,'expected_norm_squared_for_S_5_over_2':expected}
 if expected:
  assert abs(norm2-expected)<1e-9;w/=np.sqrt(norm2);hw=target.mv(w);energy=float(np.vdot(w,hw).real);res=float(la.norm(hw-E*w));assert abs(energy-E)<1e-10 and res<1e-9
  row.update(energy=energy,energy_relative_to_original=energy-E,eigen_residual_at_original_energy=res)
 else:assert norm2<1e-18
 out.append(row);print(json.dumps(row),flush=True)
 source=target;v=w
result={'Jp':jp,'spin_verified':2.5,'checks':out,'analysis_seconds':time.perf_counter()-t,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'scope':'Static total-spin ladder identities and equal-energy multiplet check only; no spectrum or frequency-dependent calculation.'}
(ROOT/'spin_multiplet_validation.json').write_text(json.dumps(result,indent=2))
