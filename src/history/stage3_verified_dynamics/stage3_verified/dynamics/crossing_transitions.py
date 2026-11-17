import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,time,resource,numpy as np,scipy.linalg as la
from ed import Model
from cyclic_model import CyclicModel
from observables import convert_basis,spin_seed,total_spin_check,inversion_vector
from sectors import point_permutations
ROOT=pathlib.Path(__file__).resolve().parent;DATA=ROOT.parent/'crossing'
def pack(a):return {'real':np.asarray(a).real.tolist(),'imag':np.asarray(a).imag.tolist()}
def run(j):
 t=time.perf_counter();m=Model(primitive=False,jp=j);m.sector([0,0]);states=[]
 for r in [0,3]:
  name=f'p{j:.9f}_n00_r{r}_e0';d=json.loads((DATA/(name+'.json')).read_text());z=np.load(DATA/(name+'.npz'));c=CyclicModel(jp=j,cell_n=(0,0),angular=r);v=convert_basis(c,m,z['vector']);del c
  states.append((d,v));assert la.norm(m.mv(v)-d['energy']*v)<1e-8
 states.sort(key=lambda x:x[0]['energy']);(g,v),(ex,w)=states;amps=[];checks={};mom={}
 for ch in ['n','m','H']:
  b=spin_seed(m,m,v,[0,0],ch);amp=np.vdot(w,b);amps.append(amp);h=m.mv(b)-g['energy']*b;elastic=abs(np.vdot(v,b))**2
  mom[ch]={'M0':float(np.vdot(b,b).real),'M1':float(np.vdot(b,h).real),'M2':float(np.vdot(h,h).real),'elastic_weight':float(elastic),'inelastic_M0':float(np.vdot(b,b).real-elastic)}
 total=spin_seed(m,m,v,[0,0],'total');checks['total_z_conservation_seed_error']=float(la.norm(total-v/6));checks['m_h_opposite_branch_amplitude_max']=float(max(abs(amps[1]),abs(amps[2])))
 for d,x in states:
  F=point_permutations()[6];R=point_permutations()[1]
  checks[d['name']]={'H_residual':float(la.norm(m.mv(x)-d['energy']*x)),'spin':total_spin_check(m,x),'reflection_residual':float(la.norm(inversion_vector(m,x,F)-(1 if d['character']==0 else -1)*x))}
 a=np.array(amps);result={'Jp':j,'source_M':.5,'ground':g['name'],'excited':ex['name'],'source_spin':g['spin'],'excited_spin':ex['spin'],'omega':ex['energy']-g['energy'],'channels':['n','m','h'],'amplitudes':pack(a),'residue_matrix':pack(np.outer(a.conj(),a)),'moments':mom,'checks':checks,'wall_seconds':time.perf_counter()-t,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024}
 assert checks['total_z_conservation_seed_error']<1e-9 and checks['m_h_opposite_branch_amplitude_max']<1e-8
 (ROOT/f'crossing_{j:.7f}.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
if __name__=='__main__':
 for j in [.3875,.3921875]:run(j)
