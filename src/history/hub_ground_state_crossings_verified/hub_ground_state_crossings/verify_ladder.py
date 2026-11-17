"""Static global-spin ladder check; no frequency-dependent or spectral calculation."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,numpy as np,scipy.linalg as la,argparse,time,resource
from cyclic_model import CyclicModel
from ed import Model
from observables import convert_basis,spin_seed
from refine_crossings import summarize
ROOT=pathlib.Path(__file__).resolve().parent

def verify(name):
 rec=json.loads((ROOT/(name+'.json')).read_text());assert rec['verified'];proof=summarize(rec['Jp']);assert proof['global_ground_identified'] and proof['winner']['name']==name
 z=np.load(ROOT/(name+'.npz'));j=rec['Jp'];S=rec['spin'];E=rec['energy'];src=CyclicModel(jp=j,cell_n=tuple(z['cell_n']),angular=int(z['angular']));m=Model(up=14,primitive=False,jp=j);m.sector(z['q']);v=convert_basis(src,m,z['vector']);del src;start=time.perf_counter();checks=[]
 for up in range(15,int(round(13.5+S+1))+1):
  target=Model(up=up,primitive=False,jp=j);target.sector(z['q']);w=3*spin_seed(m,target,v,[0,0],'total','+');M=m.up-13.5;norm=float(np.vdot(w,w).real);expected=S*(S+1)-M*(M+1);out={'source_M':M,'target_M':M+1,'norm_squared':norm,'expected':expected,'target_dimension':target.nb}
  if expected>1e-10:
   assert abs(norm-expected)<1e-8;w/=np.sqrt(norm);hv=target.mv(w);res=float(la.norm(hv-E*w));assert res<1e-8;out['energy_residual_at_original_E']=res
  else:assert norm<1e-18
  checks.append(out);m=target;v=w
 result={'name':name,'Jp':j,'spin':S,'ladder_checks':checks,'wall_seconds':time.perf_counter()-start,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'scope':'Static SU(2) identities only; all ground searches remain in Nup=14.'};(ROOT/('ladder_'+name+'.json')).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('name');a=p.parse_args();verify(a.name)
