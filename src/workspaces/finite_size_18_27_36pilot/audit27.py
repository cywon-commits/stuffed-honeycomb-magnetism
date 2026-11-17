import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,sys,json,time,resource,numpy as np,scipy.linalg as la
ROOT=pathlib.Path(__file__).resolve().parent;D=ROOT/'reference27';sys.path.insert(0,str(D/'code'))
from cyclic_model import CyclicModel
from observables import total_spin_check

def subsystem(m,v,ix):
 old=(m.bonds,m.flips,m.js,m.diag)
 try:
  pairs=[(int(a),int(b),2.) for i,a in enumerate(ix) for b in ix[i+1:]];m.bonds=pairs;m.flips=np.array([(1<<a)|(1<<b) for a,b,j in pairs],np.uint32);m.js=np.full(len(pairs),2.);m.diag=np.full(m.nb,.75*len(ix))
  for a,b,j in pairs:m.diag+=2*(.25-.5*(((m.masks>>a)^(m.masks>>b))&1))
  out=m.mv(v);mean=np.vdot(v,out).real;var=la.norm(out)**2-mean**2;return float(mean),float(max(0,var))
 finally:m.bonds,m.flips,m.js,m.diag=old
results=[];start=time.perf_counter();old=json.loads((ROOT/'reference27_results.json').read_text())
for rec in json.loads((D/'audit_inputs.json').read_text()):
 m=CyclicModel(jp=rec['r'],cell_n=tuple(rec['cell_n']),angular=rec['angular']);z=np.load(D/'checkpoints'/rec['file']);v=z['vector'];v/=la.norm(v);e=float(z['energy']);res=float(la.norm(m.mv(v)-e*v));spin=total_spin_check(m,v);H,Hvar=subsystem(m,v,list(range(2,27,3)));R,Rvar=subsystem(m,v,[i for i in range(27) if i%3!=2]);reference=next(x for x in old if abs(x['r']-rec['r'])<1e-9);herr=abs(H-reference['observables']['H2']);rerr=abs(R-reference['observables']['R2']);assert res<1e-8 and spin['S2_eigen_residual']<1e-7 and max(herr,rerr)<1e-7
 results.append({**rec,'energy':e,'fresh_H_residual':res,'fresh_spin':spin,'H2_operator':H,'R2_operator':R,'H2_variance':Hvar,'R2_variance':Rvar,'H2_vs_ingested_G_error':herr,'R2_vs_ingested_G_error':rerr});print(results[-1],flush=True);del m
(ROOT/'reference27_checkpoint_audit.json').write_text(json.dumps({'status':'PASS','rows':results,'wall_seconds':time.perf_counter()-start,'peak_RSS_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024},indent=2))
import subprocess
subprocess.run([sys.executable,str(ROOT/'preanalysis.py')],check=True)
