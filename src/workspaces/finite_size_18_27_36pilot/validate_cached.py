import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,csv,numpy as np,scipy.linalg as la
from geometry import ROOT,make
from ed18 import Basis
checks=[]
for shape in ['K18','C18']:
 files=sorted((ROOT/'checkpoints'/shape).glob('*.json'));z0=np.load(files[0].with_suffix('.npz'));g=make(z0['T']);m=Basis(g);dims=[]
 for q in g['momenta']:m.sector(q);dims.append(m.nb)
 assert len(dims)==6 and sum(dims)==len(m.states)==48620
 for f in files:
  p=json.loads(f.read_text());z=np.load(f.with_suffix('.npz'));m.sector(z['q']);m.set_r(float(z['r']));v=z['vector'];h=m.mv(v);res=float(la.norm(h-float(z['energy'])*v));S,S2,sres=m.spin(v);norm=float(abs(la.norm(m.expand(v))-1));G=z['Gdot'];energy=sum((1 if t==0 else float(z['r']))*G[a,b] for a,b,t in g['bonds']);assert res<1e-8 and sres<1e-8 and S==p['winner']['S'] and norm<1e-12 and abs(energy-float(z['energy']))<1e-9
  checks.append({'shape':shape,'r':float(z['r']),'fresh_checkpoint_H_residual':res,'fresh_checkpoint_S2_residual':sres,'reconstructed_norm_error':norm,'energy_from_saved_G_error':abs(energy-float(z['energy']))})
for p in csv.DictReader((ROOT/'crossing_brackets_18.csv').open()):
 assert float(p['width'])<=.005 and float(p['delta_E_hi_minus_lo_branch_at_rlo'])>100*float(p['maximum_endpoint_H_residual']) and float(p['delta_E_hi_minus_lo_branch_at_rhi'])<-100*float(p['maximum_endpoint_H_residual'])
from scan18 import write_csv
write_csv(ROOT/'fresh_checkpoint_verification.csv',checks);(ROOT/'cached_reproduction_validation.json').write_text(json.dumps({'status':'PASS','checkpoints_freshly_reconstructed_and_verified':len(checks),'all_momentum_dimensions_sum_to_full_Sz0_48620':True,'max_fresh_H_residual':max(x['fresh_checkpoint_H_residual'] for x in checks),'max_fresh_S2_residual':max(x['fresh_checkpoint_S2_residual'] for x in checks)},indent=2));print('PASS',len(checks),'fresh checkpoint checks and sector completeness.')
