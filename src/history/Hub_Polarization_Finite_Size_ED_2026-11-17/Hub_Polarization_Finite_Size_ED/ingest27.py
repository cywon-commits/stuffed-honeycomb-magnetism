import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,pathlib,csv,shutil,hashlib,numpy as np
from geometry import ROOT,make
from observables import values
OLD=ROOT.parent;D=ROOT/'reference27';D.mkdir(exist_ok=True)
def run():
 if (D/'input_records.json').exists():
  g=make([[3,0],[0,3]]);records=[]
  for rec in json.loads((D/'input_records.json').read_text()):
   source=D/rec['source_file'];assert hashlib.sha256(source.read_bytes()).hexdigest()==rec['source_sha256']
   G=np.load(D/f"correlations_r{rec['r']:.9f}.npz")['Gdot'];o=values(g,G,rec['S']);err=abs(sum((1 if t==0 else rec['r'])*G[a,b] for a,b,t in g['bonds'])-rec['energy']);assert err<1e-7 and o['spin_sum_rule_error']<1e-7;records.append({**rec,'observables':o})
  (ROOT/'reference27_results.json').write_text(json.dumps(records,indent=2));print('PASS cached verified 27-spin sources, energy and spin-sum checks.',flush=True);return
 g=make([[3,0],[0,3]]);perm=[3*((int(x)%3)*3+(int(y)%3))+s for x,y in g['cells'] for s in range(3)];records=[]
 def add(r,E,S,k,hres,sres,G,path,meta,source):
  G=G[np.ix_(perm,perm)];o=values(g,G,S);err=abs(sum((1 if t==0 else r)*G[a,b] for a,b,t in g['bonds'])-E);assert err<1e-7 and o['spin_sum_rule_error']<1e-7,(source,err,o['spin_sum_rule_error']);target=D/(path.name);shutil.copy2(path,target);shutil.copy2(meta,D/meta.name);np.savez_compressed(D/f'correlations_r{r:.9f}.npz',Gdot=G,r=r,T=g['T']);records.append({'shape':'N27','N':27,'Ncell':9,'r':r,'energy':E,'S':S,'k1':k[0],'k2':k[1],'H_residual':hres,'S2_residual':sres,'energy_from_G_error':err,'provenance':source,'source_file':path.name,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'observables':o})
 for r in [.2,.5]:
  meta=OLD/'static'/f'static_j{r:g}.json';m=json.loads(meta.read_text())['summary'];p=OLD/'static'/f'static_j{r:g}_realspace.npz';G=np.load(p)['G_selected'];add(r,m['energy'],m['S_total'],m['cell_ground_momentum'],m['energy_residual'],m['S2_eigen_residual'],G,p,meta,'Re-ingested verified 27-spin checkpoint-derived real-space matrix')
 for meta in sorted((OLD/'crossing').glob('obs_*_e0.json')):
  m=json.loads(meta.read_text())['summary']
  if m['energy_above_global_ground']>1e-8 or m['role']!='GROUND':continue
  r=m['Jp'];p=meta.with_name(meta.stem+'_realspace.npz');z=np.load(p);G=z['Gdot'];level=next(x for x in csv.DictReader((OLD/'crossing/all_levels.csv').open()) if x['name']==m['name']);add(r,m['energy'],m['spin'],(np.array(m['cell_n'])/3).tolist(),m['residual'],float(level['S2_residual']),G,p,meta,'Verified crossing sweep: ground-state role and energy checked against archived levels')
 for r in [.6,.7]:
  p=OLD/'triangle'/f'raw_observables_j{r:g}.npz';meta=OLD/'triangle'/f'observables_j{r:g}.json';m=json.loads(meta.read_text());add(r,m['energy'],m['validation']['spin']['nearest_half_integer_spin'],m['cell_momentum'],m['validation']['H_residual'],m['validation']['spin']['S2_eigen_residual'],np.load(p)['Gdot'],p,meta,'Verified additional ground-sector scan, static correlation matrix')
 # Preserve exact measured 27-spin brackets, not prose estimates.
 shutil.copy2(OLD/'crossing/crossing_brackets.csv',D/'crossing_brackets_27.csv');shutil.copy2(OLD/'crossing/ground_samples.csv',D/'ground_samples_27.csv');shutil.copy2(OLD/'triangle/triangle_sector_comparison.csv',D/'later_sector_coverage_27.csv')
 audits=[]
 for r,cp,n,a in [(.2,OLD/'static/j0.2_n12_r0.npz',[1,2],0),(.35,OLD/'crossing/p0.350000000_n00_r0_e0.npz',[0,0],0),(.5,OLD/'static/j0.5_n00_r3.npz',[0,0],3),(.6,OLD/'triangle/p0.600000000_n00_r0_e0.npz',[0,0],0),(.7,OLD/'triangle/p0.700000000_n00_r3_e0.npz',[0,0],3)]:
  shutil.copy2(cp,D/'checkpoints'/cp.name);audits.append({'r':r,'file':cp.name,'cell_n':n,'angular':a,'sha256':hashlib.sha256(cp.read_bytes()).hexdigest()})
 (D/'audit_inputs.json').write_text(json.dumps(audits,indent=2));records=sorted({round(x['r'],9):x for x in records}.values(),key=lambda x:x['r']);(ROOT/'reference27_results.json').write_text(json.dumps(records,indent=2));print('PASS',len(records),'verified 27-spin ground correlation matrices; energy and total-spin sum rules reproduced.',flush=True)
if __name__=='__main__':run()
