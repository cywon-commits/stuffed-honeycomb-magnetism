"""Adaptive, cold-start global ground-sector refinement with concavity exclusions."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,time,numpy as np
import cross_scan as cs
ROOT=cs.ROOT
FAMILIES=cs.PLAN

def root_record(j,n,r,index=0):
 path=ROOT/(cs.tag(j,n,r,index)+'.json')
 if path.exists():
  try:d=json.loads(path.read_text())
  except json.JSONDecodeError:return None
  return d if d.get('verified') else None
 if j in [.2,.5] and index==0:
  d=cs.old_record(j,n,r);p=json.loads((cs.OLD/f'j{j:g}_n{n[0]}{n[1]}_r{r}.json').read_text());profiles=json.loads((ROOT/'endpoint_spin_profiles.json').read_text());spin=next(x for x in profiles if x['Jp']==j and x['cell_n']==list(n) and x['character']==r)
  return {'name':f'endpoint_j{j:g}_n{n[0]}{n[1]}_r{r}','Jp':j,'cell_n':list(n),'character':r,'physical_root_index':0,'energy':p['energy'],'residual':p['residual'],'spin':spin['nearest_half_integer_spin'],'spin_check':spin,'verified':True,'dimension':p['dimension'],'reference_endpoint':True}
 return None

def all_j():
 vals={.2,.5}
 for p in ROOT.glob('p*_e0.json'):
  try:d=json.loads(p.read_text())
  except json.JSONDecodeError:continue
  if d.get('verified'):vals.add(d['Jp'])
 return sorted(vals)
def add_parity(d):
 cachefile=ROOT/'reflection_labels.json'
 cache=json.loads(cachefile.read_text()) if cachefile.exists() else {}
 if d['name'] in cache:d.update(cache[d['name']]);return d
 n=tuple(d['cell_n']);r=d['character'];parity=None;res=None
 if n==(0,1):parity=1 if r==0 else -1
 elif (n==(0,0) and r in [0,3]) or (n==(1,2) and r==0):
  from observables import inversion_vector
  from sectors import point_permutations
  import scipy.linalg as la
  if d.get('reference_endpoint'):checkpoint=cs.OLD/f"j{d['Jp']:g}_n{n[0]}{n[1]}_r{r}.npz"
  else:checkpoint=ROOT/(d['name']+'.npz')
  if checkpoint.exists():
   m,_,_=cs.get_model(d['Jp'],n,r);v=np.load(checkpoint)['vector'];F=point_permutations()[6];fv=inversion_vector(m,v,F);ov=np.vdot(v,fv);sign=1 if ov.real>=0 else -1;res=float(la.norm(fv-sign*v));parity=sign if res<1e-6 else None
 info={'reflection_parity':parity,'reflection_eigen_residual':res};cache[d['name']]=info;cachefile.write_text(json.dumps(cache,indent=2));d.update(info);return d

def label(d):
 d=add_parity(d);return (tuple(d['cell_n']),d['character'],d['spin'],d.get('reflection_parity'))

def summarize(j,save=True):
 entries=[];actual=[]
 for n,r in FAMILIES:
  d=root_record(j,n,r)
  if d:entries.append({'family':[list(n),r],'kind':'MEASURED_EIGENPAIR','record':d});actual.append(d)
  else:
   bound,ends=cs.lower_chord(j,n,r);entries.append({'family':[list(n),r],'kind':'CONCAVITY_LOWER_BOUND','energy_lower_bound':bound,'bounds_endpoint_couplings':ends})
 if not actual:return None
 best=add_parity(min(actual,key=lambda d:d['energy']));unresolved=[x for x in entries if x['kind']=='CONCAVITY_LOWER_BOUND' and x['energy_lower_bound']<=best['energy']+1e-8]
 sectors={}
 for n in [(0,0),(0,1),(1,2)]:
  measured=[d for d in actual if tuple(d['cell_n'])==n];sbest=min(measured,key=lambda d:d['energy']) if measured else None
  missing=[x for x in entries if tuple(x['family'][0])==n and x['kind']=='CONCAVITY_LOWER_BOUND'];lb=min([d['energy']-max(1e-9,5*d['residual']) for d in measured]+[x['energy_lower_bound'] for x in missing])
  sectors[str(n)]={'cell_n':n,'lowest_measured':sbest,'energy_lower_bound':lb,'minimum_identified':bool(sbest is not None and all(x['energy_lower_bound']>sbest['energy']+1e-8 for x in missing))}
 out={'Jp':j,'global_ground_identified':not unresolved,'winner':best,'label':[list(best['cell_n']),best['character'],best['spin'],best.get('reflection_parity')],'sector_minima':sectors,'family_evidence':entries,'unresolved_families':unresolved,'all_sector_minima_identified':all(x['minimum_identified'] for x in sectors.values()),'bound_scope':'Conditional numerical exclusions: lowest-root identification at measured endpoints is supported by independent cold starts, direct residuals and validated ED, not an interval-arithmetic proof.'}
 if save:(ROOT/f'point_{j:.9f}.json').write_text(json.dumps(out,indent=2))
 return out

def search_point(j,full_sector_minima=False):
 j=round(float(j),9)
 # Cold Gamma solves also audit for spin changes inside the same spatial sector.
 for n,r in [((0,0),0),((0,0),3)]:cs.solve_one(j,n,r,ncv=64)
 for step in range(len(FAMILIES)):
  d=summarize(j);unresolved=d['unresolved_families']
  if full_sector_minima:
   for x in d['family_evidence']:
    if x['kind']=='CONCAVITY_LOWER_BOUND':
     sn=d['sector_minima'][str(tuple(x['family'][0]))];upper=sn['lowest_measured']['energy'] if sn['lowest_measured'] else np.inf
     if x['energy_lower_bound']<=upper+1e-8 and x not in unresolved:unresolved.append(x)
  if not unresolved:break
  e=min(unresolved,key=lambda x:x['energy_lower_bound']);n,r=e['family'];cs.solve_one(j,tuple(n),r,ncv=64)
 d=summarize(j);assert d['global_ground_identified'],d['unresolved_families'];print(json.dumps({'point_verified':j,'energy':d['winner']['energy'],'label':d['label']}),flush=True);return d

def upper_lines():
 lines=[]
 for p in ROOT.glob('p*_e*.json'):
  try:d=json.loads(p.read_text())
  except json.JSONDecodeError:continue
  if d.get('verified') and 'slope_dE_dJp' in d:lines.append((d['energy']-d['Jp']*d['slope_dE_dJp']+1e-9,d['slope_dE_dJp']))
 for j in [.2,.5]:
  d=json.loads((cs.OLD/f'static_j{j:g}.json').read_text())['summary'];slope=d['hub_rim_energy_from_Fourier']/j;lines.append((d['energy']-j*slope+1e-9,slope))
 return lines

def interval_exclusions(a,b,winner_family):
 lines=upper_lines();grid={a,b}
 # All affine-upper-envelope breakpoints, plus every lower-chord knot, suffice.
 for i,(c,s) in enumerate(lines):
  for d,t in lines[:i]:
   if abs(s-t)>1e-12:
    x=(d-c)/(s-t)
    if a<x<b:grid.add(float(x))
 for n,r in FAMILIES:
  for d in cs.points(n,r):
   if a<d['Jp']<b:grid.add(d['Jp'])
 grid=np.array(sorted(grid));upper=np.min(np.array([c+s*grid for c,s in lines]),axis=0);unresolved=[];details=[]
 for n,r in FAMILIES:
  if (tuple(n),r)==winner_family:continue
  pts=cs.points(n,r);xs=np.array([p['Jp'] for p in pts]);ys=np.array([p['energy']-max(1e-9,5*p['residual']) for p in pts]);margin=float(np.min(np.interp(grid,xs,ys)-upper));details.append({'cell_n':n,'character':r,'minimum_margin':margin})
  if margin<=1e-8:unresolved.append((n,r))
 return {'a':a,'b':b,'excluded_other_families':not unresolved,'unresolved_families':unresolved,'family_margins':details,'sampling_points_for_piecewise_linear_bounds':len(grid)}

def ensure_spin_branch(j,n,r,spin,parity=None,max_roots=4):
 for idx in range(max_roots):
  rec=root_record(j,n,r,idx)
  if rec is None or (rec.get('reference_endpoint') and (rec['spin']!=spin or (parity is not None and add_parity(rec).get('reflection_parity')!=parity))):rec=cs.solve_one(j,n,r,index=idx,ncv=64)
  if not rec.get('verified'):return None
  if rec.get('verified') and rec['spin']==spin and (parity is None or add_parity(rec).get('reflection_parity')==parity):return rec
 return None

def main(tol=.004999,max_new_points=32,strict_exclusions=False):
 for j in [.275,.35,.425]:
  while not (ROOT/f'initial_j{j:g}.json').exists():time.sleep(20)
  assert summarize(j)['global_ground_identified']
 audited=set();new=0
 while True:
  summaries=[summarize(j) for j in all_j()];summaries=[d for d in summaries if d and d['global_ground_identified']];summaries.sort(key=lambda d:d['Jp']);jobs=[];certs=[]
  for left,right in zip(summaries[:-1],summaries[1:]):
   a,b=left['Jp'],right['Jp'];different=label(left['winner'])!=label(right['winner'])
   if different:
    if b-a>tol:jobs.append((0,b-a,a,b,'different_quantum_numbers'))
   else:
    family=(tuple(left['winner']['cell_n']),left['winner']['character']);cert=interval_exclusions(a,b,family);certs.append(cert)
    if strict_exclusions and not cert['excluded_other_families'] and b-a>tol:jobs.append((1,b-a,a,b,'unexcluded_alternative_spatial_family'))
    elif b-a>.03751 and (a,b) not in audited:jobs.append((2,b-a,a,b,'same_family_spin_midpoint_audit'))
  (ROOT/'interval_exclusions.json').write_text(json.dumps(certs,indent=2))
  if not jobs or new>=max_new_points:break
  jobs.sort(key=lambda x:(x[0],-x[1]));_,width,a,b,why=jobs[0];mid=round((a+b)/2,9);audited.add((a,b));print(json.dumps({'refining':[a,b],'midpoint':mid,'reason':why}),flush=True);search_point(mid);new+=1
 summaries=[summarize(j) for j in all_j()];summaries=sorted([d for d in summaries if d and d['global_ground_identified']],key=lambda d:d['Jp']);crossings=[]
 for left,right in zip(summaries[:-1],summaries[1:]):
  A,B=left['winner'],right['winner'];a,b=left['Jp'],right['Jp']
  if label(A)==label(B):continue
  item={'a':a,'b':b,'width':b-a,'left_quantum_numbers':[A['cell_n'],A['character'],A['spin'],A.get('reflection_parity')],'right_quantum_numbers':[B['cell_n'],B['character'],B['spin'],B.get('reflection_parity')],'target_width_met':b-a<=tol}
  # Actual competing fixed-spin eigenpairs on both sides, not linear extrapolated Ritz values.
  branches=[]
  for n,r,s,parity in [label(A),label(B)]:
   for j in [a,b]:
    rec=root_record(j,n,r)
    if rec is None or rec['spin']!=s or (parity is not None and add_parity(rec).get('reflection_parity')!=parity):rec=ensure_spin_branch(j,n,r,s,parity)
    if rec:branches.append({'Jp':j,'cell_n':n,'character':r,'spin':s,'reflection_parity':parity,'name':rec['name'],'energy':rec['energy'],'residual':rec['residual'],'slope_dE_dJp':rec.get('slope_dE_dJp')})
  item['competing_branch_records']=branches
  for j,key in [(a,'gap_right_minus_left_at_a'),(b,'gap_right_minus_left_at_b')]:
   av=next((d['energy'] for d in branches if d['Jp']==j and (tuple(d['cell_n']),d['character'],d['spin'],d.get('reflection_parity'))==label(A)),None);bv=next((d['energy'] for d in branches if d['Jp']==j and (tuple(d['cell_n']),d['character'],d['spin'],d.get('reflection_parity'))==label(B)),None)
   if av is not None and bv is not None:item[key]=bv-av
  crossings.append(item)
 out={'target_width':tol,'new_point_count_this_process':new,'refinement_coupling_count':len([d for d in summaries if d['Jp'] not in [.2,.275,.35,.425,.5]]),'strict_exclusions_requested':strict_exclusions,'new_point_limit_reached':new>=max_new_points,'crossings':crossings,'verified_couplings':[d['Jp'] for d in summaries],'all_bracket_widths_met':all(d['target_width_met'] for d in crossings),'limit_note':'Spatial-family concavity exclusions are recorded continuously between samples. Same-spatial-family spin changes or islands narrower than the sampled resolution are not rigorously excluded.'};(ROOT/'crossing_intervals.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--strict-exclusions',action='store_true');p.add_argument('--max-new-points',type=int,default=32);p.add_argument('--target-width',type=float,default=.004999);a=p.parse_args();main(a.target_width,a.max_new_points,a.strict_exclusions)
