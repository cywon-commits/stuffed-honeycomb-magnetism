"""Tighten measured brackets using the already solved secant diagnostic points."""
import json,csv
from geometry import ROOT
from scan18 import write_csv
with (ROOT/'crossing_brackets_18.csv').open() as f:rows=list(csv.DictReader(f))
for row in rows:
 if row.get('tightened_using_measured_intersection_point')=='True':continue
 a,b=float(row['r_lo']),float(row['r_hi']);m=float(row['secant_intersection_estimate']);d=float(row['actual_branch_gap_at_estimate']);p=json.loads((ROOT/'checkpoints'/row['shape']/f'r{m:.9f}.json').read_text());res=p['winner']['H_residual'];assert abs(d)>100*res
 row['original_coarse_r_lo']=a;row['original_coarse_r_hi']=b
 if d>0:row['r_lo']=m;row['delta_E_hi_minus_lo_branch_at_rlo']=d
 else:row['r_hi']=m;row['delta_E_hi_minus_lo_branch_at_rhi']=d
 row['width']=float(row['r_hi'])-float(row['r_lo']);row['tightened_using_measured_intersection_point']=True
 for key in ['k_lo','k_hi']:
  q=json.loads(row[key]);row[key]=json.dumps([x[-2:] for x in q])
write_csv(ROOT/'crossing_brackets_18.csv',rows)
# Re-evaluate the residual envelope for the final, newly tightened endpoints.
with (ROOT/'crossing_brackets_18.csv').open() as f:rows=list(csv.DictReader(f))
for row in rows:
 residuals=[]
 for name in ['r_lo','r_hi']:
  p=json.loads((ROOT/'checkpoints'/row['shape']/f"r{float(row[name]):.9f}.json").read_text())
  for skey,qkey in [('S_lo','k_lo'),('S_hi','k_hi')]:
   S=int(float(row[skey]));qs=json.loads(row[qkey]);cand=[x for x in p['levels'] if x['S']==S and any(abs(x['k1']-q[0])+abs(x['k2']-q[1])<1e-8 for q in qs)];assert cand;residuals.append(min(cand,key=lambda x:x['energy'])['H_residual'])
 row['maximum_endpoint_H_residual']=max(residuals)
write_csv(ROOT/'crossing_brackets_18.csv',rows)
