"""Analyze representative ground states and both competing branches near each crossing."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json
import cross_scan as cs
from refine_crossings import summarize
from analyze_crossing import analyze
ROOT=cs.ROOT

def ensure_checkpoint(rec):
 file=ROOT/(rec['name']+'.npz')
 if not file.exists():cs.solve_one(rec['Jp'],tuple(rec['cell_n']),rec.get('character',rec.get('angular')),index=rec.get('physical_root_index',0))
 assert file.exists();return rec['name']
def run(interval_file="crossing_intervals.json"):
 selected=set()
 for j in [.2375,.275,.35,.425]:
  proof=summarize(j)
  if proof and proof['global_ground_identified']:selected.add(ensure_checkpoint(proof['winner']))
 if (ROOT/interval_file).exists():
  intervals=json.loads((ROOT/interval_file).read_text())
  for item in intervals['crossings']:
   for j in [item['a'],item['b']]:
    proof=summarize(j);assert proof['global_ground_identified'];selected.add(ensure_checkpoint(proof['winner']))
   # Both eigenstate branches at one common coupling reveal a near-crossing contrast.
   for d in item['competing_branch_records']:
    if d['Jp']==item['a']:
     rec=json.loads((ROOT/(d['name']+'.json')).read_text());selected.add(ensure_checkpoint(rec))
 for name in sorted(selected):analyze(name)
 print(json.dumps({'static_states_completed':sorted(selected)}),flush=True)
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--intervals',default='crossing_intervals.json');a=p.parse_args();run(a.intervals)
