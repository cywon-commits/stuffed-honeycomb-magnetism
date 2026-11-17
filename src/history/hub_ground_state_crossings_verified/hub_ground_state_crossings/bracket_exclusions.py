"""Conditional concavity exclusion of third spatial families inside crossing brackets."""
import json
from refine_crossings import ROOT,interval_exclusions

def run():
 file=ROOT/'crossing_intervals.json'
 if not file.exists():return
 out=[]
 for c in json.loads(file.read_text())['crossings']:
  d=interval_exclusions(c['a'],c['b'],None);pairs={(tuple(q[0]),q[1]) for q in [c['left_quantum_numbers'],c['right_quantum_numbers']]};other=[x for x in d['family_margins'] if (tuple(x['cell_n']),x['character']) not in pairs];out.append({'a':c['a'],'b':c['b'],'other_family_margins':other,'all_third_spatial_families_excluded':all(x['minimum_margin']>1e-8 for x in other),'scope':'Conditional on verified numerical lowest-root identification; no exclusion of unsampled spin changes within the two competing spatial families.'})
 (ROOT/'bracket_third_family_exclusions.json').write_text(json.dumps(out,indent=2))
if __name__=='__main__':run()
