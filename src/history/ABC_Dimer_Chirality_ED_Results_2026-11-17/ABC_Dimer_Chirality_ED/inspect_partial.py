import json,numpy as np,pathlib
R=pathlib.Path(__file__).resolve().parent
for f in sorted(R.glob('*_spectrum.json')):
 d=json.loads(f.read_text());p=[x for x in d['poles'] if x['verified_eigenpair'] and x['scalar_accessible'] and x['total_seed_weight']>1e-12]
 print(f.name,'verified scalar',len(p),'minomega',min([x['omega'] for x in p],default=None),'maxres',max([x['H_residual'] for x in p],default=0))
