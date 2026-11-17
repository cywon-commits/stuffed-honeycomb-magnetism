import json
from experiment import ROOT,dump
old=[]
for line in (ROOT/'j1_G_chi_solve.log').read_text().splitlines():
 try:r=json.loads(line)
 except json.JSONDecodeError:continue
 if 'index' in r and r['scalar_accessible']:old.append(r)
new=json.loads((ROOT/'j1_G_a0_p-1_spectrum.json').read_text())['poles'];comparisons=[]
for a in old:
 b=min(new,key=lambda p:abs(p['energy']-a['energy']));comparisons.append({'omega':a['omega'],'old_requested_pairs':8,'expanded_requested_pairs':len(new),'energy_difference':abs(a['energy']-b['energy']),'total_weight_difference':abs(a['total_seed_weight']-b['total_seed_weight']),'old_direct_residual':a['H_residual'],'expanded_direct_residual':b['H_residual']})
assert max(x['energy_difference'] for x in comparisons)<1e-9 and max(x['total_weight_difference'] for x in comparisons)<1e-8
dump(ROOT/'convergence_stability.json',{'status':'PASS','comparisons':comparisons});print('PASS stable verified Γ chirality poles in expanded solve.')
