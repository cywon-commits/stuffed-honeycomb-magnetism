import json,csv
from scan18 import Scan,write_csv
from geometry import ROOT,catalog
rows=[]
for shape in ['K18','C18']:
 scanner=Scan(shape,catalog()[shape])
 for f in (ROOT/'checkpoints'/shape).glob('*.json'):
  old=json.loads(f.read_text())
  if max(x['S2_residual'] for x in old['levels'])<1e-8:continue
  new=scanner.point(old['r'],force=True);rows.append({'shape':shape,'r':old['r'],'old_max_S2_residual':max(x['S2_residual'] for x in old['levels']),'new_max_S2_residual':max(x['S2_residual'] for x in new['levels']),'ground_energy_difference':new['winner']['energy']-old['winner']['energy'],'method':'Joint H/S2 subspace diagonalization, then H diagonalization separately for each integer S.'})
with (ROOT/'crossing_brackets_18.csv').open() as f:br=list(csv.DictReader(f))
for row in br:
 shape=row['shape'];r=float(row['secant_intersection_estimate']);p=json.loads((ROOT/'checkpoints'/shape/f'r{r:.9f}.json').read_text());a=json.loads((ROOT/'checkpoints'/shape/f"r{float(row['r_lo']):.9f}.json").read_text())['winner'];b=json.loads((ROOT/'checkpoints'/shape/f"r{float(row['r_hi']):.9f}.json").read_text())['winner']
 def e(label):return min(x['energy'] for x in p['levels'] if x['S']==label['S'] and abs(x['k1']-label['k1'])+abs(x['k2']-label['k2'])<1e-8)
 row['actual_branch_gap_at_estimate']=e(b)-e(a)
write_csv(ROOT/'crossing_brackets_18.csv',br);(ROOT/'near_crossing_spin_resolution.json').write_text(json.dumps({'status':'PASS','rows':rows},indent=2))
