"""Reflection labels and same-symmetry overlaps for tracked physical eigenvectors."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,numpy as np,scipy.linalg as la
from cross_scan import get_model,ROOT
from observables import inversion_vector
from sectors import point_permutations
cachefile=ROOT/'branch_symmetry_and_overlaps.json'
cached={d['name']:d for d in json.loads(cachefile.read_text())} if cachefile.exists() else {}
records=[]
for path in ROOT.glob('p*_e*.json'):
 d=json.loads(path.read_text())
 if not d.get('verified') or not path.with_suffix('.npz').exists():continue
 if d['name'] in cached:
  records.append(cached[d['name']]);continue
 n=tuple(d['cell_n']);r=d['character']
 if not ((n==(0,0) and r in [0,3]) or (n==(1,2) and r==0)):continue
 m,_,_=get_model(d['Jp'],n,r);v=np.load(path.with_suffix('.npz'))['vector'];F=point_permutations()[6];fv=inversion_vector(m,v,F);ov=np.vdot(v,fv);sign=1 if ov.real>=0 else -1;res=la.norm(fv-sign*v)
 rec={'name':d['name'],'Jp':d['Jp'],'cell_n':n,'character':r,'spin':d['spin'],'reflection_overlap_real':float(ov.real),'reflection_eigen_residual':float(res),'reflection_parity':sign if res<1e-6 else None}
 records.append(rec)
# Overlaps compare identical reduced basis conventions within one spatial family.
for rec in records:
 candidates=[p for p in records if p['Jp']<rec['Jp'] and p['cell_n']==rec['cell_n'] and p['character']==rec['character'] and p['spin']==rec['spin'] and p['reflection_parity']==rec['reflection_parity']]
 if candidates:
  jprev=max(p['Jp'] for p in candidates);b=np.load(ROOT/(rec['name']+'.npz'))['vector'];matches=[]
  for prior in candidates:
   if prior['Jp']==jprev:
    a=np.load(ROOT/(prior['name']+'.npz'))['vector'];matches.append((float(abs(np.vdot(a,b))),prior))
  overlap,prior=max(matches,key=lambda x:x[0]);rec['previous_matching_branch']=prior['name'];rec['overlap_with_previous']=overlap;rec['overlap_tracking_note']='Largest overlap among same-quantum-number states at nearest earlier measured coupling; incomplete excited-state coverage.'
(ROOT/'branch_symmetry_and_overlaps.json').write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2))
