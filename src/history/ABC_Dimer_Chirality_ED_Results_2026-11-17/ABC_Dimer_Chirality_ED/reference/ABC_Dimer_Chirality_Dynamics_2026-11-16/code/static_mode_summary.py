#!/usr/bin/env python3
"""Recompute 27spin isotropic M0 and known ground-manifold deflation by physical Q.
Uses ONLY supplied static Gram matrices, not a dynamical ED calculation.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
import pathlib,sys,csv,numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]
SRC=ROOT/'reference/ABC_triangle_ED/triangle'
sys.path.insert(0,str(SRC))
from ensemble import fourier,unfold
z=np.load(SRC/'raw_observables_j1.npz')
g=np.load(SRC/'isotropic_ground_manifold.npz')
rows=[]
for kind,ofs in [('D',0),('chi',108)]:
 slices={
   'total':z['ensemble_full_connected'][ofs:ofs+54,ofs:ofs+54],
   'ground':g['known_ground_connected'][ofs:ofs+54,ofs:ofs+54],
   'excited':g['excited_connected'][ofs:ofs+54,ofs:ofs+54],
 }
 maps={label:unfold(fourier(C),kind)[0] for label,C in slices.items()}
 for i in range(27):
  t=maps['total'][i]
  E,U=np.linalg.eigh(t['matrix'])
  for index in range(2):
   u=U[:,index]
   strengths={label:float(np.vdot(u,maps[label][i]['matrix']@u).real) for label in ['total','ground','excited']}
   assert abs(strengths['total']-strengths['ground']-strengths['excited'])<3e-8
   row={'operator':kind,'Q1':float(t['Q'][0]),'Q2':float(t['Q'][1]),
        'cell_n1':int(t['cell_n'][0]),'cell_n2':int(t['cell_n'][1]),
        'fold':int(t['fold']),'eigenmode':index,
        'M0_total':strengths['total'],'M0_ground_pair':strengths['ground'],
        'M0_positive_energy':strengths['excited']}
   rows.append(row)
with (ROOT/'results/27spin_static_mode_targets.csv').open('w',newline='') as f:
 wr=csv.DictWriter(f,fieldnames=list(rows[0]));wr.writeheader();wr.writerows(rows)
for kind in ['D','chi']:
 top=sorted([r for r in rows if r['operator']==kind],key=lambda x:x['M0_positive_energy'],reverse=True)[:4]
 for r in top: print(kind,r)
print('PASS: 27 momentum * 2 modes * 2 operators, zero-energy projector subtraction validated')
