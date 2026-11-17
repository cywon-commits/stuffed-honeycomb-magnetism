import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,numpy as np,scipy.linalg as la
from experiment import ROOT,geometry_weights
from group_poles import unpack
from ensemble import unfold,fourier
from build_data import csvwrite
R=json.loads((ROOT/'responses.json').read_text());ref=np.load(ROOT/'inputs/raw_observables_j1.npz');rows=[];cells=[];vectors={}
def entry(name,kind,mode,u):
 r=R[name];M=r.get('moments',r);m=[float(np.vdot(u,unpack(M[f'M{p}'])@u).real) for p in range(3)];w=sum(float(np.vdot(u,unpack(p['residue'])@u).real) for p in r['poles']);return {'response':name,'operator':kind,'mode':mode,'M0_positive_energy':m[0],'M1':m[1],'M2':m[2],'zero_energy_weight':float(np.vdot(u,unpack(M['M0_elastic'])@u).real),'verified_weight':w,'verified_fraction':w/m[0] if m[0]>1e-12 else None,'unidentified_weight':m[0]-w}
for name,Q in [('isotropic_G_pair_average',[0,0]),('isotropic_Q_pair_average',[8/9,7/9]),('isotropic_K_pair_average',[1/3,2/3]),('j09_G',[0,0]),('j09_q',[8/9,7/9]),('j07_G',[0,0])]:
 if name not in R:continue
 P=geometry_weights(np.array(Q)) if len(R[name]['channel_names'])==18 else np.eye(6)
 for kind,a,off in [('D',0,0),('chi',4,108)]:
  data=unfold(fourier(ref['ensemble_full_connected'][off:off+54,off:off+54]),kind)[0];refQ=min(data,key=lambda x:la.norm(x['Q']-Q))
  for mode in range(2):
   u=np.zeros(6,complex);u[a:a+2]=refQ['eigenvectors'][:,mode];u=P@u;vectors[name+f'_{kind}_{mode}']=u;rows.append(entry(name,kind,mode,u))
for name in ['j09_G','j09_q','j07_G']:
 if name not in R:continue
 M=R[name]['moments'];X=unpack(M['M0'])
 for a,kind in enumerate(['D','Ddag','chi']):
  lo=6*a;e,U=la.eigh((X[lo:lo+6,lo:lo+6]+X[lo:lo+6,lo:lo+6].conj().T)/2)
  for mode in range(6):
   u=np.zeros(18,complex);u[lo:lo+6]=U[:,mode];vectors[name+f'_cell_{kind}_{mode}']=u;cells.append(entry(name,kind,mode,u))
csvwrite('face_mode_coverage.csv',rows);csvwrite('cell_mode_coverage.csv',cells);np.savez_compressed(ROOT/'face_mode_vectors.npz',**vectors);print('Face-mode coverage written.')
