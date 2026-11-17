import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,scipy.linalg as la
from experiment import ROOT,geometry_weights,pack,dump
from ensemble import unfold,fourier
z=np.load(ROOT/'fresh_local_sum_rule_data.npz');A=z['ground_cell_amplitudes'];G=A.conj().T@A;perm=np.r_[np.arange(6,12),np.arange(6),np.arange(12,18)];sign=np.r_[np.ones(12),-np.ones(6)];Gp=G[np.ix_(perm,perm)].conj()*sign[:,None]*sign[None,:];old=np.load(ROOT/'inputs/isotropic_ground_manifold.npz')['known_ground_connected'];records=[]
for Q in [[0,0],[1/3,2/3],[2/3,1/3]]:
 W=geometry_weights(np.array(Q));selected=W.conj().T@G@W;partner=W.conj().T@Gp@W;avg=(selected+partner)/2
 for kind,a,off in [('D',0,0),('chi',4,108)]:
  expected=unfold(fourier(old[off:off+54,off:off+54]),kind)[0];d=min(expected,key=lambda x:la.norm(x['Q']-Q));error=float(la.norm(avg[a:a+2,a:a+2]-d['matrix']));assert error<1e-7;records.append({'physical_Q':Q,'kind':kind,'selected_K_elastic':pack(selected[a:a+2,a:a+2]),'selected_Kprime_elastic':pack(partner[a:a+2,a:a+2]),'orbital_average_elastic':pack(avg[a:a+2,a:a+2]),'orbital_average_eigenvalues':la.eigvalsh(avg[a:a+2,a:a+2]).tolist(),'reference_error':error})
dump(ROOT/'ground_projection_mode_checks.json',{'status':'PASS','records':records,'source':'Fresh cell seeds and RDMs using refined K/Kprime ground pair; no finite-energy eigensolver or quadrature used.'});print('PASS Γ/K/Kprime ground projections; max error',max(r['reference_error'] for r in records))
