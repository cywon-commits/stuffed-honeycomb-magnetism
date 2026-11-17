import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,numpy as np,scipy.linalg as la
from experiment import ROOT,pack,dump
from group_poles import unpack
from symmetry_partners import V
from ensemble import maps
from sectors import point_permutations
from build_data import csvwrite
R=json.loads((ROOT/'responses.json').read_text());data={};rows=[];valid=[]
for name in ['isotropic_Q_pair_average','isotropic_K_pair_average']:
 if name not in R:continue
 r=R[name];Q=np.array(r['Q']);v0=V(Q)
 if name=='isotropic_Q_pair_average':allmaps=json.loads((ROOT/'symmetry_partner_maps.json').read_text())['records'];chosen={tuple(np.round(t['physical_Q'],10)):t for t in allmaps if t['initial_source']=='K'};targets=[(t['physical_Q'],unpack(t['operator_transform']),t['source_after_map']) for t in chosen.values()]
 else:
  p=point_permutations()[3];ix,c,_=maps(p);vp=np.zeros_like(v0);vp[ix]=c[:,None]*v0;qp=np.array([2/3,1/3]);T=V(qp).conj().T@vp;err=la.norm(V(qp)@T-vp);assert err<1e-10;valid.append({'physical_K_to_Kprime_operator_error':float(err)});targets=[(Q,np.eye(6),'K'),(qp,T,'Kprime')]
 for iq,(q,T,mappedK) in enumerate(targets):
  for i,p in enumerate(r['poles']):
   src=mappedK if p.get('source') in ['j1_Q','j1_K'] else ('Kprime' if mappedK=='K' else 'K');gs=np.array([1/3,2/3]) if src=='K' else np.array([2/3,1/3]);target=(gs+q)%1;celln=np.rint(np.mod([target[0]+target[1],2*target[0]-target[1]],1)*3).astype(int)%3
   key=name+f'_star{iq}_pole{i}';M=T@unpack(p['residue'])@T.conj().T;data[key]=M;rows.append({'residue_id':key,'response':name,'Q1':float(q[0]),'Q2':float(q[1]),'omega_over_J':p['omega'],'S_final':.5,'mapped_source':src,'target_physical_q1':float(target[0]),'target_physical_q2':float(target[1]),'target_cell_n1':int(celln[0]),'target_cell_n2':int(celln[1]),'direct_residual_in_original_sector':p['max_direct_H_residual'],'residue_trace':float(np.trace(M).real),'provenance':'Exact unitary symmetry map of independently converged pole; not a separate eigensolve','ground_orbital_weights':'1/2 K + 1/2 Kprime included'})
np.savez_compressed(ROOT/'symmetry_partner_residues.npz',**data);csvwrite('symmetry_partner_poles.csv',rows);dump(ROOT/'symmetry_partner_export_validation.json',{'status':'PASS','rows':len(rows),'checks':valid});print('Exported',len(rows),'symmetry-derived pole residues.')
