"""Cached reproducible scientific tables; reads measured moments and verified subspaces."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,csv,pathlib,numpy as np,scipy.linalg as la
from experiment import ROOT,pack,dump,geometry_weights
from group_poles import unpack
from ensemble import unfold,fourier

def csvwrite(path,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with (ROOT/path).open('w') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def canonical(groups,name):
 gs=[g for g in groups if g['name']==name and g['verified_eigenpair'] and g['scalar_accessible'] and g['residue_trace']>1e-14]
 if name=='j1_G':gs=[g for g in gs if ('_a0_' in g['tag']) or (g['tag']==name and np.trace(unpack(g['residue'])[:4,:4]).real>1e-12)]
 return gs

def merged(groups,name):
 out=[]
 for g in sorted(canonical(groups,name),key=lambda x:x['omega']):
  if out and abs(out[-1]['omega']-g['omega'])<1e-8:
   p=out[-1];p['R']+=unpack(g['residue']);p['multiplicity']+=g['multiplicity'];p['max_direct_H_residual']=max(p['max_direct_H_residual'],g['max_direct_H_residual']);p['max_S2_residual']=max(p['max_S2_residual'],g['max_S2_residual']);p['tags'].append(g['tag'])
  else:out.append({'omega':g['omega'],'energy':g['energy'],'R':unpack(g['residue']),'multiplicity':g['multiplicity'],'max_direct_H_residual':g['max_direct_H_residual'],'max_S2_residual':g['max_S2_residual'],'tags':[g['tag']]})
 return out

def main():
 groups=json.loads((ROOT/'degenerate_pole_groups.json').read_text())['groups'];mom={}
 for f in sorted(ROOT.glob('*_moments.json')):
  d=json.loads(f.read_text())
  if 'name' in d:mom[d['name']]=d
 static=[];proj=[];cost=[];rawpol=[];residues={};responses={};convergence=[]
 for name,d in mom.items():
  channels=d['channel_names'];n=len(channels)
  for a in range(n):
   for b in range(n):
    for key in ['Mraw','M0_elastic','M0','M1','M2']:
     z=unpack(d[key])[a,b];static.append({'experiment':name,'Jprime_over_J':d['j'],'Q_or_cell_k1':d['Q'][0],'Q_or_cell_k2':d['Q'][1],'channel_a':channels[a],'channel_b':channels[b],'moment':key,'real':float(z.real),'imag':float(z.imag)})
    x=unpack(d['Mraw'])[a,b];e=unpack(d['M0_elastic'])[a,b];pos=unpack(d['M0'])[a,b];proj.append({'experiment':name,'channel_a':channels[a],'channel_b':channels[b],'M0_raw_real':x.real,'M0_raw_imag':x.imag,'zero_energy_real':e.real,'zero_energy_imag':e.imag,'positive_energy_real':pos.real,'positive_energy_imag':pos.imag})
  rows=merged(groups,name);responses[name]={'moments':{k:pack(unpack(d[k])) for k in ['Mraw','M0_elastic','M0','M1','M2']},'poles':[],'channel_names':channels,'j':d['j'],'Q':d['Q'],'target_q':d['target_q']}
  for i,p in enumerate(rows):
   key=name+f'_pole_{i}';residues[key]=p['R'];r={'pole_id':key,'experiment':name,'Jprime_over_J':d['j'],'S_ground':.5,'source_label':('Kprime' if d['source_partner'] else 'K') if d['j']==1 else 'cell_Gamma_C6_character3','transfer_1':d['Q'][0],'transfer_2':d['Q'][1],'transfer_convention':'physical primitive' if d['j']==1 else 'three-site cell','target_cell_n1':int(round(3*(d['target_q'][0]+d['target_q'][1])))%3,'target_cell_n2':int(round(3*(2*d['target_q'][0]-d['target_q'][1])))%3,'target_q_meaning':'primitive momentum' if d['j']==1 else 'cell-character reciprocal-coordinate root','target_q1':d['target_q'][0],'target_q2':d['target_q'][1],'target_cell_k1':(d['target_q'][0]+d['target_q'][1])%1,'target_cell_k2':(2*d['target_q'][0]-d['target_q'][1])%1,'omega_over_J':p['omega'],'E_final_over_J':p['energy'],'S_final':.5,'multiplicity':p['multiplicity'],'direct_H_residual':p['max_direct_H_residual'],'S2_residual':p['max_S2_residual'],'total_channel_trace':float(np.trace(p['R']).real),'symmetry_blocks':';'.join(p['tags']),'status':'VERIFIED_SCALAR_EIGENSPACE_RESIDUE'};rawpol.append(r);responses[name]['poles'].append({**{k:v for k,v in p.items() if k!='R'},'residue':pack(p['R'])})
  for p in range(3):
   exact=unpack(d[f'M{p}']);partial=sum((r['omega']**p*r['R'] for r in rows),np.zeros_like(exact));rem=(exact-partial);responses[name][f'M{p}_verified_partial']=pack(partial);responses[name][f'M{p}_unresolved']=pack(rem)
  cost.append({'experiment':name,'phase':'seed_and_exact_moments','dimension':d['dimension'],'H_calls':d['H_calls'],'wall_seconds':d['wall_seconds'],'peak_process_RSS_MiB':d['peak_RSS_MiB']})
 for f in sorted(ROOT.glob('*_spectrum.json')):
  d=json.loads(f.read_text());tag=d.get('tag',d['name']);npoles=len(d['poles']);its=d.get('iterative_H_calls',d['H_calls_including_spin_checks']-2*npoles);cost.append({'experiment':tag,'phase':d.get('method','iterative_eigensolve'),'dimension':d['dimension'],'H_calls':d['H_calls_including_spin_checks'],'Arnoldi_H_iterations':its,'iteration_convention':'one H action per Arnoldi step; restart cycles not exposed by scipy','CPU_seconds':d.get('CPU_seconds'),'wall_seconds':d['wall_seconds'],'peak_process_RSS_MiB':d['peak_RSS_MiB']})
  for p in d['poles']:convergence.append({'tag':tag,'eigenpair_index':p['index'],'omega':p['omega'],'spin':p['spin']['nearest_half_integer_spin'],'direct_H_residual':p['H_residual'],'S2_residual':p['spin']['S2_eigen_residual'],'verified':p['verified_eigenpair'],'scalar_accessible':p['scalar_accessible'],'weight':p['total_seed_weight']})
 for f in sorted(ROOT.glob('*_quadrature.json')):
  d=json.loads(f.read_text());cost.append({'experiment':d['name'],'phase':'BLOCK_KRYLOV_QUADRATURE_NOT_POLES','H_calls':d['H_calls'],'Arnoldi_H_iterations':d['H_calls'],'CPU_seconds':d.get('CPU_seconds'),'wall_seconds':d['wall_seconds'],'peak_process_RSS_MiB':d['peak_RSS_MiB']})
 # Symmetry-restored isotropic pair average. Scalars are M-independent; average only orbital sources.
 T=np.eye(6,dtype=complex)[[2,3,0,1,4,5]];T[4:]*=-1
 if 'j1_G' in responses:
  d=responses['j1_G'];r={'j':1,'Q':[0,0],'channel_names':d['channel_names'],'poles':[],'ensemble':'equal verified K/Kprime orbital weights; fixed M=1/2'}
  for key in ['Mraw','M0_elastic','M0','M1','M2']:
   X=unpack(d['moments'][key]);r[key]=pack((X+T@X.conj()@T.conj().T)/2)
  for p in d['poles']:
   X=unpack(p['residue']);r['poles'].append({**{k:v for k,v in p.items() if k!='residue'},'residue':pack((X+T@X.conj()@T.conj().T)/2)})
  responses['isotropic_G_pair_average']=r
 if 'j1_Q' in responses and 'j1_Qprime' in responses:
  r={'j':1,'Q':[8/9,7/9],'channel_names':responses['j1_Q']['channel_names'],'poles':[],'ensemble':'equal verified K/Kprime orbital weights; sources deflated separately'}
  for key in ['Mraw','M0_elastic','M0','M1','M2']:r[key]=pack((unpack(responses['j1_Q']['moments'][key])+unpack(responses['j1_Qprime']['moments'][key]))/2)
  entries=[]
  for name in ['j1_Q','j1_Qprime']:
   for p in responses[name]['poles']:entries.append({**{k:v for k,v in p.items() if k!='residue'},'residue':pack(unpack(p['residue'])/2),'source':name})
  r['poles']=sorted(entries,key=lambda x:x['omega']);responses['isotropic_Q_pair_average']=r
 if 'j1_K' in responses and 'j1_KprimeSource' in responses:
  r={'j':1,'Q':[1/3,2/3],'channel_names':responses['j1_K']['channel_names'],'poles':[],'ensemble':'equal verified K/Kprime orbital weights; sources deflated separately; partial verified poles'}
  for key in ['Mraw','M0_elastic','M0','M1','M2']:r[key]=pack((unpack(responses['j1_K']['moments'][key])+unpack(responses['j1_KprimeSource']['moments'][key]))/2)
  for name in ['j1_K','j1_KprimeSource']:
   for pole in responses[name]['poles']:r['poles'].append({**{k:v for k,v in pole.items() if k!='residue'},'residue':pack(unpack(pole['residue'])/2),'source':name})
  r['poles']=sorted(r['poles'],key=lambda x:x['omega']);responses['isotropic_K_pair_average']=r
 ensemble_rows=[];crossrows=[]
 for name,r in responses.items():
  if not name.startswith('isotropic_'):continue
  for i,pole in enumerate(r['poles']):
   key=name+f'_pole_{i}';residues[key]=unpack(pole['residue']);ensemble_rows.append({'pole_id':key,'response':name,'omega':pole['omega'],'S_final':.5,'H_residual':pole['max_direct_H_residual'],'S2_residual':pole['max_S2_residual'],'weight_factor':'orbital ground-pair average already included','source':pole.get('source','K plus exact conjugate partner')})
   for a in range(6):
    for b in range(6):crossrows.append({'pole_id':key,'omega':pole['omega'],'channel_a':r['channel_names'][a],'channel_b':r['channel_names'][b],'residue_real':residues[key][a,b].real,'residue_imag':residues[key][a,b].imag})
 csvwrite('ensemble_poles.csv',ensemble_rows);csvwrite('complex_residues.csv',crossrows)
 # Fixed physical face form factors; no primitive block diagonalization away from isotropy.
 modes=[];bounds=[];mode_poles=[];ref=np.load(ROOT/'inputs/raw_observables_j1.npz')
 def reference_u(Q,kind):
  off=0 if kind=='D' else 108;data=unfold(fourier(ref['ensemble_full_connected'][off:off+54,off:off+54]),kind)[0];d=min(data,key=lambda r:la.norm(r['Q']-Q));return d['eigenvectors'][:,-1]
 configurations=[('isotropic_G_pair_average',np.array([0,0])),('isotropic_Q_pair_average',np.array([8/9,7/9])),('isotropic_K_pair_average',np.array([1/3,2/3])),('j09_G',np.array([0,0])),('j09_q',np.array([8/9,7/9])),('j07_G',np.array([0,0]))]
 for name,Q in configurations:
  if name not in responses:continue
  r=responses[name];mm=r.get('moments',r);size=len(r['channel_names']);P=geometry_weights(Q) if size==18 else np.eye(6)
  for kind,off in [('D',0),('chi',4)]:
   u=np.zeros(6,complex);u[off:off+2]=reference_u(Q,kind);u=P@u;residues[name+'_'+kind+'_reference_mode']=u;mvals=[float(np.vdot(u,unpack(mm[f'M{p}'])@u).real) for p in range(3)];acc=[0.,0.,0.];ps=[]
   for p in r['poles']:
    wgt_raw=float(np.vdot(u,unpack(p['residue'])@u).real);assert wgt_raw>=-1e-12;wgt=max(0.,wgt_raw);acc=[v+p['omega']**k*wgt for k,v in enumerate(acc)];ps.append((p['omega'],wgt));mode_poles.append({'response':name,'mode':kind,'omega':p['omega'],'weight':wgt,'H_residual':p['max_direct_H_residual'],'S2_residual':p['max_S2_residual'],'status':'VERIFIED_EIGENSPACE_OR_EXACT_SYMMETRY_PARTNER'})
   remaining=[a-b for a,b in zip(mvals,acc)];modes.append({'response':name,'Jprime_over_J':r['j'],'probe_Q1':Q[0],'probe_Q2':Q[1],'mode':kind,'M0_positive_energy':mvals[0],'M1':mvals[1],'M2':mvals[2],'centroid':mvals[1]/mvals[0],'resolved_M0':acc[0],'resolved_fraction':acc[0]/mvals[0],'unresolved_M0':remaining[0],'unresolved_M1':remaining[1],'unresolved_M2':remaining[2],'elastic_M0':max(0.,float(np.vdot(u,unpack(mm['M0_elastic'])@u).real)),'mode_definition':'fixed normalized isotropic reference leading face eigenvector, embedded in full cell fields when j!=1'})
   for cutoff in [.1,.5,1,2,3,4,6]:
    lower=max(0.,remaining[0]-remaining[1]/cutoff,remaining[0]-remaining[2]/cutoff**2);above=min(remaining[0],remaining[1]/cutoff,remaining[2]/cutoff**2);bounds.append({'response':name,'mode':kind,'cutoff_over_J':cutoff,'unidentified_total_weight':remaining[0],'unidentified_below_cutoff_lower_bound':lower,'unidentified_above_cutoff_upper_bound':above,'assumption':'positive deflated spectral measure; floating point verified moments, not interval arithmetic'})
 all_remainders={}
 for name,r in responses.items():
  M=r.get('moments',r);nc=len(r['channel_names']);all_remainders[name]={}
  for power in range(3):
   exact=unpack(M[f'M{power}']);partial=sum((p['omega']**power*unpack(p['residue']) for p in r['poles']),np.zeros((nc,nc),complex));all_remainders[name][f'M{power}']=pack(exact-partial)
  r['unresolved_moment_matrices']=all_remainders[name]
 dump(ROOT/'spectral_moment_remainders.json',all_remainders)
 csvwrite('static_moments.csv',static);csvwrite('ground_projection.csv',proj);csvwrite('poles.csv',rawpol);csvwrite('mode_convergence.csv',convergence);csvwrite('computational_costs.csv',cost);csvwrite('mode_summary.csv',modes);csvwrite('mode_poles.csv',mode_poles);csvwrite('remaining_weight_bounds.csv',bounds);dump(ROOT/'responses.json',responses);dump(ROOT/'static_moments.json',mom);np.savez_compressed(ROOT/'residues.npz',**residues);print('Exported',len(rawpol),'verified scalar eigenspaces,',len(static),'complex moment entries')
if __name__=='__main__':main()
