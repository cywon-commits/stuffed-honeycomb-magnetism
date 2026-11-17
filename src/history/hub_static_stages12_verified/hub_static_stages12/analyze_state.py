"""Compute all nine full sublattice static correlation matrices from a verified vector."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import argparse,json,pathlib,time,resource,csv
import numpy as np,scipy.linalg as la
from ed import Model
from cyclic_model import CyclicModel
from coverage import verified_ground
from observables import *
from sectors import cell_stars,point_permutations
ROOT=pathlib.Path(__file__).resolve().parent

def pack(a):a=np.asarray(a);return {'real':a.real.tolist(),'imag':a.imag.tolist()}
def analyze(jp,path,primitive=False):
 start=time.perf_counter()
 if jp!=1:
  proof=verified_ground(jp)
  if proof is None:raise ValueError('Complete converged momentum/character coverage is required before labeling this state a ground state.')
  if pathlib.Path(path).stem!=proof['ground_name']:raise ValueError('Checkpoint does not match the independently determined ground-state block.')
 z=np.load(path);
 if 'jp' in z:assert abs(float(z['jp'])-jp)<1e-12
 m=Model(primitive=False,jp=jp);q=np.array(z['q']);m.sector(q)
 if primitive or ('cyclic' in z and bool(z['cyclic'])):
  if primitive:src=Model();src.sector(q)
  else:src=CyclicModel(jp=jp,cell_n=tuple(z['cell_n']),angular=int(z['angular']))
  v=convert_basis(src,m,z['vector']);del src
 else:v=np.ascontiguousarray(z['vector'])
 norm=la.norm(v);assert abs(norm-1)<1e-10;v/=norm
 hv=m.mv(v);e=float(np.vdot(v,hv).real);res=float(la.norm(hv-e*v));assert res<1e-8,(jp,e,res)
 if jp==1:
  reference=json.loads((ROOT/'isotropic_ground_verification.json').read_text())
  if abs(e-reference['energy'])>1e-7:raise ValueError('Isotropic vector energy does not match the previously verified global ground state.')
 spin=total_spin_check(m,v);assert spin['S2_eigen_residual']<1e-5,spin
 stats=pair_expectations(m,v);G=stats['Gdot'];Gz=stats['Gzz'];mag=stats['magnetization'];Nc=m.L**2
 selected=correlation_matrices(G,m.L);selectedz=correlation_matrices(Gz,m.L);connected={k:C.copy() for k,C in selected.items()};connected[(0,0)]-=Nc*np.outer(mag,mag)
 selectedzconn={k:C.copy() for k,C in selectedz.items()};selectedzconn[(0,0)]-=Nc*np.outer(mag,mag)
 perms=point_permutations(jp=jp);Gavg=sum(G[np.ix_(p,p)] for p in perms)/len(perms);Cavg=correlation_matrices(Gavg,m.L)
 _,inv=cell_stars(jp=jp);cellk=np.array([q[0]+q[1],2*q[0]-q[1]])%1
 invinfo={'averaging':'SU(2) spin-multiplet scalar average and 12-element exact point-group average','single_state_cell_k':cellk.tolist()}
 if np.max(abs(np.exp(4j*np.pi*(m.gs@m.q))-1))<1e-8:
  iv=inversion_vector(m,v,inv);ov=np.vdot(v,iv);invinfo.update(overlap=pack(ov),parity_eigen_residual=float(min(la.norm(iv-v),la.norm(iv+v))))
 else:invinfo['single_state_inversion_partner']='Different cell momentum; orthogonal by translation symmetry.'
 n=np.array([1,-1,0.]);h=np.array([0,0,1.]);unif=np.array([1,1,0.]);ones=np.ones(3)
 raw0=Cavg[(0,0)];S2=float(G.sum());M=float(Nc*mag.sum());varz=float(Gz.sum()-M*M)
 energy_reconstructed=sum(j*G[a,b] for a,b,j in m.bonds)
 checks={'source_norm_error':float(abs(norm-1)),'energy_reconstruction_error':float(abs(energy_reconstructed-e)),'S2_pair_vs_action_error':float(abs(S2-spin['S2'])),'Sz_total':M,'Sz_variance':varz,'connected_Czz_total_Gamma':float(abs(ones@selectedzconn[(0,0)]@ones)),'dot_total_Gamma_vs_S2_error':float(abs(ones@raw0@ones-S2/Nc)),'sublattice_onsite_sum_rule_error':float(max(abs(sum(Cavg[k][a,a].real for k in Cavg)-.75*Nc) for a in range(3))),'Hermiticity_max':float(max(la.norm(C-C.conj().T) for C in Cavg.values())),'minimum_C_eigenvalue':float(min(la.eigvalsh(C)[0] for C in Cavg.values())),'inversion_nh_Gamma_abs':float(abs(n@raw0@h)),'inversion_nm_Gamma_abs':float(abs(n@raw0@unif)),'selected_inversion_nh_Gamma_abs':float(abs(n@selected[(0,0)]@h))}
 for key in ['energy_reconstruction_error','S2_pair_vs_action_error','connected_Czz_total_Gamma','dot_total_Gamma_vs_S2_error','sublattice_onsite_sum_rule_error','Hermiticity_max','inversion_nh_Gamma_abs','inversion_nm_Gamma_abs']:
  assert checks[key]<2e-7,(key,checks[key])
 assert checks['minimum_C_eigenvalue']>-1e-8 and abs(varz)<1e-7 and abs(M-.5)<1e-10
 wK=coefficients([0,0],'physical_K');wKp=coefficients([0,0],'physical_Kprime')
 rim=[G[a,b] for a,b,j in m.bonds if a%3!=2 and b%3!=2];hubrim=[G[a,b] for a,b,j in m.bonds if a%3==2 or b%3==2]
 hubs=list(range(2,m.n,3));hubpairs=[G[a,b] for ii,a in enumerate(hubs) for b in hubs[ii+1:]]
 summary={'Jp_over_J':jp,'energy':e,'energy_residual':res,'S_total':spin['nearest_half_integer_spin'],'S2':S2,'S2_eigen_residual':spin['S2_eigen_residual'],'cell_ground_momentum':cellk.tolist(),'cell_basis_dimension':m.nb,'rim_neel_m2':float((n@raw0@n).real*Nc/(2*Nc)**2),'rim_neel_connected_selected_m2':float((n@connected[(0,0)]@n).real*Nc/(2*Nc)**2),'hub_S_Gamma':float(raw0[2,2].real),'hub_total_spin_squared':float(Nc*raw0[2,2].real),'rim_hub_uniform_Gamma':float((unif@raw0@h).real),'C_AH_Gamma':float(raw0[0,2].real),'C_BH_Gamma':float(raw0[1,2].real),'rim_hub_neel_Gamma_abs':float(abs(n@raw0@h)),'physical_K_120_intensity':float(np.vdot(wK,raw0@wK).real),'physical_Kprime_120_intensity':float(np.vdot(wKp,raw0@wKp).real),'ferromagnetic_m2':S2/m.n**2,'rim_NN_mean':float(np.mean(rim)),'hub_rim_bond_mean':float(np.mean(hubrim)),'hub_hub_all_distinct_mean':float(np.mean(hubpairs)),'single_state_sublattice_Sz_per_site':mag.tolist(),'pair_kernel_seconds':stats['seconds'],'analysis_seconds':time.perf_counter()-start,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'state_convention':'Primary dot matrices: point-group and spin-multiplet averaged. Also saved: selected Sz=+1/2 raw, connected and zz matrices.'}
 rows=[];modes=[]
 for k,C in Cavg.items():
  vals,vecs=la.eigh(C)
  for j,ev in enumerate(vals):
   mode=vecs[:,j];pivot=np.argmax(abs(mode));mode*=np.exp(-1j*np.angle(mode[pivot]));modes.append({'n1':k[0],'n2':k[1],'eigenmode':j,'eigenvalue':float(ev),'hub_weight':float(abs(mode[2])**2),'rim_neel_weight':float(abs(mode[0]-mode[1])**2/2),'rim_uniform_weight':float(abs(mode[0]+mode[1])**2/2),'vector':pack(mode)})
  for label,mats in [('SU2_pointgroup_average',Cavg),('selected_raw_dot',selected),('selected_connected_dot',connected),('selected_raw_zz',selectedz),('selected_connected_zz',selectedzconn)]:
   for a in range(3):
    for b in range(3):rows.append([jp,k[0],k[1],k[0]/3,k[1]/3,label,'ABH'[a],'ABH'[b],mats[k][a,b].real,mats[k][a,b].imag])
 output={'summary':summary,'validation':checks,'spin_check':spin,'inversion':invinfo,'C_SU2_pointgroup':{str(k):pack(C) for k,C in Cavg.items()},'C_selected_raw':{str(k):pack(C) for k,C in selected.items()},'C_selected_connected':{str(k):pack(C) for k,C in connected.items()},'correlation_eigenmodes':modes}
 name=f'static_j{jp:g}';(ROOT/(name+'.json')).write_text(json.dumps(output,indent=2))
 np.savez_compressed(ROOT/(name+'_realspace.npz'),G_selected=G,Gzz_selected=Gz,G_pointgroup=Gavg,magnetization=mag)
 with (ROOT/(name+'_C.csv')).open('w') as f:
  writer=csv.writer(f);writer.writerow(['Jp_over_J','n1','n2','cell_k1','cell_k2','convention','alpha','beta','real','imag']);writer.writerows(rows)
 print(json.dumps(output['summary'],indent=2));print(json.dumps(checks,indent=2))
 return output
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--jp',type=float,required=True);p.add_argument('--file',required=True);p.add_argument('--primitive',action='store_true');a=p.parse_args();analyze(a.jp,a.file,a.primitive)
