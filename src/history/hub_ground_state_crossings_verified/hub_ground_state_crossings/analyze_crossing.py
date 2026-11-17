"""Static correlations for a verified ground state or explicitly labeled competing eigenstate."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,time,resource,csv,argparse,numpy as np,scipy.linalg as la
import cross_scan as cs
from ed import Model
from cyclic_model import CyclicModel
from observables import convert_basis,pair_expectations,correlation_matrices,total_spin_check
from sectors import point_permutations
from refine_crossings import summarize
ROOT=cs.ROOT

def pack(x):a=np.asarray(x);return {'real':a.real.tolist(),'imag':a.imag.tolist()}
def analyze(name,force=False):
 output=ROOT/('obs_'+name+'.json')
 if output.exists() and not force:return json.loads(output.read_text())
 rec=json.loads((ROOT/(name+'.json')).read_text());assert rec['verified'];j=rec['Jp'];proof=summarize(j);assert proof['global_ground_identified']
 z=np.load(ROOT/(name+'.npz'));start=time.perf_counter();m=Model(primitive=False,jp=j);m.sector(z['q']);src=CyclicModel(jp=j,cell_n=tuple(z['cell_n']),angular=int(z['angular']));v=convert_basis(src,m,z['vector']);del src
 assert abs(la.norm(v)-1)<1e-10;v/=la.norm(v);h=m.mv(v);energy=float(np.vdot(v,h).real);res=float(la.norm(h-energy*v));assert res<1e-8
 spin=total_spin_check(m,v);assert spin['S2_eigen_residual']<1e-5 and spin['nearest_half_integer_spin']==rec['spin']
 stat=pair_expectations(m,v);G=stat['Gdot'];Gz=stat['Gzz'];mag=stat['magnetization'];perms=point_permutations(jp=j);Gavg=sum(G[np.ix_(p,p)] for p in perms)/12;C=correlation_matrices(Gavg,3);raw0=C[(0,0)];n=np.array([1,-1,0.]);u=np.array([1,1,0.]);hh=np.array([0,0,1.]);weight=np.exp(-2j*np.pi*np.array([0,1/3,2/3]))/np.sqrt(3)
 selected=correlation_matrices(G,3);conn=selected[(0,0)]-9*np.outer(mag,mag)
 winner=proof['winner'];is_ground=abs(energy-winner['energy'])<1e-8 and (rec['cell_n'],rec['character'],rec['spin'])==(winner['cell_n'],winner['character'],winner['spin'])
 summary={'name':name,'Jp':j,'role':'GROUND' if is_ground else 'COMPETING_VERIFIED_EIGENSTATE','energy':energy,'energy_above_global_ground':energy-winner['energy'],'residual':res,'spin':rec['spin'],'cell_n':rec['cell_n'],'character':rec['character'],'rim_neel_m2':float((n@raw0@n).real/36),'rim_neel_selected_connected_m2':float((n@conn@n).real/36),'hub_S_Gamma':float(raw0[2,2].real),'hub_S_cell01':float(C[(0,1)][2,2].real),'hub_S_cell12':float(C[(1,2)][2,2].real),'rim_hub_uniform_Gamma':float((u@raw0@hh).real),'rim_hub_neel_Gamma_abs':float(abs(n@raw0@hh)),'physical_K_intensity':float(np.vdot(weight,raw0@weight).real),'total_moment_squared':float(spin['S2']/729),'selected_magnetization_per_site':mag.tolist(),'S2_residual':spin['S2_eigen_residual'],'analysis_seconds':time.perf_counter()-start,'pair_kernel_seconds':stat['seconds'],'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'ensemble':'Scalar correlations averaged over the exact point group and the spin multiplet; selected fixed-M connected rim value is separate.'}
 checks={'energy_reconstruction_error':float(abs(sum(b*G[a,d] for a,d,b in m.bonds)-energy)),'S2_pair_vs_action_error':float(abs(G.sum()-spin['S2'])),'Sz_total':float(9*sum(mag)),'Sz_variance':float(Gz.sum()-(9*sum(mag))**2),'Hermiticity_max':float(max(la.norm(x-x.conj().T) for x in C.values())),'minimum_C_eigenvalue':float(min(la.eigvalsh(x)[0] for x in C.values())),'onsite_sum_rule_error':float(max(abs(sum(M[a,a].real for M in C.values())-6.75) for a in range(3))),'inversion_nh_Gamma_abs':summary['rim_hub_neel_Gamma_abs'],'slope_from_correlations_error':float(abs(sum(G[a,b] for a,b,_ in m.bonds if a%3==2 or b%3==2)-rec['slope_dE_dJp']))}
 assert checks['slope_from_correlations_error']<1e-8 and checks['energy_reconstruction_error']<1e-8 and checks['S2_pair_vs_action_error']<1e-7 and abs(checks['Sz_total']-.5)<1e-9 and abs(checks['Sz_variance'])<1e-7 and checks['Hermiticity_max']<1e-8 and checks['minimum_C_eigenvalue']>-1e-8 and checks['onsite_sum_rule_error']<1e-8 and checks['inversion_nh_Gamma_abs']<1e-8
 result={'summary':summary,'validation':checks,'C':{str(k):pack(M) for k,M in C.items()}};output.write_text(json.dumps(result,indent=2));np.savez_compressed(ROOT/('obs_'+name+'_realspace.npz'),Gdot=G,Gzz=Gz,G_pointgroup=Gavg,magnetization=mag)
 with (ROOT/('obs_'+name+'_C.csv')).open('w') as f:
  w=csv.writer(f);w.writerow(['Jp','state','role','n1','n2','alpha','beta','real','imag'])
  for (a,b),M in C.items():
   for i in range(3):
    for k in range(3):w.writerow([j,name,summary['role'],a,b,'ABH'[i],'ABH'[k],M[i,k].real,M[i,k].imag])
 print(json.dumps(summary,indent=2));return result
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('name');a=p.parse_args();analyze(a.name)
