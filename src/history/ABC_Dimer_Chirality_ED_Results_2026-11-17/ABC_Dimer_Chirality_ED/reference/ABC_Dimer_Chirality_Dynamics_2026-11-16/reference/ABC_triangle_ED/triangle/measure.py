import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,time,resource,argparse,csv,itertools,numpy as np,scipy.linalg as la
from scipy.linalg.blas import zgemm
from ed import Model
from cyclic_model import CyclicModel
from observables import convert_basis,cell_q,total_spin_check,pair_expectations
from local_geometry import faces,matrices,D,chi,P,I,B
from ensemble import unpack_rho,realspace,ensemble,fourier,unfold
from native import seeds
ROOT=pathlib.Path(__file__).resolve().parent

def pack(x):return {'real':np.asarray(x).real.tolist(),'imag':np.asarray(x).imag.tolist()}
def source(j):
 if j==1:
  path=ROOT.parent/'K_primitive.npz';m=Model(primitive=True,jp=1);z=np.load(path);m.sector(z['q']);label='verified_primitive_K'
 elif j in [.2,.5]:
  r=0 if j==.2 else 3;n=(1,2) if j==.2 else (0,0);path=ROOT.parent/f'static/j{j:g}_n{n[0]}{n[1]}_r{r}.npz';z=np.load(path);m=CyclicModel(jp=j,cell_n=n,angular=r);label=path.stem
 else:
  g=json.loads((ROOT/f'ground_j{j:g}.json').read_text());assert g['global_ground_identified'];rec=g['winner'];path=ROOT/(rec['name']+'.npz');z=np.load(path);m=CyclicModel(jp=j,cell_n=tuple(rec['cell_n']),angular=rec['character']);label=path.stem
 q=np.asarray(z['q']);k=np.mod([q[0]+q[1],2*q[0]-q[1]],1);k[np.isclose(k,1)]=0;target=Model(primitive=False,jp=j);target.sector(cell_q(k));v=convert_basis(m,target,z['vector']);del m;v/=la.norm(v);return target,v,float(z['energy']),k,label

def write_csv(name,rows):
 with (ROOT/name).open('w') as f:
  keys=list(dict.fromkeys(k for r in rows for k in r));wr=csv.DictWriter(f,fieldnames=keys);wr.writeheader();wr.writerows(rows)
def run(j,force=False):
 output=ROOT/f'observables_j{j:g}.json'
 if output.exists() and not force:return
 start=time.perf_counter();s,v,E,k0,label=source(j);hv=s.mv(v);res=float(la.norm(hv-E*v));spin=total_spin_check(s,v);assert res<1e-7 and spin['S2_eigen_residual']<1e-5;target=Model(primitive=False,jp=j);G=[];rho=None;mean=None;kernel_time=0
 for n1,n2 in itertools.product(range(3),repeat=2):
  q=np.array([n1,n2])/3;target.sector(cell_q((k0+q)%1));t=time.perf_counter();out,r=seeds(s,target,v,q,rdm=n1==n2==0);kernel_time+=time.perf_counter()-t;gram=zgemm(1.,out.T,out.T,trans_a=2);assert la.norm(gram-gram.conj().T)<1e-8 and la.eigvalsh(gram)[0]>-1e-8;G.append(gram)
  if n1==n2==0:rho=r;mean=out@v.conj()/3
  del out
  print(json.dumps({'Jp':j,'q':[n1,n2],'elapsed_seconds':time.perf_counter()-start}),flush=True)
 G=np.array(G);raw=realspace(G);means=np.concatenate([np.tile(mean[a*6:a*6+6],9) for a in range(3)]);selected=raw-means.conj()[:,None]*means[None,:];avg,avgmeans,avgrho=ensemble(raw,means,rho,j);cov=avg-avgmeans.conj()[:,None]*avgmeans[None,:];SD=fourier(cov[:54,:54]);SC=fourier(cov[108:,108:]);checks={'H_residual':res,'spin':spin,'selected_covariance_Hermiticity':float(la.norm(selected-selected.conj().T)),'ensemble_covariance_Hermiticity':float(la.norm(cov-cov.conj().T)),'selected_covariance_minimum_eigenvalue':float(la.eigvalsh((selected+selected.conj().T)/2)[0]),'ensemble_covariance_minimum_eigenvalue':float(la.eigvalsh((cov+cov.conj().T)/2)[0])};assert min(checks['selected_covariance_minimum_eigenvalue'],checks['ensemble_covariance_minimum_eigenvalue'])>-1e-8
 local=[];matrows=[];identity_errors=[];rhos_twirl=[]
 for which,rs in [('selected',rho),('ensemble',avgrho)]:
  for mu,r in enumerate(rs):
   bonds,c,tw,sig,summary=unpack_rho(r);rhos_twirl.append(tw);metric={'trace_error':float(abs(np.trace(r)-1)),'Hermiticity_error':float(la.norm(r-r.conj().T)),'minimum_rho_eigenvalue':float(la.eigvalsh((r+r.conj().T)/2)[0])};assert metric['trace_error']<1e-8 and metric['Hermiticity_error']<1e-8 and metric['minimum_rho_eigenvalue']>-1e-8
   d=matrices()[mu,0];cmean=(-1)**mu*c;local.append({'Jp':j,'state':label,'ensemble':which,'orientation_mu':mu,'AB':bonds[0],'BC':bonds[1],'CA':bonds[2],'D_mean_real':float(np.trace(r@d).real),'D_mean_imag':float(np.trace(r@d).imag),'chi_CCW_mean':cmean,**summary,**metric});matrows.append({'ensemble':which,'mu':mu,'rho':pack(r),'SU2_twirl':pack(tw),'unnormalized_doublet_flavor':pack(sig),'conditional_doublet_flavor':pack(sig/summary['P_doublet'])})
   if which=='selected':
    diagD=raw[mu,mu];diagDD=raw[54+mu,54+mu];diagC=raw[108+mu,108+mu];p=summary['P_doublet'];identity_errors.append(max(abs(diagC-3*p/16),abs((diagD+diagDD)/2-9*p/8),abs(diagD-diagDD-3*np.sqrt(3)*cmean),abs(mean[mu]-np.trace(r@d)),abs(mean[12+mu]-cmean)))
   else:
    assert abs(avgmeans[mu]-np.trace(r@d))<1e-8 and abs(avgmeans[108+mu]-cmean)<1e-8
 checks['maximum_local_27spin_identity_error']=float(max(identity_errors));assert checks['maximum_local_27spin_identity_error']<1e-8
 stat=pair_expectations(s,v);Gspin=stat['Gdot'];coords,bonds_list=s.coords,s.bonds;checks['energy_from_two_point_error']=float(abs(sum(J*Gspin[a,b] for a,b,J in bonds_list)-E));assert checks['energy_from_two_point_error']<1e-7
 f,o,centers,verts=faces();pairrows=[]
 for t in range(54):
  for u in range(54):
   shared=len(set(f[t])&set(f[u]));category='identical' if t==u else ('overlap_'+str(shared) if shared else 'disjoint');delta=centers[u]-centers[t];dist=min(np.linalg.norm([d[0]+d[1]/2,np.sqrt(3)*d[1]/2]) for a,b in itertools.product([-1,0,1],repeat=2) for d in [delta+a*np.array([3,3])+b*np.array([6,-3])]);pairrows.append({'Jp':j,'t':t,'tprime':u,'category':category,'shared_sites':shared,'centroid_distance':dist,'D_real':cov[t,u].real,'D_imag':cov[t,u].imag,'chi_real':cov[108+t,108+u].real,'chi_imag':cov[108+t,108+u].imag,'selected_D_real':selected[t,u].real,'selected_D_imag':selected[t,u].imag,'selected_chi_real':selected[108+t,108+u].real,'selected_chi_imag':selected[108+t,108+u].imag})
 structure=[]
 for kind,S in [('D',SD),('chi',SC)]:
  for iq,n in enumerate(itertools.product(range(3),repeat=2)):
   e,U=la.eigh(S[iq]);structure.append({'kind':kind,'cell_n':n,'matrix':pack(S[iq]),'eigenvalues':e.tolist(),'eigenvectors':pack(U)})
 physical=[]
 if j==1:
  for kind,S in [('D',SD),('chi',SC)]:
   data,valid=unfold(S,kind);valid['geometric_dimer_gauge']=pack(valid['geometric_dimer_gauge']);checks['primitive_'+kind]=valid;assert valid['unfolding_unitarity_error']<1e-8 and valid['different_primitive_momentum_block_error']<1e-7
   for r in data:
    Q=r['Q'];phase=np.exp(-2j*np.pi*(coords@Q))/np.sqrt(27);mag=float(np.vdot(phase,Gspin@phase).real);physical.append({'kind':kind,'physical_Q':Q.tolist(),'cell_n':r['cell_n'],'matrix':pack(r['matrix']),'eigenvalues':r['eigenvalues'].tolist(),'eigenvectors':pack(r['eigenvectors']),'magnetic_Sdot_Q':mag})
 write_csv(f'local_rdm_j{j:g}.csv',local);write_csv(f'intertriangle_j{j:g}.csv',pairrows);np.savez_compressed(ROOT/f'raw_observables_j{j:g}.npz',selected_rho_ABC=rho,ensemble_rho_ABC=avgrho,selected_operator_Gram_q=G,selected_full_raw_Gram=raw,selected_full_connected=selected,ensemble_full_raw_Gram=avg,ensemble_full_connected=cov,selected_means=means,ensemble_means=avgmeans,D_structure=SD,chi_structure=SC,Gdot=Gspin,faces_ABC=f,faces_CCW=o,face_centroids=centers)
 result={'Jp':j,'source':label,'energy':E,'cell_momentum':k0.tolist(),'local':local,'rdm_matrices':matrows,'cell_structure_factors':structure,'isotropic_primitive_structure_factors':physical,'validation':checks,'wall_seconds':time.perf_counter()-start,'kernel_seconds':kernel_time,'peak_rss_MiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,'ensemble':'Selected fixed-M rho and scalar correlations retained; ensemble uses exact spin twirl, point-group and conjugate-source average; isotropy also primitive-translation cosets.'};output.write_text(json.dumps(result,indent=2));print(json.dumps({'completed':j,'validation':checks,'wall_seconds':result['wall_seconds']}),flush=True)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('j',type=float);p.add_argument('--force',action='store_true');a=p.parse_args();run(a.j,a.force)
