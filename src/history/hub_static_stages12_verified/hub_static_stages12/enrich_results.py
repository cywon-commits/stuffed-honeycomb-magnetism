"""Derived normalized modes, channel matrices and additional positivity checks."""
import pathlib,json,csv,numpy as np,scipy.linalg as la
from observables import coefficients,correlation_matrices
ROOT=pathlib.Path(__file__).resolve().parent
def unpack(z):return np.array(z['real'])+1j*np.array(z['imag'])
def enrich(path):
 d=json.loads(path.read_text());jp=d['summary']['Jp_over_J'];name=f'static_j{jp:g}';mats={k:unpack(v) for k,v in d['C_SU2_pointgroup'].items()};C=mats['(0, 0)'];n=np.array([1,-1,0.])/np.sqrt(2);m=np.array([1,1,0.])/np.sqrt(2);h=np.array([0,0,1.]);B=np.column_stack([m,h]);sub=B.T@C@B;ev,v=la.eigh(sub)
 s=d['summary'];s.update(rim_neel_normalized_mode_eigenvalue=float(np.vdot(n,C@n).real),uniform_hub_mode_lower_eigenvalue=float(ev[0]),uniform_hub_mode_upper_eigenvalue=float(ev[1]),upper_uniform_hub_mode_hub_weight=float(abs(v[1,1])**2),uniform_hub_correlation_coefficient=float((np.sqrt(2)*sub[0,1].real)/np.sqrt((2*sub[0,0].real)*sub[1,1].real)))
 for key,label in [('(0, 1)','cell01'),('(1, 2)','cell12')]:s['hub_S_'+label]=float(mats[key][2,2].real)
 sh={k:float(v[2,2].real) for k,v in mats.items()};maximum=max(sh.values());s['hub_peak_cell_n']=[list(map(int,k.strip('()').split(','))) for k,value in sh.items() if maximum-value<1e-8];s['hub_peak_S']=maximum
 # Compare eigenprojectors rather than arbitrary vectors inside a degenerate eigenspace.
 groups=[]
 for key,Ck in mats.items():
  e,U=la.eigh(Ck);todo=list(range(3))
  while todo:
   first=todo.pop(0);inds=[first]+[j for j in todo if abs(e[j]-e[first])<1e-8];todo=[j for j in todo if j not in inds];P=U[:,inds]@U[:,inds].conj().T
   groups.append({'cell_n':key,'indices':inds,'dimension':len(inds),'eigenvalue_mean':float(np.mean(e[inds])),'hub_projector_trace':float(P[2,2].real),'rim_neel_projector_trace':float(np.vdot(n,P@n).real),'projector_real':P.real.tolist(),'projector_imag':P.imag.tolist()})
 d['degenerate_correlation_eigenspaces']=groups
 G=np.load(ROOT/(name+'_realspace.npz'));zzraw=correlation_matrices(G['Gzz_selected'],3);zz={k:C.copy() for k,C in zzraw.items()};mag=G['magnetization'];zz[(0,0)]-=9*np.outer(mag,mag)
 sets=[mats,{k:unpack(v) for k,v in d['C_selected_raw'].items()},{k:unpack(v) for k,v in d['C_selected_connected'].items()},zzraw,zz]
 d['validation']['all_conventions_Hermiticity_max']=float(max(la.norm(M-M.conj().T) for S in sets for M in S.values()));d['validation']['all_conventions_minimum_eigenvalue']=float(min(la.eigvalsh(M)[0] for S in sets for M in S.values()))
 assert d['validation']['all_conventions_Hermiticity_max']<1e-8 and d['validation']['all_conventions_minimum_eigenvalue']>-1e-8
 with (ROOT/(name+'_channels.csv')).open('w') as f:
  w=csv.writer(f);w.writerow(['Jp_over_J','cell_n1','cell_n2','channel_left','channel_right','real','imag','convention'])
  for a in range(3):
   for b in range(3):
    Ck=mats[str((a,b))];channels=['n','m','H','Phi']+(['physical_K','physical_Kprime'] if (a,b)==(0,0) else [])
    for left in channels:
     for right in channels:
      val=np.vdot(coefficients([a/3,b/3],left),Ck@coefficients([a/3,b/3],right));w.writerow([jp,a,b,left,right,val.real,val.imag,'SU2_pointgroup_average'])
 erim=0j;ehub=0j
 for a in range(3):
  for b in range(3):
   Ck=mats[str((a,b))];k=np.array([a,b])/3;ff=1+np.exp(-2j*np.pi*k[0])+np.exp(-2j*np.pi*k[1]);phi=coefficients(k,'Phi');erim+=ff*Ck[0,1];ehub+=jp*(phi[0]*Ck[2,0]+phi[1]*Ck[2,1])
 d['validation']['Fourier_energy_reconstruction_error']=float(abs(erim+ehub-s['energy']));d['validation']['Fourier_energy_imaginary_abs']=float(abs((erim+ehub).imag));assert abs(erim+ehub-s['energy'])<1e-8
 s['rim_total_spin_squared']=float(9*np.vdot(np.array([1,1,0.]),C@np.array([1,1,0.])).real);s['rim_hub_total_dot']=float(9*np.vdot(np.array([1,1,0.]),C@np.array([0,0,1.])).real);s['rim_energy_from_Fourier']=float(erim.real);s['hub_rim_energy_from_Fourier']=float(ehub.real)
 d['correlation_eigenmode_note']='Individual vectors inside degenerate eigenspaces are basis dependent. Compare saved projectors or inversion-resolved Gamma modes.'
 path.write_text(json.dumps(d,indent=2));return s
if __name__=='__main__':
 for p in ROOT.glob('static_j*.json'):
  if p.name.endswith('_realspace.json'):continue
  print(json.dumps(enrich(p),indent=2))
