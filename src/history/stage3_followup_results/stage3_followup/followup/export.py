import pathlib,json,csv,numpy as np,scipy.linalg as la
ROOT=pathlib.Path(__file__).resolve().parent
def mat(d):return np.array(d['real'])+1j*np.array(d['imag'])
def pack(x):return {'real':np.asarray(x).real.tolist(),'imag':np.asarray(x).imag.tolist()}
def write(n,rows):
 if rows:
  keys=list(dict.fromkeys(k for r in rows for k in r))
  with (ROOT/n).open('w') as f:w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows({k:json.dumps(v) if isinstance(v,(dict,list,tuple)) else v for k,v in r.items()} for r in rows)
def collect():
 base=json.loads((ROOT.parent/'dynamics/weak_spectrum.json').read_text());records=[]
 for i,d in enumerate(base['poles']):records.append({**d,'name':f'known_2375_r{i}','Jp':.2375,'reflection_character':i,'reflection_parity':(-1)**i,'origin':'PREVIOUS_VERIFIED_STATE'})
 for p in ROOT.glob('j0.*_r*_root*.json'):
  d=json.loads(p.read_text())
  if d.get('verified'):records.append({**d,'origin':'NEW_DEFLATED_SOLVE'})
 p=ROOT/'matched_j0.2.json'
 if p.exists():
  ref=json.loads(p.read_text());records.extend({**d,'origin':'PREVIOUS_VERIFIED_EIGENSTATE_FRESH_TARGET_RESIDUE'} for d in ref['poles'])
 else:ref=None
 records.sort(key=lambda d:(d['Jp'],d['energy'],d['name']))
 rows=[];res=[]
 for d in records:
  R=mat(d['residue_matrix']);assert np.linalg.eigvalsh((R+R.conj().T)/2)[0]>-1e-9
  rows.append({'name':d['name'],'Jp':d['Jp'],'omega_J':d['omega'],'omega_over_Jp_squared':d['omega']/d['Jp']**2,'energy_J':d['energy'],'S':d['spin']['nearest_half_integer_spin'],'reflection_parity':d['reflection_parity'],'H_residual':d['H_residual'],'S2_residual':d['spin']['S2_eigen_residual'],'R_nn':R[0,0].real,'R_mm':R[1,1].real,'R_hh':R[2,2].real,'R_nh_real':R[0,2].real,'R_nh_imag':R[0,2].imag,'R_nn_over_R_hh':R[0,0].real/R[2,2].real if R[2,2].real>1e-25 else None,'origin':d['origin']})
  for a,an in enumerate(['n','m','h']):
   for b,bn in enumerate(['n','m','h']):res.append({'name':d['name'],'Jp':d['Jp'],'omega_J':d['omega'],'channel':an+bn,'real':R[a,b].real,'imag':R[a,b].imag})
 rows.sort(key=lambda x:(x['Jp'],x['omega_J']));write('verified_poles.csv',rows);write('complex_residues.csv',res);(ROOT/'verified_poles.json').write_text(json.dumps(records,indent=2))
 bounds=[];remainder={}
 for j,d in [(.2375,base),(.2,ref)]:
  if d is None:continue
  known=[r for r in records if r['Jp']==j];M=[mat(d['moments']['1'][f'M{p}']) for p in range(3)];sums=[sum((r['omega']**p*mat(r['residue_matrix']) for r in known),np.zeros((3,3),complex)) for p in range(3)];R=[M[p]-sums[p] for p in range(3)];mins=[float(np.linalg.eigvalsh((x+x.conj().T)/2)[0]) for x in R];assert min(mins)>-1e-8
  remainder[str(j)]={'exact_moments':{f'M{p}':pack(M[p]) for p in range(3)},'verified_pole_moment_sums':{f'M{p}':pack(sums[p]) for p in range(3)},'unresolved_moments':{f'M{p}':pack(R[p]) for p in range(3)},'minimum_remainder_eigenvalues':mins,'minimum_block_Hankel_eigenvalue':float(np.linalg.eigvalsh(np.block([[R[0],R[1]],[R[1].conj().T,R[2]]]))[0])}
  E0=d['E0'];Wmax=.25*(27+54*j)-E0;remain0=R[0][2,2].real;remain1=R[1][2,2].real;remain2=R[2][2,2].real
  for threshold in [.01,.04,.05,.1]:
   v=sum(mat(r['residue_matrix'])[2,2].real for r in known if r['omega']<threshold);first=max(0,remain0-remain1/threshold);second=max(0,remain0-((Wmax+threshold)*remain1-remain2)/(Wmax*threshold));lower=max(first,second);tail_lower=max(0,(remain1-threshold*remain0)/(Wmax-threshold),(remain2-threshold*remain1)/(Wmax*(Wmax-threshold)))
   if remain0>1e-15 and remain1/remain0>threshold:
    denominator=remain2-2*threshold*remain1+threshold**2*remain0
    if denominator>1e-25:tail_lower=max(tail_lower,(remain1-threshold*remain0)**2/denominator)
   upper_remaining=max(0,remain0-tail_lower);bounds.append({'Jp':j,'threshold_J':threshold,'verified_hub_weight_below':v,'verified_fraction_total':v/M[0][2,2].real,'unresolved_total_hub_weight':remain0,'additional_weight_below_lower_bound_M1':first,'additional_weight_below_lower_bound_M1_M2_bandwidth':second,'total_weight_below_lower_bound':v+lower,'total_weight_below_upper_bound':v+upper_remaining,'simple_upper_bound_from_remainder_mass':v+remain0,'conservative_total_upper_bound_plus_1e-8':min(v+remain0,v+upper_remaining+1e-8),'conservative_total_lower_bound_minus_1e-8':max(v,v+lower-1e-8),'maximum_possible_excitation_J':Wmax,'bound_note':'Positivity and moment bounds conditioned on verified numerical input, with a displayed roundoff margin; not additional poles.'})
 write('cumulative_weight_bounds.csv',bounds);(ROOT/'spectral_remainders.json').write_text(json.dumps(remainder,indent=2))
 required=[ROOT.parent/f'static/j0.2_n01_r{r}.npz' for r in [0,1]]+[ROOT.parent/f'crossing/p0.237500000_n01_r{r}_e0.npz' for r in [0,1]]+[ROOT/(d['name']+'.npz') for d in records if d['origin']=='NEW_DEFLATED_SOLVE']
 if not all(p.exists() for p in required):
  if not (ROOT/'matched_branch_overlaps.json').exists():raise FileNotFoundError('Use the full checkpoint archive for vector-overlap recomputation.')
  return rows,bounds,records
 overlaps=[]
 if ref:
  for r in [0,1]:
   a=np.load(ROOT.parent/f'static/j0.2_n01_r{r}.npz')['vector'];b=np.load(ROOT.parent/f'crossing/p0.237500000_n01_r{r}_e0.npz')['vector'];overlaps.append({'reflection_character':r,'spin':.5,'abs_overlap':float(abs(np.vdot(a,b))),'matching':'Same reflection/S labels and direct common-basis cyclic-sector eigenvector overlap; not sorted-index-only matching.'})
 from scipy.optimize import linear_sum_assignment
 labels=sorted({(d['reflection_character'],d['spin']['nearest_half_integer_spin']) for d in records if d['origin']=='NEW_DEFLATED_SOLVE'})
 for reflection,spin in labels:
  A=[d for d in records if d['Jp']==.2 and d['origin']=='NEW_DEFLATED_SOLVE' and d['reflection_character']==reflection and d['spin']['nearest_half_integer_spin']==spin];B=[d for d in records if d['Jp']==.2375 and d['origin']=='NEW_DEFLATED_SOLVE' and d['reflection_character']==reflection and d['spin']['nearest_half_integer_spin']==spin]
  if not A or not B:continue
  scores=np.array([[float(abs(np.vdot(np.load(ROOT/(x['name']+'.npz'))['vector'],np.load(ROOT/(y['name']+'.npz'))['vector']))) for y in B] for x in A]);ii,jj=linear_sum_assignment(-scores)
  for i,k in zip(ii,jj):overlaps.append({'source':A[i]['name'],'target':B[k]['name'],'reflection_character':reflection,'spin':spin,'abs_overlap':float(scores[i,k]),'status':'MATCHED' if scores[i,k]>.8 else 'TENTATIVE','matching':'Maximum-total-overlap one-to-one assignment within equal spin and reflection; not sorted-index matching.'})
 (ROOT/'matched_branch_overlaps.json').write_text(json.dumps(overlaps,indent=2));return rows,bounds,records
if __name__=='__main__':
 r,b,_=collect();print('Verified poles',len(r));print(json.dumps(b,indent=2))
