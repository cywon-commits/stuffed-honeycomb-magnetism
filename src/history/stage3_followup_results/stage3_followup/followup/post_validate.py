import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import pathlib,json,numpy as np,scipy.linalg as la
from direct_inputs import inputs
from observables import total_spin_check,inversion_vector
from sectors import point_permutations
from export import collect,mat
ROOT=pathlib.Path(__file__).resolve().parent
def run():
 out=[];groups=[];direct_moments={j:[np.zeros((3,3),complex) for _ in range(3)] for j in [.2,.2375]}
 for j in [.2375,.2]:
  for r in [0,1]:
   files=[]
   for p in sorted(ROOT.glob(f'j{j:g}_r{r}_root*.json')):
    d=json.loads(p.read_text())
    if d.get('verified'):files.append((p,d))
   if not files:continue
   m,B,q,E0=inputs(j,r);Q=[q];HB=np.column_stack([m.mv(B[:,i])-E0*B[:,i] for i in range(3)])
   for order,M in enumerate([B.conj().T@B,B.conj().T@HB,HB.conj().T@HB]):direct_moments[j][order]+=M
   for p,d in files:
    u=np.load(p.with_suffix('.npz'))['vector'];old=u.copy()
    for repeat in range(2):
     for q in Q:u-=q*np.vdot(q,u)
    u/=la.norm(u);change=float(la.norm(u-old));res=float(la.norm(m.mv(u)-d['energy']*u));spin=total_spin_check(m,u);a=u.conj()@B;original=np.array(d['amplitudes']['real'])+1j*np.array(d['amplitudes']['imag']);aerr=float(la.norm(a-original));fr=float(la.norm(inversion_vector(m,u,point_permutations()[7])-(-1)**r*u));assert change<1e-8 and res<1e-8 and spin['S2_eigen_residual']<1e-8 and aerr<1e-9 and fr<1e-8;Q.append(u);out.append({'name':d['name'],'explicit_two_pass_reorthogonalization_change':change,'post_orthogonalization_H_residual':res,'S2_residual':spin['S2_eigen_residual'],'independent_bit_seed_amplitude_error':aerr,'reflection_eigen_residual':fr})
   Gram=np.array([[np.vdot(a,b) for b in Q] for a in Q]);err=float(np.max(abs(Gram-np.eye(len(Q)))));assert err<1e-9;groups.append({'Jp':j,'reflection_character':r,'physical_vectors_including_known':len(Q),'Gram_matrix_error':err})
 rows,bounds,records=collect();rem=json.loads((ROOT/'spectral_remainders.json').read_text())
 for d in rem.values():assert min(d['minimum_remainder_eigenvalues'])>-1e-8 and d['minimum_block_Hankel_eigenvalue']>-1e-8
 for b in bounds:assert b['total_weight_below_lower_bound']<=b['total_weight_below_upper_bound']+1e-8
 small=json.loads((ROOT/'small_validation.json').read_text());assert small['status']=='PASS'
 moment_errors={}
 for j in [.2,.2375]:
  f=ROOT/'matched_j0.2.json' if j==.2 else ROOT.parent/'dynamics/weak_moments.json';reference=json.loads(f.read_text())['moments']['1'];errors=[float(la.norm(direct_moments[j][i]-mat(reference[f'M{i}']))) for i in range(3)];assert max(errors)<1e-8;moment_errors[str(j)]=errors
 result={'independent_reflection_sum_M0_M1_M2_matrix_norm_errors':moment_errors,'status':'PASS','post_checked_new_states':out,'orthogonality_groups':groups,'all_verified_pole_count':len(rows),'verified_additional_2375_pole_count':sum(d['Jp']==.2375 and d['origin']=='NEW_DEFLATED_SOLVE' for d in rows),'pole_and_remainder_PSD':True,'matrix_moment_block_Hankel_PSD':True,'max_H_residual_new_after_reorth':max(x['post_orthogonalization_H_residual'] for x in out),'max_S2_residual_new':max(x['S2_residual'] for x in out),'max_independent_seed_amplitude_error':max(x['independent_bit_seed_amplitude_error'] for x in out)}
 (ROOT/'final_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':run()
