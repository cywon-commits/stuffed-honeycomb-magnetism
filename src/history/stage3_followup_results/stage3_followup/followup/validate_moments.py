import pathlib,json,numpy as np
from ed import geometry
from observables import coefficients,site_weights
from search import pack
ROOT=pathlib.Path(__file__).resolve().parent
def unpack(x):return np.array(x['real'])+1j*np.array(x['imag'])
def run():
 rows=[];out={}
 for j in [.2,.2375]:
  if j==.2:
   file=ROOT/'matched_j0.2.json';real=ROOT.parent/'static/static_j0.2_realspace.npz'
  else:file=ROOT.parent/'dynamics/weak_moments.json';real=ROOT.parent/'crossing/obs_p0.237500000_n12_r0_e0_realspace.npz'
  if not file.exists():continue
  d=json.loads(file.read_text());z=np.load(real);G=z['Gdot'] if 'Gdot' in z else z['G_selected'];Gz=z['Gzz'] if 'Gzz' in z else z['Gzz_selected'];_,bonds=geometry(3,j);out[str(j)]={}
  for sign in [1,-1]:
   k=[0,sign/3];cf=coefficients(k,'Phi');C=np.column_stack([np.eye(3),[(cf[0]-cf[1])/2,(cf[0]+cf[1])/2,0]]);M=[C.conj().T@unpack(d['moments'][str(sign)][f'M{p}'])@C for p in range(3)];weights=np.array([site_weights(3,k,coefficients(k,ch)) for ch in ['n','m','H','Phi']]);static=weights.conj()@Gz@weights.T;fsum=np.zeros((4,4),complex)
   for a,b,J in bonds:
    diff=weights[:,a]-weights[:,b];fsum-=.5*J*(G[a,b]-Gz[a,b])*np.outer(diff.conj(),diff)
   # Symmetrized exact first moment requires reverse transpose, then basis transformation at forward k.
   f=(unpack(d['moments'][str(sign)]['M1'])+unpack(d['moments'][str(-sign)]['M1']).T)/2;f=C.conj().T@f@C
   e0=float(np.max(abs(M[0]-static)));e1=float(np.max(abs(f-fsum)));assert max(e0,e1)<1e-8
   out[str(j)][str(sign)]={f'M{p}':pack(M[p]) for p in range(3)}
   rows.append({'Jp':j,'probe_sign':sign,'M0_static_error':e0,'M1_bond_fsum_error':e1,'moment_Hermiticity_errors':[float(np.max(abs(x-x.conj().T))) for x in M],'minimum_moment_eigenvalues':[float(np.linalg.eigvalsh((x+x.conj().T)/2)[0]) for x in M]})
 (ROOT/'four_channel_exact_moments.json').write_text(json.dumps(out,indent=2));(ROOT/'moment_validation.json').write_text(json.dumps(rows,indent=2));print(rows)
if __name__=='__main__':run()
