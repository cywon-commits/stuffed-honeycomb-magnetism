"""Basis-independent degenerate residues using the measured eigenvector Gram matrix.
For nonorthogonal eigenvectors V, spectral projector is V(V†V)^-1 V†.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,numpy as np,scipy.linalg as la
from experiment import ROOT,pack,dump

def unpack(x):return np.array(x['real'])+1j*np.array(x['imag'])
def run():
 allgroups=[];diags=[];small={};oldsmall=np.load(ROOT/'pole_subspace_data.npz') if (ROOT/'pole_subspace_data.npz').exists() else None
 for f in sorted(ROOT.glob('*_spectrum.json')):
  rec=json.loads(f.read_text());tag=rec.get('tag',rec['name']);path=ROOT/(tag+'_poles.npz')
  if path.exists():
   data=np.load(path);A=data['amplitudes']
   if len(A)!=len(rec['poles']):continue
   cache=ROOT/(tag+'_eigenvector_Gram.npz')
   if cache.exists() and np.load(cache)['Gram'].shape==(len(A),len(A)):G=np.load(cache)['Gram']
   else:
    V=data['vectors'];G=la.blas.zgemm(1.,V,V,trans_a=2);np.savez_compressed(cache,Gram=G);del V
  elif oldsmall is not None and tag+'_Gram' in oldsmall:
   G=oldsmall[tag+'_Gram'];A=oldsmall[tag+'_amplitudes'];assert len(A)==len(rec['poles'])
  else:continue
  small[tag+'_Gram']=G;small[tag+'_amplitudes']=A;groups=[];used=set()
  for i,p in enumerate(rec['poles']):
   if i in used:continue
   ix=[k for k,q in enumerate(rec['poles']) if k not in used and abs(q['energy']-p['energy'])<1e-8 and q['spin']['nearest_half_integer_spin']==p['spin']['nearest_half_integer_spin']];used.update(ix);GG=G[np.ix_(ix,ix)];ev=la.eigvalsh(GG);assert ev.min()>1e-8;R=A[ix].conj().T@la.solve(GG,A[ix],assume_a='her');R=(R+R.conj().T)/2;ok=all(rec['poles'][i]['verified_eigenpair'] for i in ix);scalar=all(rec['poles'][i]['scalar_accessible'] for i in ix);row={'source_file':f.name,'tag':tag,'name':rec['name'],'indices':ix,'multiplicity':len(ix),'energy':float(np.mean([rec['poles'][i]['energy'] for i in ix])),'omega':float(np.mean([rec['poles'][i]['omega'] for i in ix])),'spin':p['spin']['nearest_half_integer_spin'],'verified_eigenpair':ok,'scalar_accessible':scalar,'max_direct_H_residual':max(rec['poles'][i]['H_residual'] for i in ix),'max_S2_residual':max(rec['poles'][i]['spin']['S2_eigen_residual'] for i in ix),'angular':rec.get('angular'),'reflection':rec.get('reflection'),'residue':pack(R),'residue_trace':float(np.trace(R).real),'Gram_eigenvalues':ev.tolist(),'status':'VERIFIED_DEGENERATE_SCALAR_RESIDUE' if ok and scalar else ('VERIFIED_SPIN_FORBIDDEN_SUBSPACE' if ok else 'UNCONVERGED_NOT_PHYSICAL_POLE')};groups.append(row)
  # Authoritative corrections to moment sums; never sum nonorthogonal vectors naively.
  nc=A.shape[1];sums=[sum((g['omega']**p*unpack(g['residue']) for g in groups if g['verified_eigenpair'] and g['scalar_accessible']),np.zeros((nc,nc),complex)) for p in range(3)];mom=json.loads((ROOT/(rec['name']+'_moments.json')).read_text());rec['partial_moments']=[pack(x) for x in sums];rec['remainder_moments']=[pack(unpack(mom[f'M{p}'])-sums[p]) for p in range(3)];rec['residue_sum_convention']='Authoritative sums use degenerate-subspace V(V†V)^-1 V†. Individual solver vectors may be nonorthogonal inside exact degeneracies; individual residue matrices are basis-dependent.';dump(f,rec);allgroups.extend(groups);diags.append({'tag':tag,'maximum_Gram_offdiagonal':float(np.max(abs(G-np.eye(len(G))))),'subspace_groups':len(groups)})
 if oldsmall is not None:oldsmall.close()
 np.savez_compressed(ROOT/'pole_subspace_data.npz',**small);dump(ROOT/'degenerate_pole_groups.json',{'groups':allgroups,'orthogonality_diagnostics':diags});print('Grouped',len(allgroups),'eigenspaces from',len(diags),'solves')
if __name__=='__main__':run()
