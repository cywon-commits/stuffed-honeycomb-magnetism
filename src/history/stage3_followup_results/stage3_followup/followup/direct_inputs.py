"""Direct bit-diagonal Fourier seeds with exact target-reflection projection."""
import pathlib,json,numpy as np,scipy.linalg as la
from cyclic_model import CyclicModel
from observables import full_samples,site_weights,coefficients
from sectors import point_permutations
ROOT=pathlib.Path(__file__).resolve().parent
def preimage(m,perm):
 lo=np.zeros(1<<m.split,np.uint32);hi=np.zeros(1<<(m.n-m.split),np.uint32)
 for table,offset in [(lo,0),(hi,m.split)]:
  for mask in range(1,len(table)):
   bit=mask&-mask;table[mask]=table[mask^bit]|(1<<int(perm[offset+bit.bit_length()-1]))
 return lo[m.masks&((1<<m.split)-1)]|hi[m.masks>>m.split]
def inputs(j,r):
 data=ROOT.parent/('static' if j==.2 else 'crossing');srcname='j0.2_n12_r0' if j==.2 else 'p0.237500000_n12_r0_e0';knownname=f'j0.2_n01_r{r}' if j==.2 else f'p0.237500000_n01_r{r}_e0';z=np.load(data/(srcname+'.npz'));src=CyclicModel(jp=j,cell_n=(1,2),angular=0);m=CyclicModel(jp=j,cell_n=(1,0),angular=r);pre=preimage(m,point_permutations()[7]);psi=full_samples(src,z['vector'],m.masks);psif=full_samples(src,z['vector'],pre);B=[]
 for ch in ['n','m','H']:
  w=site_weights(3,[0,1/3],coefficients([0,1/3],ch));d=sum(w[a]*(((m.masks>>a)&1)-.5) for a in range(27));df=sum(w[a]*(((pre>>a)&1)-.5) for a in range(27));B.append(.5*(d*psi+(-1)**r*df*psif)*np.sqrt(m.lengths.astype(float)))
 E0=float(z['energy']);del src;known=np.load(data/(knownname+'.npz'));old=CyclicModel(jp=j,cell_n=(0,1),angular=r);gs=m.gs[:9];perms=m.perms[:9];oldchars=np.exp(-2j*np.pi*(gs@old.q));newchars=np.exp(-2j*np.pi*(gs@m.q));lookup={tuple(p):i for i,p in enumerate(perms)}
 for p in point_permutations():
  inv=np.argsort(p)
  if all(abs(oldchars[lookup[tuple(inv[t[p]])]]-newchars[i])<1e-8 for i,t in enumerate(perms)):break
 else:raise RuntimeError('No star map')
 q=full_samples(old,known['vector'],preimage(m,inv))*np.sqrt(m.lengths.astype(float));assert abs(la.norm(q)-1)<1e-9;B=np.column_stack(B);amp=q.conj()@B;R=np.outer(amp.conj(),amp)
 if j==.2:ref=json.loads((ROOT/'matched_j0.2.json').read_text())['poles'][r]
 else:ref=json.loads((ROOT.parent/'dynamics/weak_spectrum.json').read_text())['poles'][r]
 refR=np.array(ref['residue_matrix']['real'])+1j*np.array(ref['residue_matrix']['imag']);err=float(la.norm(R-refR));assert err<1e-9
 (ROOT/f'direct_seed_j{j:g}_r{r}_validation.json').write_text(json.dumps({'independent_bit_diagonal_projected_seed_known_residue_error':err,'known_norm_error':float(abs(la.norm(q)-1)),'status':'PASS'},indent=2));return m,B,q,E0
