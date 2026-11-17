import numpy as np,scipy.linalg as la,itertools
from local_geometry import faces,w,I,B,D,chi,P,flavor
from ed import geometry,translations
from sectors import point_permutations

def unpack_rho(rho):
 bonds=np.array([np.trace(rho@x).real for x in B]);c=float(np.trace(rho@chi).real);twirl=I/8+2/3*sum(g*x for g,x in zip(bonds,B))+4/3*c*chi;sig=np.array([[sum(np.vdot(flavor[i,m],rho@flavor[j,m]) for m in range(2)) for j in range(2)] for i in range(2)]);p=float(np.trace(rho@P).real);normsig=sig/p;pol=np.array([2*normsig[0,1].real,-2*normsig[0,1].imag,(normsig[0,0]-normsig[1,1]).real]);return bonds,c,twirl,sig,{'P_doublet':p,'P_quartet':1-p,'chirality':c,'chirality_variance':float(np.trace(rho@chi@chi).real-c*c),'flavor_purity':float(np.trace(normsig@normsig).real),'flavor_pseudospin':pol.tolist(),'rho_purity':float(np.trace(rho@rho).real)}
def maps(p):
 f,o,_,_=faces();lookup={frozenset(t):i for i,t in enumerate(f)};idx=np.empty(162,int);coef=np.ones(162,complex);rho_maps=[]
 for t in range(54):
  j=lookup[frozenset(p[f[t]])];pos=[list(o[j]).index(v) for v in p[o[t]]];k=pos[0];proper=(pos[1]-pos[0])%3==1;phase=w**(-k) if proper else w**(k-1);idx[t]=j if proper else 54+j;coef[t]=phase;idx[54+t]=54+j if proper else j;coef[54+t]=phase.conjugate();idx[108+t]=108+j;coef[108+t]=1 if proper else -1
  slots=[list(f[j]).index(v) for v in p[f[t]]];bmap=np.array([sum(((a>>u)&1)<<slots[u] for u in range(3)) for a in range(8)]);rho_maps.append((j,bmap))
 return idx,coef,rho_maps

def realspace(Gq):
 R=[(x,y) for x in range(3) for y in range(3)];q=np.array(R)/3;out=np.zeros((162,162),complex)
 for i,a in enumerate(R):
  ii=[k*54+i*6+mu for k in range(3) for mu in range(6)]
  for j,b in enumerate(R):
   jj=[k*54+j*6+mu for k in range(3) for mu in range(6)];weights=np.exp(2j*np.pi*(q@np.subtract(b,a)))/9;out[np.ix_(ii,jj)]=np.einsum('q,qab->ab',weights,Gq)
 return out

def ensemble(raw,means,rhos,j):
 tridx=np.r_[np.arange(54,108),np.arange(54),np.arange(108,162)];trsign=np.r_[np.ones(108),-np.ones(54)];base=(raw+raw[np.ix_(tridx,tridx)].conj()*trsign[:,None]*trsign[None,:])/2;base_m=(means+means[tridx].conj()*trsign)/2;twirl=np.array([unpack_rho(r)[2] for r in rhos]);twirl=(twirl+twirl.conj())/2
 group=list(point_permutations())
 if j==1:
  coord,_=geometry(3);_,ts=translations(3,coord,True);group=[p[t] for p in group for t in ts[:3]]
 out=np.zeros_like(raw);mout=np.zeros_like(means);rout=np.zeros((6,8,8),complex)
 for p in group:
  ix,c,rm=maps(p);out+=c.conj()[:,None]*c[None,:]*base[np.ix_(ix,ix)];mout+=c*base_m[ix]
  for mu in range(6):
   t,b=rm[mu];rout[mu]+=twirl[t%6][np.ix_(b,b)]
 out/=len(group);mout/=len(group);rout/=len(group);return out,mout,rout

def fourier(C):
 R=np.array([(x,y) for x in range(3) for y in range(3)]);out=[]
 for n in R:
  V=np.zeros((54,6),complex)
  for r in range(9):
   for mu in range(6):V[r*6+mu,mu]=np.exp(-2j*np.pi*(R[r]@n)/3)/3
  out.append(V.conj().T@C@V)
 return np.array(out)

def unfold(S,kind):
 _,_,centers,vertices=faces();c=centers[:6];geo=np.ones(6,complex)
 if kind=='D':
  for mu,pts in enumerate(vertices[:6]):
   order=sorted(range(3),key=lambda i:(round(pts[i,1],10),round(pts[i,0]+pts[i,1]/2,10)));geo[mu]=w**(-order[0])
 rows=[];off=0.;unit=0.
 for qid,(n1,n2) in enumerate(itertools.product(range(3),repeat=2)):
  base=np.array([(n1+n2)/9,(2*n1-n2)/9]);W=np.zeros((6,6),complex);Qs=[]
  for fold in range(3):
   Q=np.mod(base+fold*np.array([1/3,2/3]),1);Qs.append(Q)
   for mu in range(6):W[mu,2*fold+mu%2]=geo[mu]*np.exp(-2j*np.pi*(c[mu]@Q))/np.sqrt(3)
  unit=max(unit,float(la.norm(W.conj().T@W-np.eye(6))));U=W.conj().T@S[qid]@W
  for f,Q in enumerate(Qs):
   block=U[2*f:2*f+2,2*f:2*f+2];e,v=la.eigh(block);rows.append({'Q':Q,'cell_n':[n1,n2],'fold':f,'matrix':block,'eigenvalues':e,'eigenvectors':v})
   for g in range(3):
    if g!=f:off=max(off,float(la.norm(U[2*f:2*f+2,2*g:2*g+2])))
 return rows,{'unfolding_unitarity_error':unit,'different_primitive_momentum_block_error':off,'geometric_dimer_gauge':geo}

def mixture_decomposition(raw,means,j):
 tridx=np.r_[np.arange(54,108),np.arange(54),np.arange(108,162)];sign=np.r_[np.ones(108),-np.ones(54)];partner=means[tridx].conj()*sign;moment=(means.conj()[:,None]*means[None,:]+partner.conj()[:,None]*partner[None,:])/2;group=list(point_permutations())
 if j==1:
  coords,_=geometry(3);_,ts=translations(3,coords,True);group=[p[t] for p in group for t in ts[:3]]
 second=np.zeros((162,162),complex)
 for p in group:
  ix,c,_=maps(p);second+=c.conj()[:,None]*c[None,:]*moment[np.ix_(ix,ix)]
 second/=len(group);return second
