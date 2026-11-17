"""Independent full-bit-state checks of reduced spin seeds and static matrices."""
import numpy as np,scipy.linalg as la,pathlib,json
from itertools import combinations
from ed import Model,geometry
from observables import *
from sectors import cell_stars
ROOT=pathlib.Path(__file__).resolve().parent

def states(n,up):
 ss=sorted(sum(1<<p for p in ps) for ps in combinations(range(n),up));return ss,{s:i for i,s in enumerate(ss)}
def projector(m):
 ss,idx=states(m.n,m.up);U=np.zeros((len(ss),m.nb),complex)
 for j,r in enumerate(m.masks):
  for g,p in enumerate(m.perms):
   s=sum(1<<int(p[a]) for a in range(m.n) if int(r)>>a&1);U[idx[s],j]=m.phase[g]/np.sqrt(float(m.lengths[j]))
 assert la.norm(U.conj().T@U-np.eye(m.nb))<1e-12
 return U,ss,idx

def direct_seed(n,up,psi,weights,kind):
 src,_=states(n,up);tu=up+{'z':0,'+':1,'-':-1}[kind];tar,idx=states(n,tu);out=np.zeros(len(tar),complex)
 for i,s in enumerate(src):
  for p,w in enumerate(weights):
   bit=(s>>p)&1
   if kind=='z':out[idx[s]]+=w*(bit-.5)*psi[i]
   elif (kind=='+' and not bit) or (kind=='-' and bit):out[idx[s^(1<<p)]]+=w*psi[i]
 return out

def direct_pairs(n,up,psi):
 ss,idx=states(n,up);G=np.eye(n)*.75;Gz=np.eye(n)*.25;mag=np.zeros(3)
 for z,s in enumerate(ss):
  for p in range(n):mag[p%3]+=abs(psi[z])**2*(((s>>p)&1)-.5)/(n/3)
 for a in range(n):
  for b in range(a+1,n):
   d=0j;zz=0.
   for i,s in enumerate(ss):
    anti=((s>>a)^(s>>b))&1;val=abs(psi[i])**2*(.25-.5*anti);d+=val;zz+=val
    if anti:d+=.5*np.conj(psi[i])*psi[idx[s^(1<<a)^(1<<b)]]
   assert abs(d.imag)<1e-12;G[a,b]=G[b,a]=d.real;Gz[a,b]=Gz[b,a]=zz
 return G,Gz,mag

def direct_h_action(n,up,psi,bonds):
 ss,idx=states(n,up);out=np.zeros_like(psi)
 for i,s in enumerate(ss):
  for a,b,j in bonds:
   anti=((s>>a)^(s>>b))&1;out[i]+=j*(.25-.5*anti)*psi[i]
   if anti:out[idx[s^(1<<a)^(1<<b)]]+=.5*j*psi[i]
 return out

def run():
 records=[];rng=np.random.default_rng(195)
 for L,up,ks in [(2,6,(0,0)),(2,5,(.5,0)),(3,1,(1/3,0))]:
  m=Model(L,up,False,.73);m.sector(cell_q(ks));U,ss,idx=projector(m);v=rng.normal(size=m.nb)+1j*rng.normal(size=m.nb);v/=la.norm(v);psi=U@v
  stat=pair_expectations(m,v);G,Gz,mag=direct_pairs(m.n,up,psi)
  err=max(la.norm(G-stat['Gdot']),la.norm(Gz-stat['Gzz']),la.norm(mag-stat['magnetization']));assert err<1e-10
  C=correlation_matrices(stat['Gdot'],L);worst=0.;seederr=0.;hseederr=0.
  targets={kind:Model(L,up+shift,False,.73) for kind,shift in [('z',0),('+',1),('-',-1)]}
  for a in range(L):
   for b in range(L):
    kt=np.array([a,b])/L;expected=np.zeros((3,3),complex)
    for kind,factor in [('z',1.),('+',.5),('-',.5)]:
     tar=targets[kind];tar.sector(cell_q(np.array(ks)+kt));Ut,_,_=projector(tar);seeds=[]
     for alpha in range(3):
      coeff=np.eye(3)[alpha];w=np.array([coeff[s]*np.exp(-2j*np.pi*(kt[0]*x+kt[1]*y))/L for x in range(L) for y in range(L) for s in range(3)])
      full=direct_seed(m.n,up,psi,w,kind);reduced=spin_seed(m,tar,v,kt,['A','B','H'][alpha],kind)
      er=la.norm(Ut@reduced-full);seederr=max(seederr,er);assert er<1e-10
      er=la.norm(Ut@tar.mv(reduced)-direct_h_action(m.n,tar.up,full,tar.bonds));hseederr=max(hseederr,er);assert er<1e-10
      seeds.append(full)
     expected+=factor*np.array([[np.vdot(x,y) for y in seeds] for x in seeds])
    worst=max(worst,la.norm(C[(a,b)]-expected));assert la.norm(C[(a,b)]-expected)<1e-10
    # Derive the Phi weights directly from the six-neighbor definition, not its form factors.
    w=np.zeros(m.n,complex)
    def sid(x,y,s):return 3*((x%L)*L+y%L)+s
    for x in range(L):
     for y in range(L):
      phase=np.exp(-2j*np.pi*(kt[0]*x+kt[1]*y))/L
      for xx,yy,s in [(x,y,0),(x,y,1),(x+1,y,0),(x+1,y-1,1),(x+1,y-1,0),(x,y-1,1)]:w[sid(xx,yy,s)]+=phase
    assert la.norm(w-site_weights(L,kt,coefficients(kt,'Phi')))<1e-12
  records.append({'L':L,'up':up,'source_cell_k':ks,'pair_matrix_error':float(err),'spin_seed_max_error':float(seederr),'H_seed_max_error':float(hseederr),'full_C_max_error':float(worst)})
 # Primitive K basis conversion and physical-site phases, in a small one-magnon sector.
 prim=Model(3,1,True);prim.sector([1/3,2/3]);v=rng.normal(size=prim.nb)+1j*rng.normal(size=prim.nb);v/=la.norm(v)
 cellm=Model(3,1,False);cellm.sector([0,0]);vc=convert_basis(prim,cellm,v);Up,_,_=projector(prim);Uc,_,_=projector(cellm)
 er=la.norm(Up@v-Uc@vc);assert er<1e-12
 physical=np.exp(-2j*np.pi*(prim.coords@np.array([1/3,2/3])))/np.sqrt(27)
 assert la.norm(physical-site_weights(3,[0,0],coefficients([0,0],'physical_K')))<1e-12
 # Centered A and B rim displacement sums vanish individually; staggered Taylor term is cubic.
 dA=np.array([[0,-1],[1,0],[-1,1]]);dB=-dA
 assert np.all(dA.sum(axis=0)==0) and np.all(dB.sum(axis=0)==0)
 q=np.array([.7,.2]);assert abs(np.sum((dA@q)**3)-np.sum((dB@q)**3))>1e-3
 records.append({'primitive_to_cell_error':float(er),'physical_K_weights_verified':True,'staggered_centered_form_factor_linear_term_zero':True,'leading_generic_odd_term_cubic':True})
 stars,_=cell_stars(jp=.2);assert sorted(map(len,stars))==[1,2,6]
 (ROOT/'operator_validation.json').write_text(json.dumps(records,indent=2));print(json.dumps(records,indent=2));print('ALL OBSERVABLE TESTS PASSED')
if __name__=='__main__':run()
