import numpy as np,itertools,json,pathlib
from ed import geometry
ROOT=pathlib.Path(__file__).resolve().parent
w=np.exp(2j*np.pi/3);I=np.eye(8,dtype=complex);S=np.zeros((3,3,8,8),complex)
for a in range(3):
 for b in range(8):
  S[a,2,b,b]=((b>>a)&1)-.5;S[a,0,b^(1<<a),b]=.5;S[a,1,b^(1<<a),b]=.5j if (b>>a)&1 else -.5j
B=[sum(S[a,k]@S[b,k] for k in range(3)) for a,b in [(0,1),(1,2),(2,0)]];D=B[0]+w*B[1]+w*w*B[2];chi=sum((1 if t in [(0,1,2),(1,2,0),(2,0,1)] else -1)*S[0,t[0]]@S[1,t[1]]@S[2,t[2]] for t in itertools.permutations(range(3)));P=.5*I-2/3*sum(B)
flavor=np.zeros((2,2,8),complex);flavor[0,0,[5,6]]=[1,-1];flavor[0,0]/=np.sqrt(2);flavor[0,1,[1,2]]=[1,-1];flavor[0,1]/=np.sqrt(2);flavor[1,0,3]=np.sqrt(2/3);flavor[1,0,[5,6]]=-1/np.sqrt(6);flavor[1,1,4]=-np.sqrt(2/3);flavor[1,1,[1,2]]=1/np.sqrt(6)
def sid(x,y,s,L=3):return 3*((x%L)*L+(y%L))+s
def faces(L=3):
 f=[];oriented=[];centers=[];verts=[]
 offsets=[(0,-1),(1,-1),(1,0),(0,1),(-1,1),(-1,0)]
 for x in range(L):
  for y in range(L):
   neighbors=[sid(x,y,0,L),sid(x,y,1,L),sid(x+1,y,0,L),sid(x+1,y-1,1,L),sid(x+1,y-1,0,L),sid(x,y-1,1,L)];hub=sid(x,y,2,L);C=np.array([x+2*y,x-y+1.])
   for mu in range(6):
    u,v=neighbors[mu],neighbors[(mu+1)%6];a,b=(u,v) if u%3==0 else (v,u);f.append([a,b,hub]);order=[a,b,hub] if mu%2==0 else [a,hub,b];oriented.append(order);pts={hub:C,u:C+offsets[mu],v:C+offsets[(mu+1)%6]};verts.append(np.array([pts[i] for i in order]));centers.append(sum(pts.values())/3)
 return np.array(f,np.int32),np.array(oriented,np.int32),np.array(centers),np.array(verts)
def matrices():
 ops=[]
 for mu in range(6):
  d=D if mu%2==0 else w*w*D.conj().T;c=chi if mu%2==0 else -chi;ops.append([d,d.conj().T,c])
 return np.ascontiguousarray(np.array(ops))
def verify():
 checks={'projector':np.linalg.norm(P@P-P),'chi_squared':np.linalg.norm(chi@chi-3*P/16),'dimer_symmetric':np.linalg.norm((D.conj().T@D+D@D.conj().T)/2-9*P/8),'dimer_commutator':np.linalg.norm(D.conj().T@D-D@D.conj().T-3*np.sqrt(3)*chi),'flavor_orthonormality':np.linalg.norm(flavor.reshape(4,8)@flavor.reshape(4,8).conj().T-np.eye(4))};assert max(checks.values())<2e-14
 f,o,c,v=faces();coords,bonds=geometry(3);edges={tuple(sorted((a,b))) for a,b,J in bonds};assert len(edges)==81 and len({tuple(sorted(t)) for t in f})==54 and all(all(tuple(sorted(e)) in edges for e in itertools.combinations(t,2)) for t in f)
 xy=np.stack([v[:,:,0]+v[:,:,1]/2,np.sqrt(3)*v[:,:,1]/2],axis=-1);a=xy[:,1]-xy[:,0];b=xy[:,2]-xy[:,0];areas=a[:,0]*b[:,1]-a[:,1]*b[:,0];assert min(areas)>0
 counts={k:sum(len(set(a)&set(b))==k for a in f for b in f) for k in range(4)};checks['scalar_spin_commutator_error']=max(float(np.linalg.norm(op@sum(S[:,a])-sum(S[:,a])@op)) for op in [D,chi] for a in range(3));assert checks['scalar_spin_commutator_error']<1e-13
 pairing=np.zeros((3,8),complex);pairing[0,[5,6]]=[1,-1];pairing[1,[3,5]]=[1,-1];pairing[2,[6,3]]=[1,-1];pairing/=np.sqrt(2);checks['nonorthogonal_pairing_Gram']=np.real(pairing@pairing.conj().T).tolist();checks['pairing_Gram_eigenvalues']=np.linalg.eigvalsh(pairing@pairing.conj().T).tolist()
 checks['AB_bonds']=sum({a%3,b%3}=={0,1} for a,b,J in bonds);checks['hub_bonds']=sum(2 in [a%3,b%3] for a,b,J in bonds);degree=np.bincount(np.array([(a,b) for a,b,J in bonds]).ravel(),minlength=27);assert min(degree)==max(degree)==6 and checks['AB_bonds']==27 and checks['hub_bonds']==54;checks['site_degree']=degree.tolist();gamma=lambda q:2*(np.cos(2*np.pi*q[0])+np.cos(2*np.pi*q[1])+np.cos(2*np.pi*(q[0]-q[1])));checks['NN_structure_factor_at_physical_K']=float(gamma([1/3,2/3]));checks['NN_structure_factor_at_wrong_K_11']=float(gamma([1/3,1/3]));assert abs(checks['NN_structure_factor_at_physical_K']+3)<1e-12
 checks.update(faces=len(f),bonds=len(edges),ordered_pair_overlap_counts=counts,minimum_twice_oriented_area=float(min(areas)));(ROOT/'local_geometry_validation.json').write_text(json.dumps(checks,indent=2));return checks
if __name__=='__main__':print(verify())
