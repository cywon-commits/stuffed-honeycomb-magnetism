import numpy as np,itertools,json,pathlib,csv
ROOT=pathlib.Path(__file__).resolve().parent
METRIC=np.array([[1,.5],[.5,1.]])
def key(T,x):
 a,b=T[0];c,d=T[1];det=abs(int(round(np.linalg.det(T))));return ((int(x[0])*d-int(x[1])*c)%det,(-int(x[0])*b+int(x[1])*a)%det)
def cells(T):
 N=abs(int(round(np.linalg.det(T))));R=[];lookup={};todo=[(0,0)]
 while todo:
  p=todo.pop(0);k=key(T,p)
  if k in lookup:continue
  lookup[k]=len(R);R.append(p)
  for e in [(1,0),(0,1)]:todo.append(tuple(np.array(p)+e))
  if len(R)==N:break
 assert len(R)==N;return np.array(R),lookup
def momenta(T):
 N=abs(int(round(np.linalg.det(T))));out=[]
 for a,b in itertools.product(range(N),repeat=2):
  q=np.array([a,b])/N
  if np.max(abs(np.exp(2j*np.pi*(T@q))-1))<1e-10:out.append(q)
 assert len(out)==N;return np.array(out)
def make(T):
 T=np.array(T,int);R,L=cells(T);N=len(R)
 def sid(x,y,s):return 3*L[key(T,(x,y))]+s
 bonds=[]
 for x,y in R:
  A=sid(x,y,0);bonds.extend((A,sid(xx,yy,1),0) for xx,yy in [(x,y),(x-1,y),(x,y-1)])
  bonds.extend((sid(x,y,2),sid(xx,yy,s),1) for xx,yy,s in [(x,y,0),(x,y,1),(x+1,y,0),(x+1,y-1,1),(x+1,y-1,0),(x,y-1,1)])
 edges=[tuple(sorted((a,b))) for a,b,t in bonds];unique=len(set(edges))==len(edges);perms=np.array([[sid(x+u,y+v,s) for x,y in R for s in range(3)] for u,v in R],int)
 for p in perms:assert {tuple(sorted((p[a],p[b]))) for a,b in edges}==set(edges)
 degree=np.bincount(np.array(edges).ravel(),minlength=3*N);hubrim=[sum(a==i or b==i for a,b,t in bonds if t==1) for i in range(3*N)];rimdeg=[sum(a==i or b==i for a,b,t in bonds if t==0) for i in range(3*N)];coord=np.array([(x+2*y+(s==1),x-y+(s==2)) for x,y in R for s in range(3)])
 red=T.copy()
 for _ in range(30):
  if red[1]@METRIC@red[1]<red[0]@METRIC@red[0]-1e-10:red=red[::-1].copy()
  m=int(np.rint((red[0]@METRIC@red[1])/(red[0]@METRIC@red[0])))
  if m==0:break
  red[1]-=m*red[0]
 length=np.sqrt(np.einsum('ai,ij,aj->a',red,METRIC,red)*3);angle=np.degrees(np.arccos(np.clip(red[0]@METRIC@red[1]/np.prod(length/np.sqrt(3)),-1,1)))
 shortest=min(np.sqrt((u*T[0]+v*T[1])@METRIC@(u*T[0]+v*T[1])*3) for u,v in itertools.product(range(-8,9),repeat=2) if(u,v)!=(0,0));sv=np.linalg.svd(red@np.array([[np.sqrt(3),0],[np.sqrt(3)/2,1.5]]),compute_uv=False);q=momenta(T);valid=unique and min(degree)==max(degree)==6 and all(hubrim[i]==(6 if i%3==2 else 3) for i in range(3*N)) and all(rimdeg[i]==(0 if i%3==2 else 3) for i in range(3*N))
 return {'T':T,'cells':R,'perms':perms,'bonds':bonds,'coords':coord,'momenta':q,'Ncell':N,'N':3*N,'valid':bool(valid),'Kcell_compatible':bool(np.max(abs(np.exp(2j*np.pi*(T@np.array([1/3,2/3])))-1))<1e-10),'unique_bonds':len(set(edges)),'AB_bonds':sum(t==0 for a,b,t in bonds),'hub_bonds':sum(t==1 for a,b,t in bonds),'shortest_winding_site_units':float(shortest),'reduced_T':red,'reduced_lengths':length,'reduced_angle_degrees':float(angle),'reduced_basis_aspect':float(max(length)/min(length)),'condition_number':float(sv[0]/sv[-1]),'degree':degree}
def equivalence(T,U):
 R=np.array([[0,-1],[1,1]]);F=np.array([[0,1],[1,0]])
 for refl in [np.eye(2,dtype=int),F]:
  p=refl.copy()
  for n in range(6):
   V=np.asarray(T)@p.T;M=V@np.linalg.inv(U)
   if np.max(abs(M-np.rint(M)))<1e-9:return True
   p=R@p
 return False
def catalog():
 gs=[make([[a,b],[0,6//a]]) for a in [1,2,3,6] for b in range(6//a)];k=[g for g in gs if g['valid'] and g['Kcell_compatible']];selected=make([[3,0],[2,2]]);assert k
 control=min([g for g in gs if g['valid'] and not g['Kcell_compatible']],key=lambda g:(-g['shortest_winding_site_units'],g['condition_number']))
 chosen={'K18':selected,'C18':control,'N27':make([[3,0],[0,3]]),'N36':make([[2,2],[-2,4]])};rows=[]
 for name,g in list(chosen.items())+[(f'HNF{i}',g) for i,g in enumerate(gs)]:
  rows.append({'name':name,**{x:(json.dumps(g[x].tolist()) if isinstance(g[x],np.ndarray) else g[x]) for x in ['T','Ncell','N','valid','Kcell_compatible','unique_bonds','AB_bonds','hub_bonds','shortest_winding_site_units','reduced_T','reduced_lengths','reduced_angle_degrees','reduced_basis_aspect','condition_number','momenta']}})
 for row in rows:
  T=np.array(json.loads(row['T']));P=T@np.array([[1,1],[2,-1]])
  allowed=lambda Q:bool(np.max(abs(np.exp(2j*np.pi*(P@np.asarray(Q)))-1))<1e-10)
  row['physical_triangular_K_allowed_at_r1']=allowed([1/3,2/3]);row['physical_triangular_Kprime_allowed_at_r1']=allowed([2/3,1/3]);row['exact_physical_M_points_at_r1']=json.dumps([q for q in [[.5,0],[0,.5],[.5,.5]] if allowed(q)])
 with (ROOT/'geometry_manifest.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
 proof={'number_of_index6_sublattices':len(gs),'valid_K_compatible_shapes':len(k),'all_K_compatible_shapes_pointgroup_equivalent':all(equivalence(selected['T'],g['T']) for g in k),'control_T':control['T'].tolist(),'geometry_metric':'Cell vectors (1,1),(2,-1) in primitive triangular coordinates; 60-degree metric multiplied by 3.','caution':'K-compatible here means weak-coupling CELL k=(1/3,2/3); physical triangular K folds to cell Gamma and is allowed on every ABC-cell torus.'};(ROOT/'geometry_validation.json').write_text(json.dumps(proof,indent=2));print(proof,flush=True);return chosen
if __name__=='__main__':catalog()
