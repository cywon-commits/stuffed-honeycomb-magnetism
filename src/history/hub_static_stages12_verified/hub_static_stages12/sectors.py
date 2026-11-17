"""Verify weighted point symmetries of the stuffed honeycomb and cell momentum stars."""
import numpy as np,json
from ed import geometry

def cell_stars(L=3,jp=.2):
 c,b=geometry(L,jp);h=np.array([0,1])
 def key(v):u,w=map(int,v);return ((u+2*w)%(3*L),(u-w)%(3*L))
 lookup={key(v):i for i,v in enumerate(c)}
 edges={(*sorted((a,d)),float(j)) for a,d,j in b}
 perms=[]
 for t in [lambda u,v:(-v,u+v),lambda u,v:(v,u),lambda u,v:(-u,-v)]:
  perm=np.array([lookup[key(np.array(t(*(v-h)))+h)] for v in c]);assert {(*sorted((perm[a],perm[d])),j) for a,d,j in b}==edges;perms.append(perm)
 seen=set();stars=[]
 for a in range(L):
  for d in range(L):
   if (a,d) in seen:continue
   orbit={(a,d)};todo=[(a,d)]
   while todo:
    k=np.array(todo.pop())/L;q=np.array([(k[0]+k[1])/3,(2*k[0]-k[1])/3])
    for p in [(q[0]-q[1],q[0]),(q[1],q[0])]:
     kk=(int(round(L*(p[0]+p[1])))%L,int(round(L*(2*p[0]-p[1])))%L)
     if kk not in orbit:orbit.add(kk);todo.append(kk)
   seen|=orbit;stars.append(sorted(orbit))
 return stars,perms[-1]
if __name__=='__main__':print(json.dumps(cell_stars()[0]))

def point_permutations(L=3,jp=.2):
 c,b=geometry(L,jp);h=np.array([0,1])
 def key(v):u,w=map(int,v);return ((u+2*w)%(3*L),(u-w)%(3*L))
 lookup={key(v):i for i,v in enumerate(c)};edges={(*sorted((a,d)),float(j)) for a,d,j in b};out=[]
 for reflection in [False,True]:
  for power in range(6):
   perm=[]
   for v in c:
    u,w=map(int,v-h)
    if reflection:u,w=w,u
    for _ in range(power):u,w=-w,u+w
    perm.append(lookup[key(np.array([u,w])+h)])
   perm=np.array(perm);assert {(*sorted((perm[a],perm[d])),j) for a,d,j in b}==edges;out.append(perm)
 return out
