"""Exact integer momentum stars; geometric point-group verification."""
import numpy as np,json
from ed import geometry

def momentum_stars():
 coords,bonds=geometry(3)
 def key(c):u,v=map(int,c);return ((u+2*v)%9,(u-v)%9)
 lookup={key(c):i for i,c in enumerate(coords)}
 edges={tuple(sorted((a,b))) for a,b,j in bonds}
 for transform in [lambda u,v:(-v,u+v),lambda u,v:(v,u)]:
  perm=[lookup[key(transform(*c))] for c in coords]
  assert {tuple(sorted((perm[a],perm[b]))) for a,b in edges}==edges
 qs={((a+b+3*f)%9,(2*a-b+6*f)%9) for a in range(3) for b in range(3) for f in range(3)}
 seen=set();groups=[]
 for q in sorted(qs):
  if q in seen:continue
  orbit={q};todo=[q]
  while todo:
   u,v=todo.pop()
   for p in [((u-v)%9,u%9),(v,u)]:
    assert p in qs
    if p not in orbit:orbit.add(p);todo.append(p)
  seen|=orbit;groups.append(sorted(orbit))
 assert len(seen)==27
 return groups
if __name__=='__main__':print(json.dumps(momentum_stars(),indent=2))
