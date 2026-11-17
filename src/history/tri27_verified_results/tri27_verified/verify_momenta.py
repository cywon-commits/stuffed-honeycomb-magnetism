"""Determine triangular nearest-neighbor vectors and momentum minima independently."""
import pathlib,numpy as np,json
from ed import geometry
from symmetry import momentum_stars
root=pathlib.Path(__file__).resolve().parent
coords,bonds=geometry(3);neighbors=[b if a==0 else a for a,b,j in bonds if a==0 or b==0]
deltas=set()
for j in neighbors:
 for x in range(-2,3):
  for y in range(-2,3):
   u,v=coords[j]-coords[0]+np.array([3*x+6*y,3*x-3*y])
   if u*u+v*v+u*v==1:deltas.add((int(u),int(v)))
assert deltas=={(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)}
deltas=np.array(sorted(deltas))
def sf(q):return float(np.exp(2j*np.pi*(deltas@q)).sum().real)
qs=[np.array(q)/9 for star in momentum_stars() for q in star]
vals=[sf(q) for q in qs];low=min(vals)
minima=[q.tolist() for q,v in zip(qs,vals) if abs(v-low)<1e-12]
assert len(minima)==2 and abs(low+3)<1e-12
for q in qs:assert abs(sf(q)-2*(np.cos(2*np.pi*q[0])+np.cos(2*np.pi*q[1])+np.cos(2*np.pi*(q[0]-q[1]))))<1e-12
assert abs(sf(np.array([1/3,1/3])))<1e-12
result={'primitive_vectors_cartesian':[[1.,0.],[.5,float(np.sqrt(3)/2)]],'cell_vectors_in_primitive_coordinates':[[1,1],[2,-1]],'period_vectors_in_primitive_coordinates':[[3,3],[6,-3]],'nearest_neighbor_displacements':deltas.tolist(),'structure_factor_formula':'2*[cos(2*pi*q1)+cos(2*pi*q2)+cos(2*pi*(q1-q2))]','primitive_Gamma_structure_factor':sf(np.array([0,0])),'primitive_K_structure_factor':low,'primitive_K_minima':minima,'primitive_(1/3,1/3)_structure_factor':sf(np.array([1/3,1/3])),'all_27_allowed_momenta_checked':True,'cell_fold_of_both_primitive_K_points':[0,0],'M_allowed':False}
(root/'momentum_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
