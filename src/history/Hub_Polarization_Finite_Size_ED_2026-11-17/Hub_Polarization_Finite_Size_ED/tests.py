import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import numpy as np,scipy.linalg as la,json
from geometry import ROOT,make,key
from ed18 import Basis

def nn_check(g):
 P=g['T']@np.array([[1,1],[2,-1]]);lookup={key(P,c):i for i,c in enumerate(g['coords'])};E=set()
 for i,c in enumerate(g['coords']):
  for d in [(1,0),(0,1),(1,-1),(-1,0),(0,-1),(-1,1)]:E.add(tuple(sorted((i,lookup[key(P,c+np.array(d))]))))
 assert E=={tuple(sorted((a,b))) for a,b,t in g['bonds']};return len(E)
def run():
 g=make([[2,0],[0,2]]);b=Basis(g);states=b.states;idx=b.lookup;N=len(states);rows=[]
 for r in [.37,1.]:
  H=np.zeros((N,N),complex)
  for col,m0 in enumerate(states):
   m=int(m0)
   for a,c,t in g['bonds']:
    j=1 if t==0 else r;diff=((m>>a)^(m>>c))&1;H[col,col]+=j*(.25-.5*diff)
    if diff:H[idx[m^((1<<a)|(1<<c))],col]+=.5*j
  allE=[]
  for q in g['momenta']:
   b.sector(q);b.set_r(r);P=np.column_stack([b.expand(np.eye(b.nb,dtype=complex)[:,i]) for i in range(b.nb)]);red=np.column_stack([b.mv(np.eye(b.nb,dtype=complex)[:,i]) for i in range(b.nb)]);err=float(la.norm(red-P.conj().T@H@P));assert err<1e-10;E=la.eigvalsh(red);allE.extend(E);x=np.random.default_rng(23).normal(size=b.nb)+1j*np.random.default_rng(24).normal(size=b.nb);x/=la.norm(x);comm=float(la.norm(b.spin_action(b.mv(x))-b.mv(b.spin_action(x))));assert comm<1e-10;rows.append({'r':r,'q':q.tolist(),'dimension':b.nb,'full_vs_reduced_H_error':err,'SU2_commutator_error':comm,'isometry_error':float(la.norm(P.conj().T@P-np.eye(b.nb)))})
  err=float(np.max(abs(np.sort(allE)-la.eigvalsh(H))));assert err<1e-10
 shapes=[make([[3,0],[2,2]]),make([[1,2],[0,6]]),make([[3,0],[0,3]]),make([[2,2],[-2,4]])];tri=[{'N':g['N'],'triangular_NN_edges':nn_check(g)} for g in shapes];out={'status':'PASS','independent_12spin_full_vs_all_translation_sectors':rows,'r1_independent_primitive_NN_graph_checks':tri,'full_spectrum_max_error_last_case':err};(ROOT/'tests.json').write_text(json.dumps(out,indent=2));print(out,flush=True)
if __name__=='__main__':run()
