"""Reproduce all six point-group momentum classes; skip converged result files."""
import subprocess,sys,pathlib,os,json
from symmetry import momentum_stars
root=pathlib.Path(__file__).resolve().parent
momentum_stars()
plan=[('Gamma_primitive',(0,0)),('K_primitive',(3,6)),('cell11_fold1',(2,1)),('cell11_fold2',(5,7)),('cell11_fold3',(8,4)),('other',(0,3))]
for name,q in plan:
 path=root/(name+'.json')
 if path.exists():
  r=json.loads(path.read_text())
  if r['converged'] and r['residual']<1e-8:continue
 with (root/(name+'.log')).open('w') as log:
  subprocess.run([sys.executable,str(root/'ed.py'),'--q',*[str(x/9) for x in q],'--name',name,'--maxiter','180'],env={**os.environ,'OPENBLAS_NUM_THREADS':'1'},stdout=log,stderr=subprocess.STDOUT,check=True)
