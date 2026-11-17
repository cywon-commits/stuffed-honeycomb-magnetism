"""Fresh consistent benchmarks for primitive triangular Gamma and K."""
import pathlib,subprocess,sys,os
root=pathlib.Path(__file__).resolve().parent
for name,q in [('Gamma_primitive',(0.,0.)),('K_primitive',(1/3,2/3))]:
 with (root/(name+'.log')).open('w') as log:
  subprocess.run([sys.executable,str(root/'ed.py'),'--q',*map(str,q),'--name',name,'--maxiter','180'],env={**os.environ,'OPENBLAS_NUM_THREADS':'1'},stdout=log,stderr=subprocess.STDOUT,check=True)
