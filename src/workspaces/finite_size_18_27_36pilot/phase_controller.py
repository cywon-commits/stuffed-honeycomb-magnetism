import os,time,subprocess,sys,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parent;os.chdir(ROOT)
while not all((ROOT/(s+'_crossings.json')).exists() for s in ['K18','C18']):time.sleep(2)
records=[]
for script in ['finish18.py','audit27.py','pilot36.py']:
 start=time.perf_counter();print('START',script,flush=True)
 with (ROOT/(script.replace('.py','.log'))).open('w') as f:r=subprocess.run([sys.executable,script],stdout=f,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'2'})
 records.append({'phase':script,'returncode':r.returncode,'wall_seconds':time.perf_counter()-start});(ROOT/'phase_records.json').write_text(json.dumps(records,indent=2));print('FINISH',records[-1],flush=True)
 if r.returncode:break
