import os,time,subprocess,sys,json,pathlib
ROOT=pathlib.Path(__file__).resolve().parent;os.chdir(ROOT)
while True:
 try:r=json.loads((ROOT/'queue_records.json').read_text())
 except (FileNotFoundError,json.JSONDecodeError):r=[]
 if r and r[-1]['task']=='Kprime_source_quadrature':break
 time.sleep(2)
tasks=[
 ('j07_G_seeds',['experiment.py','seeds','--name','j07_G','--j','.7','--Q','0','0']),
 ('j07_G_a3p1',['reduced_experiment.py','j07_G','--angular','3','--reflection','1','--count','6']),
 ('j07_G_a0p1',['reduced_experiment.py','j07_G','--angular','0','--reflection','1','--count','6']),
 ('j07_G_a0m1',['reduced_experiment.py','j07_G','--angular','0','--reflection','-1','--count','6']),
 ('j07_G_a3m1',['reduced_experiment.py','j07_G','--angular','3','--reflection','-1','--count','6']),
 ('j07_G_quadrature',['quadrature.py','j07_G','--dimension','90'])]
records=[]
for label,args in tasks:
 start=time.perf_counter();print('START',label,flush=True)
 with (ROOT/(label+'.log')).open('w') as log:proc=subprocess.run([sys.executable,*args],stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'2'})
 records.append({'task':label,'returncode':proc.returncode,'wall_seconds':time.perf_counter()-start});(ROOT/'extension_records.json').write_text(json.dumps(records,indent=2));print('FINISH',label,records[-1],flush=True)
print('EXTENSION COMPLETE',flush=True)
