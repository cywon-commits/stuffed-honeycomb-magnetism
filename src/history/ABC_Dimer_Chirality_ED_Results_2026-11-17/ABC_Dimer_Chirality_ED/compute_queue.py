"""Sequential priority queue with durable logs; failures are explicit."""
import os,time,subprocess,sys,pathlib,json
ROOT=pathlib.Path(__file__).resolve().parent
os.chdir(ROOT)
# Allow the already-started priority-1 chirality solve to finish.
while not (ROOT/'j1_G_a0_p-1_spectrum.json').exists():time.sleep(2)
tasks=[
 ('G_chi_orthogonal',['reduced_experiment.py','j1_G','--angular','0','--reflection','1','--count','8']),
 ('Q_mirror_plus',['reduced_experiment.py','j1_Q','--reflection','1','--count','8']),
 ('Q_mirror_minus',['reduced_experiment.py','j1_Q','--reflection','-1','--count','8']),
 ('Qprime_seeds',['experiment.py','seeds','--name','j1_Qprime','--Q',str(8/9),str(7/9),'--partner']),
 ('Qprime_mirror_plus',['reduced_experiment.py','j1_Qprime','--reflection','1','--count','8']),
 ('Qprime_mirror_minus',['reduced_experiment.py','j1_Qprime','--reflection','-1','--count','8']),
 ('G_quadrature',['quadrature.py','j1_G','--dimension','120']),
 ('Q_quadrature',['quadrature.py','j1_Q','--dimension','120']),
 ('Qprime_quadrature',['quadrature.py','j1_Qprime','--dimension','120']),
 ('j09_G_seeds',['experiment.py','seeds','--name','j09_G','--j','.9','--Q','0','0']),
 ('j09_G_a0p1',['reduced_experiment.py','j09_G','--angular','0','--reflection','1','--count','6']),
 ('j09_G_a0m1',['reduced_experiment.py','j09_G','--angular','0','--reflection','-1','--count','6']),
 ('j09_G_a3p1',['reduced_experiment.py','j09_G','--angular','3','--reflection','1','--count','6']),
 ('j09_G_a3m1',['reduced_experiment.py','j09_G','--angular','3','--reflection','-1','--count','6']),
 ('j09_G_a1',['reduced_experiment.py','j09_G','--angular','1','--reflection','none','--count','6']),
 ('j09_G_a2',['reduced_experiment.py','j09_G','--angular','2','--reflection','none','--count','6']),
 ('j09_G_a4',['reduced_experiment.py','j09_G','--angular','4','--reflection','none','--count','6']),
 ('j09_G_a5',['reduced_experiment.py','j09_G','--angular','5','--reflection','none','--count','6']),
 ('j09_G_quadrature',['quadrature.py','j09_G','--dimension','90']),
 ('j09_q_seeds',['experiment.py','seeds','--name','j09_q','--j','.9','--Q',str(2/3),'0']),
 ('j09_q_plus',['reduced_experiment.py','j09_q','--reflection','1','--count','6']),
 ('j09_q_minus',['reduced_experiment.py','j09_q','--reflection','-1','--count','6']),
 ('j09_q_quadrature',['quadrature.py','j09_q','--dimension','90']),
 ('K_seeds',['experiment.py','seeds','--name','j1_K','--Q',str(1/3),str(2/3)]),
 ('Kprime_source_seeds',['experiment.py','seeds','--name','j1_KprimeSource','--Q',str(1/3),str(2/3),'--partner']),
 ('K_quadrature',['quadrature.py','j1_K','--dimension','96']),
 ('Kprime_source_quadrature',['quadrature.py','j1_KprimeSource','--dimension','96']),
]
records=[]
for label,args in tasks:
 start=time.perf_counter();print('START',label,flush=True)
 with (ROOT/(label+'.log')).open('w') as log:
  proc=subprocess.run([sys.executable,*args],stdout=log,stderr=subprocess.STDOUT,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'2'})
 records.append({'task':label,'returncode':proc.returncode,'wall_seconds':time.perf_counter()-start});(ROOT/'queue_records.json').write_text(json.dumps(records,indent=2));print('FINISH',label,records[-1],flush=True)
 if proc.returncode:print('ERROR inspect '+label+'.log',flush=True)
print('QUEUE COMPLETE',flush=True)
