"""Reproduce cached tables/checkpoints; --compute resumes expensive calculation phases."""
import os,sys,pathlib,argparse,subprocess
ROOT=pathlib.Path(__file__).resolve().parent;os.chdir(ROOT)
p=argparse.ArgumentParser();p.add_argument('--compute',action='store_true');p.add_argument('--tests',action='store_true');a=p.parse_args()
def run(script,*args):subprocess.run([sys.executable,str(ROOT/script),*args],check=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'2','MPLCONFIGDIR':'/tmp/hub_size_mpl','XDG_CACHE_HOME':'/tmp/hub_size_cache'})
run('geometry.py');run('ingest27.py');run('validate_reference_benchmarks.py')
if a.tests:run('tests.py')
if a.compute:
 for shape in ['K18','C18']:run('scan18.py','--shape',shape)
 run('finish18.py');run('repair_near_crossing_spins.py');run('audit27.py');run('pilot36.py')
run('tighten_crossings.py')
if a.tests:run('test36pilot.py')
run('build_results.py');run('validate_cached.py');run('plot_results.py');run('report.py');print('PASS reproducible outputs and fresh checkpoint verification.')
