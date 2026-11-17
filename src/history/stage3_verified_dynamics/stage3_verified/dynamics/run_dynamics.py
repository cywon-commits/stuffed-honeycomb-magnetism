"""Reproduce unbroadened Stage 3 results; --compute invokes expensive ED."""
import pathlib,subprocess,sys,os,argparse
ROOT=pathlib.Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--compute',action='store_true');p.add_argument('--tests',action='store_true');a=p.parse_args()
def run(name,*args):subprocess.run([sys.executable,str(ROOT/name),*args],env={**os.environ,'OPENBLAS_NUM_THREADS':'1'},check=True)
if a.tests:run('test_dynamics.py')
if a.compute:
 run('crossing_transitions.py');run('check_characters.py');run('weak_spectrum.py','--moments-only');run('validate_reference.py');run('symmetry_poles.py');__import__('shutil').copy2(ROOT/'symmetry_weak_spectrum.json',ROOT/'weak_spectrum.json');run('verify_weak_symmetry.py');run('block_quadrature.py')
run('validate_reference.py');run('final_validate.py');run('build_dynamics_report.py')
