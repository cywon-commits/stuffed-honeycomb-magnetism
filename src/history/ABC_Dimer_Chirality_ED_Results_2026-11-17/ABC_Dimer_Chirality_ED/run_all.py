"""Default: rebuild cached CSV/JSON/figures and validate. --compute reruns costly ED."""
import pathlib,sys,subprocess,os,argparse
ROOT=pathlib.Path(__file__).resolve().parent
os.chdir(ROOT)
p=argparse.ArgumentParser();p.add_argument('--tests',action='store_true');p.add_argument('--compute',action='store_true');a=p.parse_args()
def run(script,*args):subprocess.run([sys.executable,str(ROOT/script),*args],check=True,env={**os.environ,'OPENBLAS_NUM_THREADS':'1','OMP_NUM_THREADS':'2','MPLCONFIGDIR':'/tmp/abc_dynamics_matplotlib','XDG_CACHE_HOME':'/tmp/abc_dynamics_cache'})
if a.tests:
 for f in ['test_triangle.py','test_fast.py','test_primitive.py']:run(f)
 run('reference/ABC_Dimer_Chirality_Dynamics_2026-11-16/code/triangle_exact_dynamics.py')
 run('reference/ABC_Dimer_Chirality_Dynamics_2026-11-16/code/static_mode_summary.py')
if a.compute:
 run('experiment.py','source');run('experiment.py','seeds','--name','j1_G','--Q','0','0');run('experiment.py','solve','--name','j1_G','--count','12','--ncv','72');run('reduced_experiment.py','j1_G','--angular','0','--reflection','-1','--count','20');run('experiment.py','seeds','--name','j1_Q','--Q',str(8/9),str(7/9));run('test_little.py');run('compute_queue.py');run('extension_queue.py');run('local_sum_rules.py');run('reuse_K_poles.py');run('validate_little_measurements.py')
run('group_poles.py');run('validate_inputs.py');run('convergence_stability.py');run('projector_ground_refs.py');run('build_data.py');run('face_coverage.py');run('verify_results.py');run('symmetry_partners.py');run('export_partners.py');run('plot_results.py')
if (ROOT/'report.py').exists():run('report.py')
print('Reproduction and validation complete.')
