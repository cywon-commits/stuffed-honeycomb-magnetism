"""Reproduce the targeted level-crossing calculation; never starts Stage 3."""
import pathlib,os,subprocess,sys,argparse
ROOT=pathlib.Path(__file__).resolve().parent
ENV={**os.environ,'OPENBLAS_NUM_THREADS':'1','MPLCONFIGDIR':str(ROOT/'.mpl'),'XDG_CACHE_HOME':str(ROOT/'.cache')}
def run(file):subprocess.run([sys.executable,str(ROOT/file)],env=ENV,check=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--compute',action='store_true',help='Run eigensolvers/refinement, reusing verified checkpoints. Omit for tables/plots from cached results.');p.add_argument('--tests',action='store_true');a=p.parse_args()
 if a.tests:
  for script in ['test_crossing.py','test_complex_deflation.py']:run(script)
 if a.compute:
  run('cross_scan.py');run('refine_crossings.py');run('postprocess.py');run('profile_parities.py')
 run('bracket_exclusions.py');run('enrich_observables.py');run('export_results.py');run('build_report.py');run('plot_crossings.py');run('validate_results.py')
if __name__=='__main__':main()
