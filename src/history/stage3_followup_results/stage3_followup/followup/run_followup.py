"""Cached reporting by default; explicit expensive recomputation with --compute."""
import pathlib,subprocess,sys,os,argparse,tempfile,shutil
ROOT=pathlib.Path(__file__).resolve().parent
p=argparse.ArgumentParser();p.add_argument('--compute',action='store_true');p.add_argument('--fresh',action='store_true');p.add_argument('--tests',action='store_true');a=p.parse_args()
def run(script,*args):subprocess.run([sys.executable,str(ROOT/script),*args],env={**os.environ,'OPENBLAS_NUM_THREADS':'1','MPLCONFIGDIR':str(ROOT/'.mpl'),'XDG_CACHE_HOME':str(ROOT/'.cache')},check=True)
if a.fresh:
 if not a.compute:p.error('--fresh requires --compute')
 backup=pathlib.Path(tempfile.mkdtemp(prefix='prior_run_',dir=ROOT))
 for file in ROOT.glob('j0.*_r*_root*'):
  if file.suffix in ['.json','.npz']:shutil.move(str(file),backup/file.name)
if a.tests:run('test_followup.py')
if a.compute:
 run('matched.py');run('validate_moments.py');run('search.py','0','--count','2');run('search.py','1','--count','2');run('search_j02.py','1','--count','1');run('search_j02.py','0','--count','1');run('search_j02_half.py','0','--count','2');run('search_j2375_threehalf.py','0','--count','3');run('search_j02.py','0','--count','3');run('partner_check.py');run('post_validate.py')
run('validate_moments.py');run('report.py')
