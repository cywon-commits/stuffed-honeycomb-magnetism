"""Reproduce static Stages 1–2; no spectral/dynamical computation is invoked."""
import pathlib,subprocess,sys,os,json,argparse
from coverage import PLAN,verified_ground
ROOT=pathlib.Path(__file__).resolve().parent
ENV={**os.environ,'OPENBLAS_NUM_THREADS':'1','MPLCONFIGDIR':str(ROOT/'.mpl'),'XDG_CACHE_HOME':str(ROOT/'.cache')}
def run(args,log=None):
 if log:
  with (ROOT/log).open('w') as f:subprocess.run([sys.executable,*map(str,args)],env=ENV,stdout=f,stderr=subprocess.STDOUT,check=True)
 else:subprocess.run([sys.executable,*map(str,args)],env=ENV,check=True)
def main():
 p=argparse.ArgumentParser();p.add_argument('--force',action='store_true');p.add_argument('--skip-tests',action='store_true');a=p.parse_args()
 if not a.skip_tests:
  for test in ['test_observables.py','test_cyclic.py','test_solver.py']:run([ROOT/test],test+'.log')
 if a.force or not (ROOT/'static_j1.json').exists():run([ROOT/'analyze_state.py','--jp','1','--file',ROOT/'isotropic_primitive_ground.npz','--primitive'],'static_j1.log')
 for jp in [.2,.5]:
  for n,ell in PLAN:
   name=f'j{jp:g}_n{n[0]}{n[1]}_r{ell}';f=ROOT/(name+'.json')
   valid=f.exists() and json.loads(f.read_text()).get('converged',False)
   if a.force or not valid:run([ROOT/'solve_sector.py','--jp',jp,'--n',*n,'--cyclic','--angular',ell,'--ncv',('32' if jp==.2 and tuple(n) in [(0,0),(0,1)] else '64'),'--seconds','600'],name+'.log')
  cov=verified_ground(jp)
  if cov is None:print(f'Jp={jp}: incomplete momentum coverage; no ground state declared.');continue
  ground=ROOT/(cov['ground_name']+'.npz')
  if not ground.exists():
   record=json.loads((ROOT/(cov['ground_name']+'.json')).read_text());run([ROOT/'solve_sector.py','--jp',jp,'--n',*record['cell_n'],'--cyclic','--angular',record['angular'],'--ncv','64'],cov['ground_name']+'.log');cov=verified_ground(jp)
   if cov is None:continue
  if a.force or not (ROOT/f'static_j{jp:g}.json').exists():run([ROOT/'analyze_state.py','--jp',jp,'--file',ground],f'static_j{jp:g}.log')
 if verified_ground(.5) is not None and (ROOT/'static_j0.5.json').exists():run([ROOT/'verify_spin_multiplet.py'],'spin_multiplet_validation.log')
 run([ROOT/'test_driver.py'],'driver_validation.log')
 run([ROOT/'enrich_results.py']);run([ROOT/'build_report.py']);run([ROOT/'plot_static.py'])
if __name__=='__main__':main()
