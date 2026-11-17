import pathlib,hashlib,json,zipfile,sys
ROOT=pathlib.Path(__file__).resolve().parent;PREFIX='stage3_followup/'
exclude={'spinhalf_refine.py','finish_run.py','archive_verification.json','final_pipeline.log','export.log'}
files=[]
for p in ROOT.iterdir():
 if p.is_file() and p.suffix in {'.py','.cpp','.json','.csv','.md','.txt','.png','.log'} and p.name not in exclude:files.append((p,'followup/'+p.name))
for name in ['weak_spectrum.json','weak_moments.json','weak_quadrature.json','CONVENTIONS.md','small_dynamics_validation.json']:
 files.append((ROOT.parent/'dynamics'/name,'dynamics/'+name))
files.extend([(ROOT.parent/'static/static_j0.2_realspace.npz','static/static_j0.2_realspace.npz'),(ROOT.parent/'crossing/obs_p0.237500000_n12_r0_e0_realspace.npz','crossing/obs_p0.237500000_n12_r0_e0_realspace.npz')])
checkpoint=[]
for folder,names in [('static',['j0.2_n12_r0','j0.2_n01_r0','j0.2_n01_r1']),('crossing',['p0.237500000_n12_r0_e0','p0.237500000_n01_r0_e0','p0.237500000_n01_r1_e0'])]:
 for name in names:
  for ext in ['.json','.npz']:checkpoint.append((ROOT.parent/folder/(name+ext),folder+'/'+name+ext))
checkpoint.extend((p,'followup/'+p.name) for p in ROOT.glob('j0.*_r*_root*.npz'))
results=[]
plans=[('stage3_followup_results.zip',files)]
if '--compact-only' not in sys.argv:plans.append(('stage3_followup_full_checkpoints.zip',files+checkpoint))
for filename,items in plans:
 manifest='\n'.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+n for p,n in sorted(items,key=lambda v:v[1]))+'\n';dest=ROOT/filename
 with zipfile.ZipFile(dest,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p,n in items:z.write(p,PREFIX+n)
  z.writestr(PREFIX+'SHA256SUMS',manifest)
 with zipfile.ZipFile(dest) as z:
  assert z.testzip() is None
  for line in manifest.splitlines():
   h,n=line.split('  ',1);assert hashlib.sha256(z.read(PREFIX+n)).hexdigest()==h
 results.append({'file':filename,'bytes':dest.stat().st_size,'sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),'entries':len(items)+1,'CRC_and_SHA_verification':'PASS'})
(ROOT/'archive_verification.json').write_text(json.dumps(results,indent=2));print(json.dumps(results,indent=2))
