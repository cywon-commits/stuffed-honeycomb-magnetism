"""Build and checksum the portable verified-result archive."""
import pathlib,json,hashlib,zipfile,shutil
ROOT=pathlib.Path(__file__).resolve().parent
OUT=ROOT/'delivery'
if OUT.exists():shutil.rmtree(OUT)
OUT.mkdir()
selected={'p0.237500000_n12_r0_e0.npz','p0.275000000_n00_r0_e0.npz','p0.350000000_n00_r0_e0.npz','p0.425000000_n00_r3_e0.npz'}
exclude={'finish_targeted.py','finish_targeted.log','preview_brackets.py','crossing_intervals_preview.json','preview.log','export.log','plots.log','final_pipeline.log','archive_verification.json'}
files=[]
for p in ROOT.iterdir():
 if not p.is_file() or p.name in exclude:continue
 if p.suffix in {'.py','.cpp','.md','.txt','.json','.csv','.svg','.png','.log'} or p.name in selected or (p.suffix=='.npz' and p.name.startswith('obs_')):files.append(p)
files+=list((ROOT/'endpoint_reference').glob('*.json'))
for p in files:
 target=OUT/p.relative_to(ROOT);target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,target)
cpinfo={'included_checkpoints':sorted(selected),'omitted_intermediate_checkpoints':'Regenerate with run_analysis.py --compute. Verified scalar records and direct residuals for every completed level are included.','previous_endpoint_ground_checkpoints':'Available in the prior static bundle linked in README.md.','no_dynamics':True}
(OUT/'checkpoint_manifest.json').write_text(json.dumps(cpinfo,indent=2))
manifest=[]
for p in sorted(OUT.rglob('*')):
 if p.is_file() and p.name!='SHA256SUMS':manifest.append(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(OUT).as_posix())
(OUT/'SHA256SUMS').write_text('\n'.join(manifest)+'\n')
zip_path=ROOT/'hub_ground_state_crossings_verified.zip'
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(OUT.rglob('*')):
  if p.is_file():z.write(p,'hub_ground_state_crossings/'+str(p.relative_to(OUT)))
with zipfile.ZipFile(zip_path) as z:assert z.testzip() is None
summary={'archive':str(zip_path),'bytes':zip_path.stat().st_size,'sha256':hashlib.sha256(zip_path.read_bytes()).hexdigest(),'entries':len(manifest)+1,'CRC_check':'PASS'}
(ROOT/'archive_verification.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
