import pathlib,json,zipfile,hashlib
ROOT=pathlib.Path(__file__).resolve().parent;DEST=ROOT/'stage3_verified_dynamics.zip';files=[]
for p in ROOT.iterdir():
 if p.is_file() and p.suffix in {'.py','.cpp','.json','.csv','.md','.txt','.png','.log'} and p.name not in {'archive_verification.json','report_build.log'}:files.append((p,'dynamics/'+p.name))
files += [(p,'dynamics/reference/'+p.name) for p in (ROOT/'reference').iterdir() if p.is_file()]
names=['p0.237500000_n12_r0_e0','p0.237500000_n01_r0_e0','p0.237500000_n01_r1_e0']+[f'p{j:.9f}_n00_r{r}_e0' for j in [.3875,.3921875] for r in [0,3]]
for n in names:
 for ext in ['.json','.npz']:files.append((ROOT.parent/'crossing'/(n+ext),'crossing/'+n+ext))
manifest='\n'.join(hashlib.sha256(p.read_bytes()).hexdigest()+'  '+name for p,name in sorted(files,key=lambda x:x[1]))+'\n'
with zipfile.ZipFile(DEST,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p,name in files:z.write(p,'stage3_verified/'+name)
 z.writestr('stage3_verified/SHA256SUMS',manifest)
with zipfile.ZipFile(DEST) as z:
 assert z.testzip() is None
 for line in manifest.splitlines():
  h,n=line.split('  ',1);assert hashlib.sha256(z.read('stage3_verified/'+n)).hexdigest()==h
result={'bytes':DEST.stat().st_size,'sha256':hashlib.sha256(DEST.read_bytes()).hexdigest(),'file_count':len(files)+1,'archive_crc_and_hash_checks':'PASS'}
(ROOT/'archive_verification.json').write_text(json.dumps(result,indent=2));print(result)
