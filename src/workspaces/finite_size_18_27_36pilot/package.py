import pathlib,hashlib,json,zipfile
ROOT=pathlib.Path(__file__).resolve().parent;OUT=ROOT/'Hub_Polarization_Finite_Size_ED_2026-11-17.zip';PREFIX='Hub_Polarization_Finite_Size_ED';files=[]
for p in sorted(ROOT.rglob('*')):
 if not p.is_file():continue
 r=p.relative_to(ROOT)
 if any(x in ['__pycache__','.git'] for x in r.parts):continue
 if p.suffix=='.zip' or p.name in ['SHA256SUMS','archive_verification.json']:continue
 if p.suffix not in ['.py','.cpp','.so','.md','.txt','.csv','.json','.npz','.png','.svg','.pdf','.log']:continue
 files.append(p)
manifest=[hashlib.sha256(p.read_bytes()).hexdigest()+'  '+p.relative_to(ROOT).as_posix() for p in files];(ROOT/'SHA256SUMS').write_text('\n'.join(manifest)+'\n');files.append(ROOT/'SHA256SUMS')
with zipfile.ZipFile(OUT,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in files:z.write(p,PREFIX+'/'+p.relative_to(ROOT).as_posix())
with zipfile.ZipFile(OUT) as z:
 assert z.testzip() is None
 for line in manifest:
  h,n=line.split('  ',1);assert hashlib.sha256(z.read(PREFIX+'/'+n)).hexdigest()==h
result={'status':'PASS','archive':OUT.name,'size_bytes':OUT.stat().st_size,'sha256':hashlib.sha256(OUT.read_bytes()).hexdigest(),'file_count':len(files),'CRC_and_all_member_SHA256_verified':True};(ROOT/'archive_verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
