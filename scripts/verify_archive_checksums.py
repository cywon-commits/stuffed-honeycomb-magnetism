#!/usr/bin/env python3
import csv, hashlib, pathlib, sys, zipfile
root=pathlib.Path(__file__).resolve().parents[1]
if '--workspace-only' in sys.argv:
 p=root/'workspaces'/'finite_size_18_27_36pilot'
 print('Latest finite-size workspace exists:',p.is_dir(),'ed18.py:',(p/'ed18.py').is_file(),'36 pilot:',(p/'pilot36.py').is_file())
 sys.exit(0 if p.is_dir() and (p/'ed18.py').is_file() else 1)
rows=list(csv.DictReader((root/'manifests'/'original_zip_manifest.csv').open()))
errors=0
for row in rows:
 path=root/'original_archives'/row['filename']
 if not path.exists():print('MISSING',path);errors+=1;continue
 sha=hashlib.sha256()
 with path.open('rb') as f:
  for block in iter(lambda:f.read(1024*1024),b''):sha.update(block)
 ok=sha.hexdigest()==row['sha256']
 print('OK' if ok else 'BAD',row['filename'])
 errors+=not ok
print('Total files:',len(rows),'errors:',errors)
sys.exit(bool(errors))
