import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json,hashlib,numpy as np,scipy.linalg as la
from experiment import ROOT,Model,geometry_weights,dump
from little_model import LittleModel
from observables import convert_basis,total_spin_check
from group_poles import unpack
records=[];completeness=[]
previous={}
if (ROOT/'little_measurement_validation.json').exists():previous={r['tag']:r for r in json.loads((ROOT/'little_measurement_validation.json').read_text())['fresh_full_translation_space_checks']}
for name,a,p in [('j1_G',0,-1),('j1_Q',0,1),('j1_Qprime',0,1),('j09_G',3,1),('j09_q',0,1),('j07_G',3,1)]:
 tag=name+f'_a{a}_p{p}';f=ROOT/(tag+'_poles.npz');mf=ROOT/(name+'_moments.json')
 if not f.exists() or not mf.exists():continue
 fingerprint=hashlib.sha256(f.read_bytes()+(ROOT/'kernel.cpp').read_bytes()+(ROOT/'little_model.py').read_bytes()).hexdigest()
 if tag in previous and previous[tag].get('fingerprint')==fingerprint:records.append(previous[tag]);continue
 d=json.loads(mf.read_text());j=d['j'];q=np.array(d['target_q']);m=LittleModel(q,a,p,primitive=j==1,jp=j);full=Model(primitive=j==1,jp=j);full.sector(q);z=np.load(f);v=z['vectors'][:,0];v/=la.norm(v);psi=convert_basis(m,full,v);Hsmall=convert_basis(m,full,m.mv(v));Hfull=full.mv(psi);err=float(la.norm(Hsmall-Hfull));res=float(la.norm(Hfull-float(z['energy'][0])*psi));spin=total_spin_check(full,psi);assert err<1e-9 and abs(la.norm(psi)-1)<1e-9 and res<1e-8;records.append({'tag':tag,'fingerprint':fingerprint,'full_translation_dimension':full.nb,'symmetry_dimension':m.nb,'H_intertwining_error':err,'full_basis_H_residual':res,'full_basis_spin':spin});del m,full,psi,Hfull,Hsmall,z
for name in ['j1_Q','j1_Qprime','j09_G','j09_q','j07_G']:
 mf=ROOT/(name+'_moments.json')
 if not mf.exists():continue
 d=json.loads(mf.read_text());M=unpack(d['M0']);sumM=np.zeros_like(M);seen=set()
 for f in ROOT.glob(name+'_a*_spectrum.json'):
  dd=json.loads(f.read_text());sumM+=unpack(dd['projected_M0']);seen.add((dd['angular'],dd['reflection']))
 for f in ROOT.glob(name+'_a*_empty.json'):
  dd=json.loads(f.read_text());sumM+=unpack(dd['M0'])
 # Only physical Gamma chirality pair is required/covered for optional j=.7.
 if name=='j07_G':P=geometry_weights(np.zeros(2))[:,4:6];err=la.norm(P.conj().T@(sumM-M)@P);scope='physical-Gamma chirality face pair only'
 else:err=la.norm(sumM-M);scope='all scalar seed channels'
 expected=8 if name=='j09_G' else (4 if name=='j07_G' else 2)
 nfiles=len(list(ROOT.glob(name+'_a*_spectrum.json')))+len(list(ROOT.glob(name+'_a*_empty.json')))
 if nfiles==expected:assert err<1e-8
 completeness.append({'name':name,'represented_blocks':nfiles,'expected_blocks':expected,'scope':scope,'M0_completeness_error':float(err),'complete':nfiles==expected})
result={'status':'PASS','fresh_full_translation_space_checks':records,'projector_completeness':completeness};dump(ROOT/'little_measurement_validation.json',result);print(json.dumps(result,indent=2))
