import pathlib,json,csv,numpy as np
ROOT=pathlib.Path(__file__).resolve().parent
rows=list(csv.DictReader((ROOT/'reference/stage3_targets_moments.csv').open()))
def unpack(x):return np.array(x['real'])+1j*np.array(x['imag'])
out=[];w=json.loads((ROOT/'weak_moments.json').read_text());m0=unpack(w['moments']['1']['M0']);mf=unpack(w['moments']['1']['M1']);mr=unpack(w['moments']['-1']['M1']);sym=(mf+mr.T)/2
for a,an in enumerate(['n','m','h']):
 for b,bn in enumerate(['n','m','h']):
  r=next(r for r in rows if float(r['Jprime_over_J'])==.2375 and r['probe_k']=='(0.0, 0.3333333333333333)' and r['channel']==an+bn)
  for label,value,key in [('M0zz',m0[a,b],'M0zz'),('M1zz_sym',sym[a,b],'M1zz_sym')]:
   ref=complex(float(r[key+'_re']),float(r[key+'_im']));out.append({'Jp':.2375,'channel':an+bn,'moment':label,'computed_real':value.real,'computed_imag':value.imag,'reference_real':ref.real,'reference_imag':ref.imag,'absolute_error':abs(value-ref)})
for j in [.3875,.3921875]:
 d=json.loads((ROOT/f'crossing_{j:.7f}.json').read_text())
 for c in ['n','m','h']:
  r=next(r for r in rows if float(r['Jprime_over_J'])==j and r['probe_k']=='(0.0, 0.0)' and r['channel']==c+c)
  for moment,key in [('M0','M0zz'),('M1','M1zz_sym')]:
   v=d['moments']['H' if c=='h' else c][moment];ref=float(r[key+'_re']);out.append({'Jp':j,'channel':c+c,'moment':key,'computed_real':v,'computed_imag':0.,'reference_real':ref,'reference_imag':0.,'absolute_error':abs(v-ref)})
with (ROOT/'moment_reference_comparison.csv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=out[0]);writer.writeheader();writer.writerows(out)
summary={'maximum_error':max(d['absolute_error'] for d in out),'entries':len(out),'second_moment_reference':'Not supplied; M2 is newly evaluated as (D B)^dagger(D B), validated by small-cluster complete spectral sums.','forward_reverse_M1_matrix_difference':float(np.max(abs(mf-mr.T)))}
assert summary['maximum_error']<1e-8
(ROOT/'moment_reference_validation.json').write_text(json.dumps(summary,indent=2));print(summary)
