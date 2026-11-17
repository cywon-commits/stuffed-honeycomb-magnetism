import numpy as np,scipy.linalg as la,json

def full_pair_action(states,lookup,psi,pairs,onsite=0.):
 out=onsite*psi.copy()
 for a,b,coef in pairs:
  ba=((states>>a)&1).astype(float)-.5;bb=((states>>b)&1).astype(float)-.5;out+=coef*ba*bb*psi;diff=ba!=bb;t=states[diff]^np.uint32((1<<a)|(1<<b));out[diff]+=.5*coef*psi[lookup[t]]
 return out

def correlations(b,v):
 psi=b.expand(v);states=b.states;G=np.eye(b.n)*.75;prob=abs(psi)**2;max_im=0.
 for a in range(b.n):
  for c in range(a+1,b.n):
   za=((states>>a)&1).astype(float)-.5;zc=((states>>c)&1).astype(float)-.5;diff=za!=zc;t=states[diff]^np.uint32((1<<a)|(1<<c));val=np.dot(prob,za*zc)+.5*np.vdot(psi[diff],psi[b.lookup[t]]);G[a,c]=G[c,a]=val.real;max_im=max(max_im,abs(val.imag))
 return G,psi,max_im

def values(g,G,S):
 nh=g['Ncell'];nr=2*nh;h=np.arange(2,g['N'],3);rim=np.array([i for i in range(g['N']) if i%3!=2]);sgn=np.array([1 if i%3==0 else(-1 if i%3==1 else 0) for i in range(g['N'])]);H2=float(G[np.ix_(h,h)].sum());R2=float(G[np.ix_(rim,rim)].sum());HR=float(G[np.ix_(h,rim)].sum());neel=float(sgn@G@sgn/nr**2);hs=[]
 for q in g['momenta']:
  phase=np.exp(-2j*np.pi*(g['cells']@q))/np.sqrt(nh);hs.append(float(np.vdot(phase,G[np.ix_(h,h)]@phase).real))
 peak=max(hs);peakq=[q.tolist() for q,s in zip(g['momenta'],hs) if abs(s-peak)<1e-8]
 return {'S':S,'m_tot':S/g['N'],'m_tot_rms':np.sqrt(S*(S+1))/g['N'],'H2':H2,'R2':R2,'HR':HR,'spin_sum_rule_error':abs(H2+R2+2*HR-S*(S+1)),'p_H':(H2-.75*nh)/(nh*(nh-1)/4),'m_R_uniform2':R2/nr**2,'m_N2':neel,'legacy_Neel_FGamma':nr**2*neel/nh,'hub_S_Gamma':hs[0],'hub_peak':peak,'hub_peak_q':json.dumps(peakq),'hub_structure_spin_summed':hs}

def checks(b,psi,G):
 groups={'H':np.arange(2,b.n,3),'R':np.array([i for i in range(b.n) if i%3!=2])};out={}
 for name,ix in groups.items():
  pp=[(int(a),int(c),2.) for k,a in enumerate(ix) for c in ix[k+1:]];v=full_pair_action(b.states,b.lookup,psi,pp,.75*len(ix));out[name+'2_operator_vs_G']=float(abs(np.vdot(psi,v).real-G[np.ix_(ix,ix)].sum()))
 return out
