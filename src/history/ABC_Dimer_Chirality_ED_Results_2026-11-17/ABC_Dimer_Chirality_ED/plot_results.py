import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/abc_dynamics_matplotlib')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/abc_dynamics_cache')
os.environ.setdefault('MPLCONFIGDIR','/tmp/abc_dynamics_matplotlib')
os.environ.setdefault('XDG_CACHE_HOME','/tmp/abc_dynamics_cache')
import json,numpy as np,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from experiment import ROOT,geometry_weights
from group_poles import unpack
from ensemble import fourier,unfold
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'figure.constrained_layout.use':True})
R=json.loads((ROOT/'responses.json').read_text());Z=np.load(ROOT/'residues.npz');T=np.eye(6,dtype=complex)[[2,3,0,1,4,5]];T[4:]*=-1

def save(fig,name):
 fig.savefig(ROOT/(name+'.png'),dpi=300);fig.savefig(ROOT/(name+'.svg'));plt.close(fig)
def quads(name,u):
 data=[]
 if name=='isotropic_G_pair_average':spec=[('j1_G',1)]
 elif name=='isotropic_Q_pair_average':spec=[('j1_Q',.5),('j1_Qprime',.5)]
 elif name=='isotropic_K_pair_average':spec=[('j1_K',.5),('j1_KprimeSource',.5)]
 else:spec=[(name,1)]
 for n,factor in spec:
  p=ROOT/(n+'_quadrature.npz')
  if not p.exists():continue
  z=np.load(p);A=z['amplitudes'];w=abs(A@u)**2
  if name=='isotropic_G_pair_average':w=(w+abs(A.conj()@T@u)**2)/2
  data.extend(zip(z['omega'],factor*w))
 return sorted(data)
def poles(name,u):return sorted((p['omega'],float(np.vdot(u,unpack(p['residue'])@u).real)) for p in R[name]['poles'])
def sticks(ax,items,color,label,marker='o'):
 if not items:return
 x,y=np.array(items).T;keep=y>1e-12;x=x[keep];y=y[keep];ax.vlines(x,0,y,color=color,lw=1.1);ax.scatter(x,y,s=14,c=color,label=label,marker=marker,zorder=4)
def draw(name,kind,ax):
 u=Z[name+'_'+kind+'_reference_mode'];sticks(ax,quads(name,u),'0.78','Approximate Krylov quadrature','x');sticks(ax,poles(name,u),'#245c9f','Verified eigenspace residues');ax.set_xlim(0,6);ax.set_ylim(bottom=0);ax.set_xlabel(r'$\omega/J$');ax.set_ylabel('Discrete spectral weight');ax.legend(fontsize=8)
 if kind=='chi':title=r'Scalar chirality, physical $\Gamma$'
 else:title=r'Dimer, physical $Q=(8/9,7/9)$'
 ax.set_title(title)
for name,kind,title in [('isotropic_G_pair_average','chi','chirality_Gamma'),('isotropic_Q_pair_average','D','dimer_finite_Q')]:
 if name not in R:continue
 fig,axes=plt.subplots(2,1,figsize=(7,6),height_ratios=[2,1]);draw(name,kind,axes[0]);u=Z[name+'_'+kind+'_reference_mode'];m=R[name].get('moments',R[name]);m0=float(np.vdot(u,unpack(m['M0'])@u).real);ps=poles(name,u)
 if ps:
  x,y=np.array(ps).T;axes[1].step(np.r_[0,x,6],np.r_[0,np.cumsum(y)/m0,np.sum(y)/m0],where='post',color='#245c9f',label='Verified accumulated weight / exact M0')
 q=quads(name,u)
 if q:
  x,y=np.array(q).T;axes[1].step(np.r_[0,x,6],np.r_[0,np.cumsum(y)/m0,np.sum(y)/m0],where='post',color='0.65',label='Quadrature estimate')
 axes[1].axhline(1,color='k',ls=':',lw=.8);axes[1].set(xlim=(0,6),ylim=(0,1.06),xlabel=r'$\omega/J$',ylabel='Weight fraction');axes[1].legend(fontsize=8);fig.suptitle('J′/J = 1; orbital ground-pair average; no broadening');save(fig,title)
# Full complex cross residues, including the orthogonal chirality face channel.
if 'isotropic_Q_pair_average' in R:
 name='isotropic_Q_pair_average';uD=Z[name+'_D_reference_mode'];raw=np.load(ROOT/'inputs/raw_observables_j1.npz');S=unfold(fourier(raw['ensemble_full_connected'][108:,108:]),'chi')[0];entry=min(S,key=lambda x:np.linalg.norm(x['Q']-[8/9,7/9]));U=entry['eigenvectors'];fig,axes=plt.subplots(2,1,figsize=(7,5),sharex=True)
 for mode in [0,1]:
  u=np.zeros(6,complex);u[4:6]=U[:,mode];xs=[];ys=[]
  for p in R[name]['poles']:xs.append(p['omega']);ys.append(np.vdot(uD,unpack(p['residue'])@u))
  for ax,part,label in [(axes[0],np.real,'Re'),(axes[1],np.imag,'Im')]:
   vals=part(ys);ax.vlines(xs,0,vals,color=['#cb6c37','#245c9f'][mode],lw=1);ax.scatter(xs,vals,s=12,c=['#cb6c37','#245c9f'][mode],label=f'χ face mode {mode}');ax.axhline(0,color='k',lw=.5);ax.set_ylabel(label+r' $R_{D,\chi}$');ax.legend(fontsize=8)
 axes[1].set(xlim=(0,3.5),xlabel=r'$\omega/J$');fig.suptitle('Verified complex cross residues; fixed CCW/geometric phase convention');save(fig,'complex_cross_residues')
# Identical physical Gamma probe, embedded into all six cell faces off isotropy.
fig,axes=plt.subplots(1,3,figsize=(11,3.5),sharex=True,sharey=True);anyplot=False
for ax,(name,label) in zip(axes,[('j07_G','0.7'),('j09_G','0.9'),('isotropic_G_pair_average','1.0')]):
 if name not in R:ax.set_visible(False);continue
 anyplot=True;u=Z[name+'_chi_reference_mode'];sticks(ax,quads(name,u),'0.8','Approximate Krylov nodes','x');sticks(ax,poles(name,u),'#245c9f','Verified poles');ax.set(title='J′/J = '+label,xlim=(0,5),xlabel=r'$\omega/J$');ax.set_ylim(bottom=0);ax.legend(fontsize=8,loc='upper right')
axes[0].set_ylabel('Uniform CCW chirality weight');fig.suptitle('Same normalized Γ chirality probe; no artificial linewidth');
if anyplot:save(fig,'chirality_coupling_evolution')
else:plt.close(fig)
# Elastic contribution is explicitly separate from all positive-frequency figures.
labels=[];vals=[]
for name,a,label in [('j1_G',4,'χ Γ, orbital pair'),('j1_K',0,'D K, selected K source'),('j1_KprimeSource',0,'D K, selected K′ source')]:
 p=ROOT/(name+'_moments.json')
 if p.exists():d=json.loads(p.read_text());M=unpack(d['M0_elastic'])[a:a+2,a:a+2];labels.append(label);vals.append(float(np.linalg.eigvalsh(M)[-1]))
if labels:
 fig,ax=plt.subplots(figsize=(7,3));ax.bar(np.arange(len(vals)),vals,color='#9b6b99');ax.set_xticks(np.arange(len(vals)),labels,rotation=10);ax.set(ylabel='Largest zero-energy face-mode weight',title='Ground-manifold contribution (ω = 0 only)');save(fig,'zero_energy_weights')
print('Plots written; no broadening used.')

if 'isotropic_K_pair_average' in R:
 name='isotropic_K_pair_average';fig,axes=plt.subplots(1,2,figsize=(10,3.5))
 for ax,kind in zip(axes,['D','chi']):
  u=Z[name+'_'+kind+'_reference_mode'];sticks(ax,quads(name,u),'0.8','Approximate quadrature','x');sticks(ax,poles(name,u),'#245c9f','Verified poles');ax.set(xlim=(0,6),xlabel=r'$\omega/J$',ylabel='Positive-energy weight',title=kind+' at physical K');ax.set_ylim(bottom=0);ax.legend(fontsize=8)
 fig.suptitle('Ground-pair deflated; K′ partner supplied by exact spatial symmetry');save(fig,'physical_K_spectra')
if 'j09_q' in R:
 fig,axes=plt.subplots(1,2,figsize=(10,3.5),sharex=True,sharey=True)
 for ax,name,label in zip(axes,['j09_q','isotropic_Q_pair_average'],['0.9','1.0']):
  u=Z[name+'_D_reference_mode'];sticks(ax,quads(name,u),'0.8','Approximate quadrature','x');sticks(ax,poles(name,u),'#245c9f','Verified poles');ax.set(xlim=(0,4),xlabel=r'$\omega/J$',title='J′/J = '+label);ax.set_ylim(bottom=0);ax.legend(fontsize=8)
 axes[0].set_ylabel('Same finite-Q dimer-probe weight');fig.suptitle('Identical physical form factor; full cell symmetry retained at 0.9');save(fig,'dimer_coupling_evolution')
if (ROOT/'ground_projection_mode_checks.json').exists():
 checks=json.loads((ROOT/'ground_projection_mode_checks.json').read_text())['records'];values=[];labels=[]
 for kind,Q,label in [('chi',[0,0],'χ Γ, subleading face mode'),('D',[1/3,2/3],'D K, subleading face mode')]:
  row=next(r for r in checks if r['kind']==kind and np.linalg.norm(np.array(r['physical_Q'])-Q)<1e-10);values.append(row['orbital_average_eigenvalues'][-1]);labels.append(label)
 fig,ax=plt.subplots(figsize=(7,3));ax.bar(range(2),values,color='#9b6b99');ax.set_xticks(range(2),labels);ax.set(ylabel='Orbital-averaged zero-energy weight',title='Ground-manifold contribution only (ω = 0)');save(fig,'zero_energy_weights')
