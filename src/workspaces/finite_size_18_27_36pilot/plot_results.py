import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ.setdefault('MPLCONFIGDIR','/tmp/hub_size_mpl');os.environ.setdefault('XDG_CACHE_HOME','/tmp/hub_size_cache')
import numpy as np,csv,json,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry import ROOT,make
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'figure.constrained_layout.use':True})
rows=list(csv.DictReader((ROOT/'size_comparison.csv').open()));styles={'K18':('#2464a0','o','18 spins: k₀-compatible rectangle'),'C18':('#d07a30','s','18 spins: more isotropic, k₀ absent'),'N27':('#39884c','^','27 spins: verified 3×3')}
def save(fig,name):
 for ext in ['png','svg','pdf']:fig.savefig(ROOT/(name+'.'+ext),dpi=300)
 plt.close(fig)
def curve(ax,field,ylabel):
 for name,(c,m,l) in styles.items():
  rr=sorted([r for r in rows if r['shape']==name and float(r['r'])<=.70000001],key=lambda r:float(r['r']));ax.plot([float(r['r']) for r in rr],[float(r[field]) for r in rr],color=c,marker=m,ms=3,lw=1,ls='None' if name=='N27' else '-',label=l)
 ax.set(xlabel="r = J′/J",ylabel=ylabel,xlim=(.09,.71));ax.grid(alpha=.15)
def marks(ax):
 for r in csv.DictReader((ROOT/'crossing_brackets_18.csv').open()):ax.axvspan(float(r['r_lo']),float(r['r_hi']),color=styles[r['shape']][0],alpha=.12)
 for r in csv.DictReader((ROOT/'reference27/crossing_brackets_27.csv').open()):ax.axvspan(float(r['a']),float(r['b']),color=styles['N27'][0],alpha=.12)
fig,axes=plt.subplots(2,2,figsize=(10,7))
for ax,f,l in zip(axes.flat,['p_H','m_N2','m_R_uniform2','m_tot'],[r'Hub pair polarization $p_H$',r'Rim Néel $m_N^2$',r'Rim uniform $R^2/N_R^2$',r'$S_{tot}/N$']):curve(ax,f,l);marks(ax)
axes[0,0].legend(fontsize=8);fig.suptitle('Finite-size hub polarization and rim compensation; bands are measured crossing brackets');save(fig,'finite_size_observables')
fig,axes=plt.subplots(1,2,figsize=(10,3.8));curve(axes[0],'S',r'$S_{tot}$');curve(axes[1],'m_tot',r'$S_{tot}/N$');marks(axes[0]);marks(axes[1]);axes[0].legend(fontsize=8);fig.suptitle('Integer versus half-integer spin quantization\nGround momenta: K18 (0,½)→Γ; C18 Γ; 27 k₀→Γ');save(fig,'ground_state_spin')
fig,axes=plt.subplots(1,2,figsize=(10,3.8));curve(axes[0],'H2',r'$\langle S_H^2\rangle$');curve(axes[1],'R2',r'$\langle S_R^2\rangle$');axes[0].legend(fontsize=8);fig.suptitle('Subsystem moments are expectation values, not conserved subsystem spin labels');save(fig,'subsystem_compensation')
fig,axes=plt.subplots(1,2,figsize=(10,3.8));curve(axes[0],'hub_S_Gamma',r'$S_h(\Gamma)$, spin summed');curve(axes[1],'hub_peak',r'$\max_q S_h(q)$, spin summed');axes[0].legend(fontsize=8);fig.suptitle('Hub structure factors on each torus’s allowed cell momenta');save(fig,'hub_structure_evolution')
# Physical-site geometry in a reduced fundamental parallelogram.
fig,axes=plt.subplots(2,2,figsize=(10,9));cart=np.array([[1,0],[.5,np.sqrt(3)/2]]);A=np.array([[1,1],[2,-1]])
for ax,(label,T) in zip(axes.flat,[('18 k₀-compatible',[[3,0],[2,2]]),('18 control',[[1,2],[0,6]]),('27',[[3,0],[0,3]]),('36 (geometry/pilot only)',[[2,2],[-2,4]])]):
 g=make(T);P=g['reduced_T']@A;V=P@cart;frac=np.mod(g['coords']@np.linalg.inv(P),1);pos=frac@V
 for a,b,t in g['bonds']:
  options=[pos[b]+u*V[0]+v*V[1] for u in [-1,0,1] for v in [-1,0,1]];j=min(options,key=lambda x:np.linalg.norm(x-pos[a]));wrapped=np.linalg.norm(j-pos[b])>1e-6;ax.plot([pos[a,0],j[0]],[pos[a,1],j[1]],color='0.72' if t else '0.35',lw=.65,ls='--' if wrapped else '-',zorder=1)
 boundary=np.array([[0,0],V[0],V.sum(axis=0),V[1],[0,0]]);ax.plot(boundary[:,0],boundary[:,1],color='k',lw=1)
 for s,c,m in [(0,'#2464a0','o'),(1,'#d07a30','s'),(2,'#824b9a','^')]:ii=np.arange(s,g['N'],3);ax.scatter(pos[ii,0],pos[ii,1],s=25,c=c,marker=m,label=['A','B','H'][s],zorder=3)
 ax.set_aspect('equal');ax.set_title(label+f"; shortest winding {g['shortest_winding_site_units']:.3f}");ax.set_xlabel('x / site spacing');ax.set_ylabel('y / site spacing')
axes[0,0].legend(ncol=3,fontsize=8);fig.suptitle('Exact nearest-neighbor graphs; dashed bonds cross a periodic boundary');save(fig,'finite_size_geometries')
# All finite-q hub intensities at representative couplings, keeping each distinct q set.
hr=list(csv.DictReader((ROOT/'hub_structure_factors.csv').open()));fig,axes=plt.subplots(1,3,figsize=(12,3.6))
for ax,r in zip(axes,[.2,.5,.7]):
 for name,(color,marker,label) in styles.items():
  qq=[x for x in hr if x['shape']==name and abs(float(x['r'])-r)<1e-8];qq.sort(key=lambda x:(float(x['k1']),float(x['k2'])))
  if qq:ax.plot(range(len(qq)),[float(x['S_h_spin_summed']) for x in qq],color=color,marker=marker,ms=4,label=label)
 ax.set(title=f'r={r}',xlabel='allowed-q index (see CSV)',ylabel='spin-summed hub structure factor');ax.grid(alpha=.15)
axes[0].legend(fontsize=7);fig.suptitle('Different allowed-q sets are not identified with one another');save(fig,'hub_q_profiles')
print('PNG/SVG/PDF figures written.')
