import os,pathlib,json
os.environ.setdefault('MPLCONFIGDIR',str(pathlib.Path(__file__).resolve().parent/'.mpl'))
os.environ.setdefault('XDG_CACHE_HOME',str(pathlib.Path(__file__).resolve().parent/'.cache'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt,numpy as np
ROOT=pathlib.Path(__file__).resolve().parent
records=sorted([json.loads(p.read_text()) for p in ROOT.glob('static_j*.json')],key=lambda d:d['summary']['Jp_over_J']);s=[d['summary'] for d in records];x=[d['Jp_over_J'] for d in s]
fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
plots=[([('rim_neel_m2','Rim Neel squared amplitude')],'m_N squared'),([('hub_S_Gamma','hub k=(0,0)'),('hub_S_cell01','hub k=(0,1/3)'),('hub_S_cell12','hub k=(1/3,2/3)')],'S_h(k)'),([('rim_hub_uniform_Gamma','uniform rim-hub cross')],'C_mh(Gamma)'),([('rim_NN_mean','rim NN'),('hub_rim_bond_mean','hub-rim bonds'),('hub_hub_all_distinct_mean','all distinct hub pairs')],'Mean spin dot product'),([('physical_K_120_intensity','physical K diagnostic')],'Physical K equal-time intensity'),([('rim_neel_normalized_mode_eigenvalue','rim Neel'),('uniform_hub_mode_lower_eigenvalue','uniform/hub lower'),('uniform_hub_mode_upper_eigenvalue','uniform/hub upper')],'Gamma correlation eigenvalue')]
for ax,(series,ylabel) in zip(axes.flat,plots):
 for key,label in series:ax.plot(x,[d[key] for d in s],'o-',label=label)
 ax.set(xlabel="J'/J",ylabel=ylabel);ax.grid(alpha=.25);ax.legend(fontsize=8)
fig.suptitle('27 spins: static correlations, verified ground states\nSU(2) and exact point-group averages; lines guide the eye')
fig.savefig(ROOT/'static_evolution.png',dpi=180);fig.savefig(ROOT/'static_evolution.svg');plt.close(fig)
fig,axes=plt.subplots(2,len(records),figsize=(4*len(records),7),squeeze=False,layout='constrained')
maxabs=max(np.max(abs(np.array(d['C_SU2_pointgroup']['(0, 0)']['real']))) for d in records)
for column,d in enumerate(records):
 jp=d['summary']['Jp_over_J'];C=np.array(d['C_SU2_pointgroup']['(0, 0)']['real']);im=axes[0,column].imshow(C,vmin=-maxabs,vmax=maxabs,cmap='coolwarm')
 axes[0,column].set(xticks=range(3),xticklabels=['A','B','H'],yticks=range(3),yticklabels=['A','B','H'],title=f"J'/J={jp:g}: C(Gamma)")
 for a in range(3):
  for b in range(3):axes[0,column].text(b,a,f'{C[a,b]:.3f}',ha='center',va='center',fontsize=9)
 sh=np.array([[d['C_SU2_pointgroup'][str((a,b))]['real'][2][2] for a in range(3)] for b in range(3)])
 im2=axes[1,column].imshow(sh,origin='lower',cmap='viridis');axes[1,column].set(xticks=range(3),xticklabels=['0','1/3','2/3'],yticks=range(3),yticklabels=['0','1/3','2/3'],xlabel='cell k1',ylabel='cell k2',title='Hub structure factor');fig.colorbar(im2,ax=axes[1,column],shrink=.8)
fig.colorbar(im,ax=list(axes[0]),shrink=.8,label='Scalar correlation matrix element');fig.savefig(ROOT/'static_matrices.png',dpi=180);fig.savefig(ROOT/'static_matrices.svg');plt.close(fig)
# Momentum-resolved mixed eigenmodes; degeneracy colors use invariant projector weights.
fig,axes=plt.subplots(1,3,figsize=(11,4),layout='constrained')
for ax,celln in zip(axes,[(0,0),(0,1),(1,2)]):
 for mode_index in range(3):
  yy=[]
  for d in records:
   Ck=d['C_SU2_pointgroup'][str(celln)];M=np.array(Ck['real'])+1j*np.array(Ck['imag']);yy.append(np.linalg.eigvalsh(M)[mode_index])
  ax.plot(x,yy,color='.75',linewidth=1)
 for d in records:
  for group in d['degenerate_correlation_eigenspaces']:
   if group['cell_n']!=str(celln):continue
   scatter=ax.scatter(d['summary']['Jp_over_J'],group['eigenvalue_mean'],c=[group['hub_projector_trace']/group['dimension']],vmin=0,vmax=1,cmap='viridis',s=45*group['dimension'],edgecolors='black',linewidths=.4)
 ax.set(xlabel="J'/J",ylabel='Correlation eigenvalue',title=f'cell k=({celln[0]}/3,{celln[1]}/3)',yscale='log');ax.grid(alpha=.25)
fig.colorbar(scatter,ax=list(axes),label='Hub weight per eigenspace dimension');fig.suptitle('Mixed sublattice correlation modes\nDegenerate spaces use projector-averaged hub weights; marker area shows multiplicity')
fig.savefig(ROOT/'hub_modes.png',dpi=180);fig.savefig(ROOT/'hub_modes.svg');plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(11,3.8),layout='constrained')
axes[0].plot(x,[d['S_total'] for d in s],'o--');axes[0].set(xlabel="J'/J",ylabel='Verified total spin S',yticks=[.5,1.5,2.5],title='Finite-cluster ground multiplet')
axes[1].plot(x,[d['ferromagnetic_m2'] for d in s],'o--');axes[1].set(xlabel="J'/J",ylabel='S(S+1) / 27 squared',title='Total-spin squared indicator')
for a,label,marker in [(0,'A','o'),(1,'B','x'),(2,'H','s')]:axes[2].plot(x,[d['single_state_sublattice_Sz_per_site'][a] for d in s],marker+'--',label=label)
axes[2].axhline(0,color='gray',linewidth=.7);axes[2].set(xlabel="J'/J",ylabel='Selected-state mean S_z per site',title='Selected Sz=+1/2 representative');axes[2].legend()
for ax in axes:ax.grid(alpha=.25)
fig.suptitle('Spin quantum numbers and partial sublattice compensation\nFinite size; dashed lines guide the eye, not a phase boundary')
fig.savefig(ROOT/'total_spin_and_polarization.png',dpi=180);fig.savefig(ROOT/'total_spin_and_polarization.svg');plt.close(fig)
