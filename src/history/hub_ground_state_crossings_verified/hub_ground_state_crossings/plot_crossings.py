import os,pathlib,json,csv,numpy as np
ROOT=pathlib.Path(__file__).resolve().parent
os.environ.setdefault('MPLCONFIGDIR',str(ROOT/'.mpl'));os.environ.setdefault('XDG_CACHE_HOME',str(ROOT/'.cache'))
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.hashsalt']='stuffed-honeycomb-27-static-crossings'
import matplotlib.pyplot as plt

def rows(name):
 with (ROOT/name).open() as f:return list(csv.DictReader(f))
g=rows('ground_samples.csv');lev=rows('all_levels.csv');obs=rows('static_observables.csv');interval=json.loads((ROOT/'crossing_intervals.json').read_text()) if (ROOT/'crossing_intervals.json').exists() else {'crossings':[]}
fig,axes=plt.subplots(2,1,figsize=(9,7),layout='constrained')
for ax,key,ylabel in [(axes[0],'spin','Verified total spin S'),(axes[1],'hub_rim_bond_mean','Mean hub-rim spin dot product')]:
 ax.plot([float(d['Jp_over_J']) for d in g],[float(d[key]) for d in g],'o--',markersize=4);ax.set(xlabel="J'/J",ylabel=ylabel);ax.grid(alpha=.25)
 for c in interval['crossings']:ax.axvspan(c['a'],c['b'],color='gray',alpha=.25)
fig.suptitle('27 spins: verified ground-sector samples\nShading: resolved crossing brackets; lines guide the eye')
fig.savefig(ROOT/'ground_sector_tracking.png',dpi=180);fig.savefig(ROOT/'ground_sector_tracking.svg',metadata={'Date':None});plt.close(fig)
# Fixed spin/spatial labels are separate curves: no interpolation across a spin change.
fig,ax=plt.subplots(figsize=(10,5),layout='constrained');base={float(d['Jp_over_J']):float(d['energy_J']) for d in lev if d['cell_n']=='[0, 0]' and d['character']=='0' and d['physical_root_index']=='0' and d['verified']=='True'}
families={}
for d in lev:
 if d['verified']!='True' or d['cell_n'] not in ['[0, 0]','[0, 1]','[1, 2]']:continue
 key=(d['cell_n'],int(d['character']),float(d['spin']),d['reflection_parity']);j=float(d['Jp_over_J']);E=float(d['energy_J'])
 if int(d['character']) not in ([0,3] if d['cell_n']=='[0, 0]' else [0]):continue
 families.setdefault(key,{})[j]=min(E,families.get(key,{}).get(j,np.inf))
for (n,r,S,p),values in families.items():
 js=sorted(j for j in values if j in base);ax.plot(js,[values[j]-base[j] for j in js],'.-',linewidth=1,label=f'n={n}/3, chi={r}, S={S:g}, F={p or "?"}')
for c in interval['crossings']:ax.axvspan(c['a'],c['b'],color='gray',alpha=.2)
ax.axhline(0,color='black',linewidth=.6);ax.set(xlabel="J'/J",ylabel='E - lowest E(Gamma, character 0) [J]',title='Tracked quantum-number branches at measured points');ax.grid(alpha=.25);ax.legend(fontsize=7,ncols=2,loc='upper left');fig.savefig(ROOT/'competing_energy_branches.png',dpi=180);fig.savefig(ROOT/'competing_energy_branches.svg',metadata={'Date':None});plt.close(fig)
fig,axes=plt.subplots(1,3,figsize=(12,4),layout='constrained');groundobs=[d for d in obs if d['role'].startswith('GROUND')]
for ax,keys,ylabel in [(axes[0],[('rim_neel_m2','Rim Neel squared amplitude')],'m_N squared'),(axes[1],[('hub_S_Gamma','Gamma'),('hub_S_cell01','(0,1/3)'),('hub_S_cell12','(1/3,2/3)')],'Hub S_h(k)'),(axes[2],[('rim_hub_uniform_Gamma','Uniform rim-hub Gamma cross')],'C_mh(Gamma)')]:
 for key,label in keys:ax.plot([float(d['Jp']) for d in groundobs],[float(d[key]) for d in groundobs],'o-',label=label,markersize=4)
 for c in interval['crossings']:ax.axvspan(c['a'],c['b'],color='gray',alpha=.2)
 ax.set(xlabel="J'/J",ylabel=ylabel);ax.grid(alpha=.25);ax.legend(fontsize=8)
fig.suptitle('Static scalar correlations: spin-multiplet and point-group averages\nFinite size; lines guide the eye');fig.savefig(ROOT/'correlations_across_crossings.png',dpi=180);fig.savefig(ROOT/'correlations_across_crossings.svg',metadata={'Date':None});plt.close(fig)
if interval['crossings']:
 fig,axes=plt.subplots(1,len(interval['crossings']),figsize=(5*len(interval['crossings']),3.7),squeeze=False,layout='constrained')
 for ax,c in zip(axes.flat,interval['crossings']):
  if 'gap_right_minus_left_at_a' in c and 'gap_right_minus_left_at_b' in c:ax.plot([c['a'],c['b']],[c['gap_right_minus_left_at_a'],c['gap_right_minus_left_at_b']],'o-')
  ax.axhline(0,color='black',linewidth=.7);ax.axvspan(c['a'],c['b'],color='gray',alpha=.2);ax.set(xlabel="J'/J",ylabel='E(right branch) - E(left branch) [J]',title=f"Bracket [{c['a']:.6f}, {c['b']:.6f}]");ax.grid(alpha=.25)
 fig.suptitle('Measured competing-branch splittings, not a complete excitation gap');fig.savefig(ROOT/'crossing_energy_gaps.png',dpi=180);fig.savefig(ROOT/'crossing_energy_gaps.svg',metadata={'Date':None});plt.close(fig)
