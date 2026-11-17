#!/usr/bin/env python3
import csv
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=Path(__file__).parent
with open(p/'triangle_27spin_existing.csv') as f:
    rows=list(csv.DictReader(f))
ground=sorted((r for r in rows if r['branch']=='ground'),key=lambda r:float(r['Jp_over_J']))
x=np.array([float(r['Jp_over_J']) for r in ground])
get=lambda key:np.array([float(r[key]) for r in ground])
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.labelsize':11,'axes.titlesize':11,'xtick.labelsize':9,'ytick.labelsize':9,'legend.fontsize':9,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(1,3,figsize=(13.1,3.65),constrained_layout=True)
ax=axes[0]
ax.scatter(x,get('P_triangle_doublet'),s=37,marker='o',label=r'Local doublet, $S_{\triangle}=1/2$',color='#17688a',zorder=3)
ax.set_ylim(.73,.91);ax.set_xlim(.16,1.04);ax.set_ylabel(r'$P_{1/2}$');ax.set_xlabel("$J'/J$")
ax.set_title('(a) Triangle doublet probability',loc='left',fontweight='bold')
ax.axhline(.75,color='0.65',lw=1,ls=':',label=r'Classical 120$^\circ$ product state')
ax.legend(loc='upper left',frameon=False,fontsize=8)
ax=axes[1]
ax.scatter(x,get('P_singlet_AB'),s=35,marker='o',color='#136d9c',label='$AB$ singlet probability')
ax.scatter(x,get('P_singlet_AC'),s=39,marker='s',color='#d26d28',label='$AC,BC$ singlet probability')
ax.set_xlim(.16,1.04);ax.set_ylim(.23,.65);ax.set_xlabel("$J'/J$");ax.set_ylabel(r'$P_{s}(ij)=1/4-\langle\mathbf{S}_i\cdot\mathbf{S}_j\rangle$')
ax.set_title('(b) Bond singlet weights',loc='left',fontweight='bold')
ax.legend(loc='center right',frameon=False,fontsize=8)
ax=axes[2]
ax.scatter(x,get('D_abs'),s=38,marker='D',color='#7548a2',label=r'$|\langle\mathcal{D}_{\triangle}\rangle|$')
ax.set_xlim(.16,1.04);ax.set_ylim(-.025,.42);ax.set_xlabel("$J'/J$");ax.set_ylabel(r'$|\langle\mathcal{D}_{\triangle}\rangle|$')
ax.set_title('(c) Three-bond inequivalence',loc='left',fontweight='bold')
ax.annotate('Equal bonds at $J\'=J$',xy=(1,0),xytext=(.63,.145),arrowprops={'arrowstyle':'->','lw':1},fontsize=9)
for ax in axes:
    for start,end in [(0.25625,0.2609375),(0.3875,0.3921875)]:
        ax.axvspan(start,end,color='#bbc6cc',alpha=.45,zorder=0)
    ax.grid(alpha=.15,axis='y',lw=.6)
fig.suptitle('27-spin ED: local triangular correlations (finite-size sectors)',size=14,y=1.12)
fig.savefig(p/'triangle_local_ED_figure.png',dpi=220,bbox_inches='tight')
fig.savefig(p/'triangle_local_ED_figure.svg',bbox_inches='tight')
print('Created local-triangle figures.')
