import pathlib,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=pathlib.Path(__file__).resolve().parent
fig,ax=plt.subplots(figsize=(7,4.5),layout='constrained')
for name in ['Gamma_primitive','K_primitive','cell11_fold1','cell11_fold2','cell11_fold3','other']:
 r=json.loads((root/(name+'.json')).read_text());h=r['history'];ax.semilogy([x['iteration'] for x in h],[x['recurrence_residual'] for x in h],'.-',label=name)
 ax.plot(r['iterations'],r['residual'],'ko',markersize=3)
ax.axhline(1e-8,color='gray',linestyle='--',label='direct-residual tolerance')
ax.set(xlabel='Lanczos iteration',ylabel='Residual norm / J',title='27-site Heisenberg model: convergence of six momentum classes')
ax.grid(alpha=.25);ax.legend(ncols=2);fig.savefig(root/'convergence.png',dpi=180);fig.savefig(root/'convergence.svg')
