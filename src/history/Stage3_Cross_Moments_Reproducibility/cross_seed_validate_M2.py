import os,time
os.environ['OPENBLAS_NUM_THREADS']='1';os.environ['OMP_NUM_THREADS']='2'
import numpy as np
from ed import Model
from cyclic_model import CyclicModel
from observables import spin_seed,convert_basis,cell_q
T=time.monotonic()
def msg(s):print(f'{time.monotonic()-T:8.2f}s {s}',flush=True)
z=np.load('j0.2_n12_r0.npz');msg('read checkpoint '+str(z.files))
src=CyclicModel(jp=.2,cell_n=(1,2),angular=0,threads=2);msg('source cyclic basis '+str(src.nb))
source=Model(primitive=False,jp=.2,threads=2);source.sector(src.q);msg('source cell basis '+str(source.nb))
v=convert_basis(src,source,z['vector']);msg('converted source norm '+str(np.linalg.norm(v)));del src
target=Model(primitive=False,jp=.2,threads=2);target.sector(source.q+cell_q((0,1/3)));msg('target cell basis '+str(target.nb)+' q '+str(target.q))
channels=['n','m','H','Phi']
b={c:spin_seed(source,target,v,(0,1/3),c,kind='z') for c in channels};msg('created seeds')
G=np.array([[np.vdot(b[a],b[c]) for c in channels] for a in channels]);msg('M0 from direct seeds: '+str(G))
E0=float(z['energy'])
Hb={c:target.mv(b[c])-E0*b[c] for c in channels};msg('all four (H-E0)|b> vectors computed')
M1=np.array([[np.vdot(b[a],Hb[c]) for c in channels] for a in channels]);
M2=np.array([[np.vdot(Hb[a],Hb[c]) for c in channels] for a in channels]);
msg('M1 matrix='+str(np.round(M1,8)))
msg('M2 matrix='+str(np.round(M2,8)))
for c in ['n','m','Phi']:
 i=channels.index(c);j=channels.index('H')
 msg(f'cross {c}-H: M0={G[i,j]:.12g}, M1={M1[i,j]:.12g}, M2={M2[i,j]:.12g}')
i=channels.index('H')
mean=M1[i,i].real/G[i,i].real;variance=M2[i,i].real/G[i,i].real-mean**2
msg(f'Hub mean={mean}, spectral std dev={np.sqrt(variance)}, M2={M2[i,i]}')
np.savez_compressed('/mnt/data/stage3_cross_moments_direct_27spin.npz',channels=channels,M0=G,M1=M1,M2=M2,energy=E0,probe_k=[0,1/3],Jp=.2)
with open('/mnt/data/stage3_cross_moments_direct_27spin.csv','w') as f:
 f.write('channel_a,channel_b,M0_re,M0_im,M1_re,M1_im,M2_re,M2_im\n')
 for i,a in enumerate(channels):
  for j,c in enumerate(channels):
   f.write(f'{a},{c},{G[i,j].real:.14g},{G[i,j].imag:.14g},{M1[i,j].real:.14g},{M1[i,j].imag:.14g},{M2[i,j].real:.14g},{M2[i,j].imag:.14g}\n')
msg('finished')
