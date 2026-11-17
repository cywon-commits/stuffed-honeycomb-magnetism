import numpy as np, math, time, os, sys, gc
from numba import njit, prange
import numba
from scipy.sparse.linalg import LinearOperator, eigsh, ArpackNoConvergence
from math import comb

N=27; KUP=14; DIM=comb(27,14)
J=1.0; Jp=1.0
numba.set_num_threads(min(3, numba.config.NUMBA_NUM_THREADS))

def log(msg):
    print(f'[{time.strftime("%H:%M:%S")}] {msg}', flush=True)

# Site convention: site=3*(3*x+y)+s, s=0(A),1(B),2(H).
def sid(x,y,s): return 3*((x%3)*3+(y%3))+s

# bonds
bonds=[]
for x in range(3):
  for y in range(3):
    # honeycomb rim edges A_R-B_R, A_R-B_{R-a1}, A_R-B_{R-a2}
    a=sid(x,y,0)
    bonds.append((a,sid(x,y,1),J))
    bonds.append((a,sid(x-1,y,1),J))
    bonds.append((a,sid(x,y-1,1),J))
# hub-to-six-vertices of the corresponding hexagon
for x in range(3):
  for y in range(3):
    h=sid(x,y,2)
    neigh=[sid(x,y,0), sid(x,y,1), sid(x+1,y,0),
           sid(x+1,y-1,1), sid(x+1,y-1,0), sid(x,y-1,1)]
    for v in neigh: bonds.append((h,v,Jp))
bonds=np.asarray(bonds,dtype=np.float64)
assert bonds.shape==(81,3), bonds.shape
assert len({(min(int(a),int(b)),max(int(a),int(b))) for a,b,_ in bonds})==81
log(f'Geometry: {len(bonds)} bonds (27 rim-rim, 54 hub-rim), J=J\'=1; Hilbert dimension={DIM:,}')

# Translation permutations on the 27 bits.
trans_lo=np.zeros((9,1<<13),dtype=np.uint32)
trans_hi=np.zeros((9,1<<14),dtype=np.uint32)
trans_xy=[]
for dx in range(3):
  for dy in range(3):
    trans_xy.append((dx,dy))
    for m in range(1<<13):
      out=0
      for p in range(13):
        if (m>>p)&1:
          c=p//3; s=p%3; x=c//3; y=c%3
          out |= 1 << sid(x+dx,y+dy,s)
      trans_lo[dx*3+dy,m]=out
    for m in range(1<<14):
      out=0
      for pp in range(14):
        if (m>>pp)&1:
          p=13+pp; c=p//3; s=p%3; x=c//3; y=c%3
          out |= 1 << sid(x+dx,y+dy,s)
      trans_hi[dx*3+dy,m]=out
log('Translation lookup tables constructed.')

# Combinadic rank lookup: rank(mask)=sum_i C(p_i,i+1), with p_i increasing.
poplo=np.zeros(1<<13,dtype=np.uint8)
ranklo=np.zeros(1<<13,dtype=np.int32)
for m in range(1<<13):
  cnt=0; rr=0
  for p in range(13):
    if (m>>p)&1:
      cnt+=1; rr+=comb(p,cnt)
  poplo[m]=cnt; ranklo[m]=rr
rankhi=np.zeros((15,1<<14),dtype=np.int32)
for kl in range(15):
  for m in range(1<<14):
    cnt=0; rr=0
    for t in range(14):
      if (m>>t)&1:
        cnt+=1; n=13+t; k=kl+cnt
        if k<=n: rr+=comb(n,k)
    rankhi[kl,m]=rr
log('Combinadic rank lookup tables constructed.')

@njit
def rank_mask(mask, poplo, ranklo, rankhi):
    mask=np.int64(mask)
    lo=mask & 8191
    hi=mask >> 13
    return ranklo[lo] + rankhi[poplo[lo],hi]

@njit
def next_comb(x):
    c=x & -x
    r=x+c
    return (((r ^ x) >> 2)//c) | r

@njit
def build_orbits(rep_of_state, shift_code, rep_masks, rep_lens, rep_valid,
                 trans_lo, trans_hi, poplo, ranklo, rankhi):
    x=np.int64((1<<14)-1)
    r=0
    nr=0
    while r < rep_of_state.size:
        best=np.int64(x)
        best_c=0
        # find lexicographically/numerically smallest translated bit pattern
        for c in range(1,9):
            t=np.int64(int(trans_lo[c,x & 8191]) | int(trans_hi[c,x >> 13]))
            if t < best:
                best=t; best_c=c
        if best_c==0:
            rid=nr
            rep_masks[rid]=x
            # Orbit size and momentum-stabilizer compatibility.
            nstab=0
            valid=True
            for c in range(9):
                t=int(trans_lo[c,x & 8191]) | int(trans_hi[c,x >> 13])
                if t==x:
                    nstab+=1
                    dx=c//3; dy=c%3
                    if (dx+dy)%3 != 0:
                        valid=False
            rep_lens[rid]=9//nstab
            rep_valid[rid]=1 if valid else 0
            rep_of_state[r]=rid
            shift_code[r]=0
            nr+=1
        else:
            rb=rank_mask(best,poplo,ranklo,rankhi)
            rid=rep_of_state[rb]
            if rid < 0:
                return -nr-1
            rep_of_state[r]=rid
            dx=best_c//3; dy=best_c%3
            # x = T_{-best_c}(best); only k=(1/3,1/3) phase is needed.
            shift_code[r]=((-dx-dy)%3)
        r+=1
        if r < rep_of_state.size:
            x=next_comb(x)
    return nr

log('Allocating basis-orbit maps (about 100 MB)...')
rep_of_state=np.full(DIM,-1,dtype=np.int32)
shift_code=np.zeros(DIM,dtype=np.uint8)
rep_masks=np.empty(DIM,dtype=np.uint32)
rep_lens=np.empty(DIM,dtype=np.uint8)
rep_valid=np.empty(DIM,dtype=np.uint8)
log('Building translation orbits...')
t0=time.time()
nrep=build_orbits(rep_of_state,shift_code,rep_masks,rep_lens,rep_valid,
                  trans_lo,trans_hi,poplo,ranklo,rankhi)
if nrep<0: raise RuntimeError(f'bad orbit lookup at {-nrep-1}')
log(f'Orbit construction complete: {nrep:,} raw orbits, {int(np.sum(rep_valid[:nrep])):,} allowed K=(1/3,1/3) orbits, elapsed {time.time()-t0:.1f}s')

# Compact the allowed momentum basis.
raw_to_basis=np.full(nrep,-1,dtype=np.int32)
nb=0
for r in range(nrep):
    if rep_valid[r]:
        raw_to_basis[r]=nb; nb+=1
basis_masks=np.empty(nb,dtype=np.uint32)
basis_lens=np.empty(nb,dtype=np.uint8)
basis_raw=np.empty(nb,dtype=np.int32)
for r in range(nrep):
    i=raw_to_basis[r]
    if i>=0:
        basis_masks[i]=rep_masks[r]
        basis_lens[i]=rep_lens[r]
        basis_raw[i]=r
log(f'K-sector dimension = {nb:,} (full Sz=1/2 dimension {DIM:,}; reduction x{DIM/nb:.3f})')

# Drop arrays no longer needed to lower memory footprint (retain state->orbit maps).
del rep_masks, rep_lens, rep_valid
# Diagonal matrix elements per representative.
@njit
def calc_diagonal(masks,bonds):
    out=np.zeros(masks.size,dtype=np.float64)
    for i in range(masks.size):
        m=int(masks[i]); e=0.0
        for b in range(bonds.shape[0]):
            u=int(bonds[b,0]); v=int(bonds[b,1]); j=bonds[b,2]
            su=1 if ((m>>u)&1) else -1
            sv=1 if ((m>>v)&1) else -1
            e += j*0.25*su*sv
        out[i]=e
    return out

diag=calc_diagonal(basis_masks,bonds)

phase=np.array([1.0+0.0j, np.exp(2j*np.pi/3), np.exp(4j*np.pi/3)],dtype=np.complex128)

@njit(parallel=True)
def matvec_numba(x,y,basis_masks,basis_lens,diag,bonds,rep_of_state,shift_code,raw_to_basis,
                 poplo,ranklo,rankhi,phase):
    # Gather form: no concurrent writes, so representative rows are parallelizable.
    # For a transition i -> j the projected matrix element is
    # (J/2)*sqrt(L_i/L_j)*exp(+i k.d); row action uses its conjugate.
    n=basis_masks.size
    for i in prange(n):
        m=int(basis_masks[i]); li=float(basis_lens[i])
        yi=diag[i]*x[i]
        for b in range(bonds.shape[0]):
            u=int(bonds[b,0]); v=int(bonds[b,1]); j=bonds[b,2]
            bu=(m>>u)&1; bv=(m>>v)&1
            if bu != bv:
                mt=m ^ (1<<u) ^ (1<<v)
                rr=rank_mask(mt,poplo,ranklo,rankhi)
                raw=rep_of_state[rr]
                it=raw_to_basis[raw]
                if it>=0:
                    lt=float(basis_lens[it])
                    ph=phase[shift_code[rr]]
                    yi += (0.5*j*math.sqrt(li/lt))*np.conjugate(ph)*x[it]
        y[i]=yi

class KOperator(LinearOperator):
    def __init__(self):
        self.tmp=np.zeros(nb,dtype=np.complex128)
        super().__init__(dtype=np.dtype(np.complex128),shape=(nb,nb))
    def _matvec(self,x):
        matvec_numba(x,self.tmp,basis_masks,basis_lens,diag,bonds,rep_of_state,shift_code,raw_to_basis,
                     poplo,ranklo,rankhi,phase)
        return self.tmp.copy()

op=KOperator()
# Warm-up compilation and matvec timing.
log('Compiling and timing one matrix-vector product...')
rng=np.random.default_rng(1234)
x=rng.normal(size=nb)+1j*rng.normal(size=nb)
x=x.astype(np.complex128); x/=np.linalg.norm(x)
y=np.empty_like(x)
t0=time.time()
matvec_numba(x,y,basis_masks,basis_lens,diag,bonds,rep_of_state,shift_code,raw_to_basis,
            poplo,ranklo,rankhi,phase)
matvec_s=time.time()-t0
log(f'One K-sector H@v: {matvec_s:.2f}s; ||Hv||={np.linalg.norm(y):.8f}')

# Hermiticity check <u|Hv> = <Hu|v>, using independent vectors.
u=rng.normal(size=nb)+1j*rng.normal(size=nb); u=u.astype(np.complex128); u/=np.linalg.norm(u)
t0=time.time()
matvec_numba(x,y,basis_masks,basis_lens,diag,bonds,rep_of_state,shift_code,raw_to_basis,
            poplo,ranklo,rankhi,phase)
hv=y.copy()
matvec_numba(u,y,basis_masks,basis_lens,diag,bonds,rep_of_state,shift_code,raw_to_basis,
            poplo,ranklo,rankhi,phase)
hu=y.copy()
lhs=np.vdot(u,hv); rhs=np.vdot(hu,x)
herm=abs(lhs-rhs)/max(1.0,abs(lhs),abs(rhs))
log(f'Hermiticity check relative residual = {herm:.3e} (2 matvecs in {time.time()-t0:.1f}s)')

# scipy.eigsh dispatches complex-Hermitian inputs through the general complex solver.
# To use the real-symmetric Lanczos path explicitly, represent H=A+iB as
# H_R=[[A,-B],[B,A]], acting on [Re(z), Im(z)]. This doubles vector length
# but keeps each Krylov vector at the same byte size as one complex K-sector vector.
class RealKOperator(LinearOperator):
    def __init__(self):
        self.tmp_complex=np.zeros(nb,dtype=np.complex128)
        super().__init__(dtype=np.dtype(np.float64),shape=(2*nb,2*nb))
    def _matvec(self,xr):
        z=xr[:nb] + 1j*xr[nb:]
        matvec_numba(z,self.tmp_complex,basis_masks,basis_lens,diag,bonds,rep_of_state,shift_code,raw_to_basis,
                     poplo,ranklo,rankhi,phase)
        return np.concatenate((self.tmp_complex.real,self.tmp_complex.imag))

op_real=RealKOperator()
log(f'Real-symmetric Lanczos representation dimension = {2*nb:,}; starting lowest-eigenvalue solve...')
if os.environ.get('TEST_ONLY') == '1':
    xr=rng.normal(size=2*nb); wr=rng.normal(size=2*nb)
    t0=time.time(); hyr=op_real.matvec(xr); hwr=op_real.matvec(wr)
    symerr=abs(np.dot(xr,hwr)-np.dot(hyr,wr))/max(1.0,abs(np.dot(xr,hwr)),abs(np.dot(hyr,wr)))
    log(f'Real symmetric embedding test: matvecs={time.time()-t0:.2f}s; symmetry residual={symerr:.3e}')
    sys.exit(0)
t0=time.time()
try:
    vals,vecs=eigsh(op_real,k=1,which='SA',ncv=16,maxiter=120,tol=1e-6)
except ArpackNoConvergence as exc:
    log(f'Lanczos did not converge; partial eigenvalues: {exc.eigenvalues}')
    raise
order=np.argsort(vals); vals=vals[order]; vecs=vecs[:,order]
log(f'Lanczos completed in {time.time()-t0:.1f}s; eigenvalue={vals[0]:.12f}')
# Convert the real doubled Ritz vector back to the complex K-sector and check residual.
z=vecs[:nb,0] + 1j*vecs[nb:,0]
matvec_numba(z,y,basis_masks,basis_lens,diag,bonds,rep_of_state,shift_code,raw_to_basis,
            poplo,ranklo,rankhi,phase)
res=np.linalg.norm(y-vals[0]*z)
log(f'Lowest-state residual ||Hv-Ev||={res:.3e}')
print('RESULT_JSON',{'full_dim':DIM,'raw_orbits':int(nrep),'K_dim':int(nb),'matvec_seconds':float(matvec_s),'hermiticity_rel':float(herm),'eigenvalue_K':float(vals[0]),'residual':float(res)})
