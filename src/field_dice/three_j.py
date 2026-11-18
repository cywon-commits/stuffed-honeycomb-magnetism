"""Three-coupling generalization: J_AB (rim A-B), J_BC (hub C - rim B), J_CA (hub C - rim A).
Ground-state total spin S at h = 0 near the honeycomb vertex, on and off the mirror J_BC = J_CA."""
import numpy as np, scipy.sparse.linalg as sla
from ed_sz import geometry, Sector
g = geometry.catalog()['K18']; N = g['N']
split = {'AB': [(a, b) for a, b, t in g['bonds'] if t == 0],
         'CA': [(a, b) for a, b, t in g['bonds'] if t == 1 and b % 3 == 0],
         'CB': [(a, b) for a, b, t in g['bonds'] if t == 1 and b % 3 == 1]}
assert all(a % 3 == 2 for a, b in split['CA'] + split['CB'])
secs = {}
for n in range(N // 2, N + 1):
    s = Sector.__new__(Sector); s.N = N; s.nup = n
    from ed_sz import states_with; s.s = states_with(N, n)
    s.H = {k: s._build(v) for k, v in split.items()}; secs[n] = s
def S_ground(jab, jbc, jca):
    E = []
    for n, s in secs.items():
        Hm = jab * s.H['AB'] + jbc * s.H['CB'] + jca * s.H['CA']
        E.append(np.linalg.eigvalsh(Hm.toarray())[0] if Hm.shape[0] < 400 else sla.eigsh(Hm, k=1, which='SA', tol=1e-11)[0][0])
    E = np.array(E); return max(m for m in range(len(E)) if E[m] < E.min() + 1e-8)
print('edge J_CA = 0 (bipartite, Lieb-Mattis predicts S = N/6 = 3):',
      [S_ground(1, x, 0) for x in [0.05, 0.2, 0.5, 1.0, 2.0]])
for j in [0.1, 0.2, 0.3]:
    row = [(d, S_ground(1, j + d / 2, j - d / 2)) for d in [0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]]
    print(f'J_AB=1, mean hub coupling j={j}: (delta, S) =', row)

def E0(jab, jbc, jca):
    return min((jab * s.H['AB'] + jbc * s.H['CB'] + jca * s.H['CA']).toarray().min() if False else
               (np.linalg.eigvalsh((jab*s.H['AB']+jbc*s.H['CB']+jca*s.H['CA']).toarray())[0] if s.H['AB'].shape[0] < 400 else
                sla.eigsh(jab*s.H['AB']+jbc*s.H['CB']+jca*s.H['CA'], k=1, which='SA', tol=1e-11)[0][0]) for s in secs.values())
eh = E0(1, 0, 0)
print('honeycomb vertex E_hex =', eh)
for name, w in [('dice (0,1,1)', (0, 1, 1)), ('triangle (1,1,1)', (1, 1, 1)), ('r=0.5 (1,.5,.5)', (1, .5, .5)), ('r=2 (.5,1,1)', (.5, 1, 1))]:
    e = E0(*w); print(f'{name}: E = {e:.6f}, lower bound sum w_k E_hex = {sum(w)*eh:.6f}, f = E/bound = {e/(sum(w)*eh):.4f}')
