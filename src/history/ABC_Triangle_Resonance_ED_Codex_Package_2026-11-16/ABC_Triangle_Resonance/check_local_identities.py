#!/usr/bin/env python3
"""Independent 8x8 spin-1/2 checks for local ABC-triangle exact identities."""
import itertools
import numpy as np
sig=(np.array([[0,1],[1,0]],complex),np.array([[0,-1j],[1j,0]],complex),np.array([[1,0],[0,-1]],complex))
I8=np.eye(8,dtype=complex)
S={}
for site in range(3):
 for axis in range(3):
  s=[np.eye(2,dtype=complex) for i in range(3)];s[site]=sig[axis]/2
  S[site,axis]=np.kron(np.kron(s[0],s[1]),s[2])
B=lambda i,j:sum(S[i,a]@S[j,a] for a in range(3))
B_AB,B_BC,B_CA=B(0,1),B(1,2),B(2,0)
S2=9*I8/4+2*(B_AB+B_BC+B_CA)
P12=(15*I8/4-S2)/3
omega=np.exp(2j*np.pi/3)
D=B_AB+omega*B_BC+omega**2*B_CA
chi=sum((+1 if (a,b,c) in ((0,1,2),(1,2,0),(2,0,1)) else -1)*S[0,a]@S[1,b]@S[2,c] for a,b,c in itertools.permutations(range(3)))
checks={
 'P12 projector':P12@P12-P12,
 'chi^2 = (3/16)P12':chi@chi-3*P12/16,
 'sym(Ddag D)=(9/8)P12':(D.conj().T@D+D@D.conj().T)/2-9*P12/8,
 '[Ddag,D]=3 sqrt(3)chi':D.conj().T@D-D@D.conj().T-3*np.sqrt(3)*chi,
 'chi Hermitian':chi-chi.conj().T,
}
for label,mat in checks.items():
 error=np.linalg.norm(mat)
 print(f'{label}: error={error:.3e}')
 assert error<2e-14
print('All exact local spin-1/2 identities verified.')
