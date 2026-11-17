#!/usr/bin/env python3
"""Independent three-spin S=1/2 ABC triangle: operator algebra and exact dynamics.

Basis: bit a+2b+4c, 1=spin up.  D=B_AB + w B_BC + w² B_CA;
chi=S_A dot (S_B cross S_C). No external 27-spin eigensolver is needed.
"""
import numpy as np, csv, json, pathlib, itertools
ROOT=pathlib.Path(__file__).resolve().parents[1]
I=np.eye(8,dtype=np.complex128); w=np.exp(2j*np.pi/3)
sigma_x=np.array([[0,1],[1,0]],complex)
sigma_y=np.array([[0,-1j],[1j,0]],complex)
sigma_z=np.diag([1.,-1.]).astype(complex)
S=np.zeros((3,3,8,8),dtype=complex)
for a in range(3):
 for axis in range(3):
  for b in range(8):
   sb=((b>>a)&1)
   if axis==2:S[a,axis,b,b]=sb-.5
   elif axis==0:S[a,axis,b^(1<<a),b]=.5
   else:S[a,axis,b^(1<<a),b]=.5j if sb else -.5j
B_AB=sum(S[0,k]@S[1,k] for k in range(3))
B_BC=sum(S[1,k]@S[2,k] for k in range(3))
B_CA=sum(S[2,k]@S[0,k] for k in range(3))
D=B_AB+w*B_BC+w*w*B_CA
chi=sum((1 if (i,j,k) in [(0,1,2),(1,2,0),(2,0,1)] else -1)*S[0,i]@S[1,j]@S[2,k] for i,j,k in itertools.permutations(range(3)))
P=.5*I-2*(B_AB+B_BC+B_CA)/3
# AB singlet / AB triplet coupled to C; M=+1/2
v0=np.zeros(8,complex);v0[[5,6]]=[1,-1];v0/=np.sqrt(2)
v1=np.zeros(8,complex);v1[3]=np.sqrt(2/3);v1[[5,6]]=-1/np.sqrt(6)
V=np.stack([v0,v1],axis=1)
assert np.linalg.norm(V.conj().T@V-np.eye(2))<1e-12

def err(A):return float(np.linalg.norm(A))
checks={
 'projector_idempotence':err(P@P-P),
 'chirality_squared':err(chi@chi-3*P/16),
 'dimer_symmetric_norm':err((D.conj().T@D+D@D.conj().T)/2-9*P/8),
 'dimer_chirality_commutator':err(D.conj().T@D-D@D.conj().T-3*np.sqrt(3)*chi),
 'nilpotent_D_square':err(D@D),
 'chirality_dimer_ladder':err(chi@D-D@chi+np.sqrt(3)*D/2),
 'D_flavor_matrix':err(V.conj().T@D@V-(-.75*sigma_z+.75j*sigma_x)),
 'chi_flavor_matrix':err(V.conj().T@chi@V-(np.sqrt(3)/4)*sigma_y),
}
# scalar spin selection: S_total^2 commutes with D and chi
Stot=np.sum(S,axis=0);S2=sum(x@x for x in Stot)
checks['scalar_D_spin_squared']=err(D@S2-S2@D)
checks['scalar_chi_spin_squared']=err(chi@S2-S2@chi)
assert max(checks.values())<3e-14,checks

# Spin-1/2 doublet ground AB-singlet for j<1; excited AB-triplet doublet.
# For j=1 the doublet flavors are exactly degenerate in isolated triangle.
rows=[];poles=[]
for j in [0.2,0.5,0.6,0.7,0.8,0.9,0.95,0.975,1.0]:
 H=B_AB+j*(B_BC+B_CA)
 E0=float(np.vdot(v0,H@v0).real)
 E1=float(np.vdot(v1,H@v1).real)
 delta=E1-E0
 assert abs(delta-(1-j))<1e-12
 assert err(H@v0-E0*v0)<2e-14 and err(H@v1-E1*v1)<2e-14
 for name,O in [('D',D),('chi',chi),('D_re',(D+D.conj().T)/2),('D_im',(D-D.conj().T)/(2j))]:
  groundamp=np.vdot(v0,O@v0)
  a=np.vdot(v1,O@v0)
  inel=abs(a)**2
  m0=float(np.vdot(v0,O.conj().T@O@v0).real)
  m1=float(np.vdot(O@v0,(H-E0*I)@(O@v0)).real)
  m2=float(np.linalg.norm((H-E0*I)@(O@v0))**2)
  assert abs(m0-(abs(groundamp)**2+inel))<1e-12
  assert abs(m1-inel*delta)<1e-12 and abs(m2-inel*delta**2)<1e-12
  rows.append(dict(j=j,operator=name,ground_energy=E0,excitation_energy=delta,elastic_weight=abs(groundamp)**2,transition_weight=inel,M0=m0,M1=m1,M2=m2,transition_amplitude_real=a.real,transition_amplitude_imag=a.imag))
 # Cross-spectral residue <g|D^dag|e><e|chi|g> and exact matrix moments
 Ad=np.vdot(v1,D@v0);Ac=np.vdot(v1,chi@v0)
 residue=Ad.conjugate()*Ac
 poles.append(dict(j=j,energy=delta,cross_D_chi_real=residue.real,cross_D_chi_imag=residue.imag))
assert abs(poles[0]['cross_D_chi_real']-3*np.sqrt(3)/16)<1e-12
with (ROOT/'results'/'isolated_triangle_spectra.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
with (ROOT/'results'/'isolated_triangle_cross.csv').open('w',newline='') as f:
 writer=csv.DictWriter(f,fieldnames=list(poles[0]));writer.writeheader();writer.writerows(poles)
with (ROOT/'results'/'local_algebra_validation.json').open('w') as f:json.dump({'max_error':max(checks.values()),'checks':checks,'flavor_basis_order':'AB singlet, AB triplet doublet; M=+1/2','D_flavor_matrix':[[{'re':a.real,'im':a.imag} for a in row] for row in V.conj().T@D@V]},f,indent=2)
print('PASS independent 8x8 algebra and isolated spectral identities, max error',max(checks.values()))
print('exact isolated triangle delta/J = 1-j, W_inel(chi)=3/16, W_inel(D)=9/16')
