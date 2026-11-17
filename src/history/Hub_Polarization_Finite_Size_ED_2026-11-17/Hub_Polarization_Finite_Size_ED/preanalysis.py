import json
from geometry import ROOT
ref=json.loads((ROOT/'reference27_results.json').read_text());lines=['# Completed Phase A/B diagnostics before the conditional 36-spin pilot','All eighteen-spin ground states use exhaustive cell momentum coverage in Sz=0, which contains every integer-S multiplet. The 27-spin data use verified checkpoint-derived spin-summed scalar correlations.','| r | shape | S | p_H | R²/N_R² | m_N² |','|---|---|---|---|---|---|']
for r in [.2,.4,.5,.6,.7]:
 for shape in ['K18','C18']:
  p=json.loads((ROOT/'checkpoints'/shape/f'r{r:.9f}.json').read_text());o=p['observables'];lines.append(f"|{r}|{shape}|{o['S']}|{o['p_H']:.8f}|{o['m_R_uniform2']:.8f}|{o['m_N2']:.8f}|")
 for p in ref:
  if abs(p['r']-r)<1e-8:
   o=p['observables'];lines.append(f"|{r}|N27|{o['S']}|{o['p_H']:.8f}|{o['m_R_uniform2']:.8f}|{o['m_N2']:.8f}|")
lines+=['The rise of hub pair polarization and subsequent increase of rim uniform correlations while total spin decreases recur in both 18-spin geometries and in 27 spins. Their crossing locations, momenta and spin plateaus differ, and the 18-spin shapes disagree at r=0.7.','All valid 18-spin k0-compatible superlattices are point-group equivalent and have shortest winding length 3 in site units. The 27/36 candidates have shortest lengths 3sqrt(3) and 6, respectively. The more isotropic eighteen-spin control lacks the weak cell momentum k0. These are substantial geometric and parity controls; no common thermodynamic transition is inferred.','The next phase is a memory-feasibility assessment and explicitly high-Sz 36-site matvec pilot only unless low-spin memory estimates permit a valid ground-state solve. It is not an attempt to assign a 36-spin ground state from a restricted sector.']
(ROOT/'PHASE_B_PREANALYSIS.md').write_text('\n\n'.join(lines));print('Phase B comparison completed before 36-spin conditional pilot.')
