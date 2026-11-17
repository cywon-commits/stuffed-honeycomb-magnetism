# Quantum Magnetism — 연구 현황 및 실행 계획

> **기준일:** 2026-11-15 (KST)  
> **연구 대상:** Spin-1/2 stuffed honeycomb antiferromagnet (honeycomb rim + hexagon-center hub)  
> **목적:** 긴 채팅을 다시 읽지 않고도 물리적 질문, 현재 증거 수준, 수치 계산의 진행 단계, 다음 행동을 판단하는 단일 기준 문서.  
> **관리 원칙:** 결과는 반드시 **확인된 수치 / 이전 대화에서 보고된 수치 / 구현 검증 / 가설 / 미실행**으로 구분한다.

## 0. 현재 상황 한눈에 보기

| 항목 | 상태 (2026-11-15) | 다음 조치 |
|---|---|---|
| 핵심 연구 질문 | 정의됨 | 선행연구 대비 독창성 검증 지속 |
| Hamiltonian 및 3×3 결합 그래프 | 코드 초안 있음 | 독립적인 그래프 검증 및 기준계 비교 |
| 18-spin ED, 3×2 torus | 이전 대화에서 수치 결과 보고됨 | 원본 실행 로그·코드·재현성 확보 |
| 27-spin, 3×3, translation reduction | 코드 및 성공 보고 있음 | Codex에서 깨끗한 환경으로 독립 재실행 |
| 27-spin matrix-free H·Lanczos | H 작용 및 에르미트성 통과 보고, eigenvalue 미수렴 | 메모리/시간 벤치마크와 잔차 기준 도입 |
| Codex 작업 | 사용자가 브라우저 Codex Cloud에 코드/지시서 전달 | 환경 사양·진행 상태 확인, 수정된 운동량 지침 전달 |
| 동적 hub susceptibility χ(q,ω) | **미계산** | small-cluster validation → Lanczos spectral method |
| magnon self-energy ReΠ, ImΠ | **미계산** | irreducible hub response/vertex/Goldstone consistency 설계 |
| thermodynamic phase diagram | **미확립** | 여러 size/geometry와 observables로 비교 |

### 당장 수정해야 할 두 가지

1. **3×3 torus에는 정확한 M 점이 없다.** 허용 운동량이 reciprocal primitive 좌표로 `(n1/3,n2/3)`이므로 `M=(1/2,0),(0,1/2),(1/2,1/2)`은 불가능하다. 직전 Codex 지시서의 ‘3개의 M sector 비교’는 변경해야 한다. M을 포함하는 클러스터에는 짝수 주기 방향이 필요하다.
2. **K 좌표를 재정의해서 확인해야 한다.** 표준 60° triangular primitive reciprocal 좌표에서 K의 대표점은 `(1/3,2/3)` 또는 `(2/3,1/3)`이다(기저 방향·원점 선택에 따라 동등한 표기 가능). **첨부된 Python 코드에서 `K=(1/3,1/3)`이라고 이름 붙인 sector가 진짜 물리적 K인지 자동으로 가정하지 말 것.** 이 코드의 translation convention 및 triangular nearest-neighbor displacement를 기준으로 구조인자/Brillouin-zone 코너를 독립적으로 검증한다. `1/3,1/3`이 실제 K가 아니라면 sector 이름과 기존 결론을 수정하고 진짜 K를 계산한다.

## 1. 연구의 핵심 물리적 질문

**Ordered rim sublattice가 생성하는 양자요동에 의해 동역학이 유도되는 quantum-disordered hub subsystem이, rim Néel long-range order를 유지하면서도 rim Goldstone magnon을 강하게 renormalize/damp할 수 있는가?**

추가로 검토할 질문:
- 약한 coupling 영역에서 hub spin은 어떤 저에너지 상관을 형성하는가? Stripe, period-three, quantum-disordered response를 어떻게 구별하는가?
- J' 증가로 rim 질서가 붕괴할 때, Goldstone damping이 원인인지, 동일한 상전이의 병행 징후인지 식별할 수 있는가?
- 실수부(Re self-energy, dispersion/stiffness 변화)와 허수부(Im self-energy, linewidth/decay)를 동시에 계산하고 Kramers–Kronig 및 spin-rotation Ward identities와 일치시키려면 어떤 근사가 필요한가?

**인과관계 주의:** 감쇠의 증가만으로 Néel order 붕괴의 원인이라고 결론 낼 수 없다. 정적 질서 지표·stiffness와 동적 damping을 독립적으로 추적해야 한다.

## 2. 모형 및 기하

\[
H=J\sum_{\langle ij\rangle_{\mathrm{rim}}}\mathbf S_i\cdot\mathbf S_j + J'\sum_{\langle iR\rangle}\mathbf S_i\cdot\mathbf T_R,\quad J,J'>0.
\]

- Rim은 honeycomb A/B 스핀, hub는 각 hexagon 중심에 위치하는 스핀.
- 각 hub는 6개의 rim spin(3 A + 3 B), 각 rim은 3개 hub에 결합; unit cell 3 spins, rim:hub = 2:1.
- Rim bonds (셀 R=(x,y)): `A_R-B_R`, `A_R-B_(R-a1)`, `A_R-B_(R-a2)`.
- **정정된 실제 hexagon six-vertex 목록**: `A_R`, `B_R`, `A_(R+a1)`, `B_(R+a1-a2)`, `A_(R+a1-a2)`, `B_(R-a2)`.
- `J'=J`: 모든 bond가 동등한 isotropic triangular-lattice AF 그래프에 대응해야 한다. graph isomorphism, coordination number 6, nearest-neighbor displacement를 각각 검증할 것.
- `J'≪J`: honeycomb rim Néel order와 허브의 유도된 저에너지 모델이 출발점.
- `J'~J`: weak-coupling effective theory extrapolation은 보장되지 않는다.
- `J'≫J`: `J→0` 극한의 bipartite imbalance 물리와 finite-J frustration을 구별해야 한다.

## 3. 이론적 접근: weak-coupling 유효 hub 모델 및 Goldstone feedback

유도된 hub interaction energy scale:
\[
K_{\mathrm{eff}}\sim J'^2/J.
\]

Weak-coupling XXZ 유효모형(표기상 `j`는 무차원):
\[
H_{\rm hub}^{\rm eff}=\frac{J'^2}{J}\sum_{RR'} \left[j_{RR'}^{xy}(T_R^xT_{R'}^x+T_R^yT_{R'}^y)+j_{RR'}^zT_R^zT_{R'}^z\right].
\]

과거 대화에서 인용된 1/S 근사 계수(논문·정의·단위 재확인 필요):

| shell | jxy | jz |
|---|---:|---:|
| NN | -0.004914 | -0.022853 |
| 2NN | +0.000867 | +0.047108 |
| 3NN | -0.000522 | -0.001019 |
| 4NN | +0.000059 | +0.004599 |
| 5NN | -0.000069 | +0.001555 |

개념적 Dyson equation:
\[
D^{-1}_{\lambda\lambda'}(\mathbf q,\omega)=D^{-1}_{0,\lambda\lambda'}(\mathbf q,\omega)-\Pi^R_{\lambda\lambda'}(\mathbf q,\omega),
\]
\[
\Pi^R_{\lambda\lambda'}(\mathbf q,\omega)\sim -J'^2 \sum_{\mu\nu}V_{\lambda\mu}(\mathbf q)\,\chi^{R,\mu\nu}_{\rm hub,irr}(\mathbf q,\omega)\,V^*_{\lambda'\nu}(\mathbf q).
\]

**주의:** 여기의 Π는 정확히 유도된 식이 아니라 schematic relation. 허브 susceptibility는 rim propagator에 대한 적절한 irreducibility를 가져야 double counting을 피한다. 실제 Goldstone 모드의 연속체/LSWT normalization을 정한 뒤 self-energy 부호, zero-frequency gaplessness, 	1/T	2 및 temperature limit, vertex form factor를 확인한다. ReΠ와 ImΠ는 아직 수치 계산되지 않았다.

## 4. 18-spin 3×2 PBC ED: 이전 대화에서 보고된 결과

**증거 등급:** 이전 대화에서 corrected hexagon list로 계산했다고 보고된 값. 이 문서 작성 중 재실행하지 않았으므로 '재현 확인 완료'와 혼동하지 말 것. 최초의 잘못된 hexagon list에 기반한 예비 데이터와 과거 미실행 12-spin ED 주장은 채택하지 않는다.

- 12 rim + 6 hub = 18 spins, `S^z=0` dimension `C(18,9)=48,620`.
- 18 rim bonds / 36 hub-rim bonds; `J=1`.
- Allowed q in hub reciprocal coordinates: `(p/3,r/2)`.
- `Γ=(0,0)`, `M=(0,1/2)` (3개의 M 중 하나만), `Q3=(1/3,0)` (**M 아님**), `Q4=(1/3,1/2)`.

과거에 보고된 static hub transverse Kubo susceptibility:

| J'/J | χΓ | χQ3 | χM | χQ4 | 보고된 최댓값 |
|---:|---:|---:|---:|---:|---|
| 0.20 | 0.00452 | 180.499 | 20.999 | 26.412 | Q3 |
| 0.30 | 0.01446 | 75.214 | 8.844 | 11.002 | Q3 |
| 0.40 | ∞ | 1.269 | 1.195 | 1.120 | Γ, spin multiplet degeneracy |
| 0.50 | ∞ | 1.298 | 1.155 | 1.039 | Γ, spin multiplet degeneracy |
| 0.80 | 1.628 | 0.472 | 0.220 | 0.221 | Γ |
| 1.00 | 1.377 | 0.331 | 0.176 | 0.174 | Γ |
| 1.50 | 0.916 | 0.230 | 0.118 | 0.115 | Γ |
| 2.00 | ∞ | 0.063 | 0.048 | 0.044 | Γ, spin multiplet degeneracy |

보고된 sector gap `E0(Sz=1)-E0(Sz=0)`:

| J'/J | Δ1/J |
|---:|---:|
| 0.20 | 0.005045 |
| 0.30 | 0.011276 |
| 0.40 | ~0 |
| 0.50 | ~0 |
| 0.80 | 0.281916 |
| 1.00 | 0.692153 |
| 1.50 | 0.011588 |
| 2.00 | ~0 |

**해석의 한계:** Γ susceptibility가 무한대인 경우는 ground-state spin multiplet degeneracy의 결과. Q3 maximum은 stripe-M order 증거가 아니며, 작은 3×2 격자는 isotropic triangular AF의 120° K order를 충분히 포함하지 못한다. 위 표만으로 thermodynamic magnetic phase transition 또는 ordered moment를 판정할 수 없다.

## 5. 27-spin 3×3 PBC: 코드 상태 및 수정 이슈

- 18 rim + 9 hub = 27 spin-1/2.
- `Sz=1/2`: `C(27,14)=20,058,300` states; `Sz=3/2`: `C(27,15)=17,383,860` states.
- translation group size 9: 비특이 orbit은 momentum sector별 약 2.23 million states.
- 코드 초안 파일명: `tri27_matrix_free_k_ed.py` (프로젝트의 Codex 작업에 전달됨).
- 이전 채팅에서 보고된 결과: `(1/3,1/3)` sector basis = 2,228,700; matrix-free `H@v` 약 2.2 s; complex Hermiticity relative error 약 `5.1e-17`; real representation check 약 `1.6e-13`; eigenvalue convergence **미완료**.
- 위 수치는 동일 Cloud에서 재측정한 수치가 아니며, 환경에 따라 runtime이 달라진다.
- **중요:** 코드에 `K=(1/3,1/3)`라고 표시돼 있어도 물리적 triangular K가 아닐 수 있다. 먼저 momentum labelling을 정정한다.
- **또 다른 중요 사항:** 3×3에서는 정확한 M 점 계산 불가. 3×2는 M의 일부를 포함하지만 K는 포함하지 않음. K와 M 둘 다 필요하면 6×6 등 적절한 주기, 혹은 비직사각 supercell 선택을 논의해야 하며 computational cost가 커진다.

### Codex에 보낼 정정 지시(최우선)

> Correction to the earlier task: A 3×3 torus does **not** contain exact triangular-lattice M points at reciprocal coordinates (1/2,0), (0,1/2), or (1/2,1/2). Do not label any 3×3 momentum sector as M. Also independently verify the actual triangular-lattice K momentum for the lattice/translation convention used by the code: in standard triangular reciprocal primitive coordinates, representative K points are (1/3,2/3) and (2/3,1/3), rather than automatically (1/3,1/3). Check the triangular graph's primitive vectors and the nearest-neighbor structure factor, correct momentum labels, and then benchmark the physically correct K and Γ sectors. Clearly record sector k=(n1/3,n2/3), dimension, eigenvalue, residual, CPU/RAM and runtime.

### 필수 numerical acceptance criteria

1. **Geometry:** bond counts 27 rim-rim / 54 hub-rim, unique undirected bonds, each site degree 6 when J=J'.
2. **Translation:** Hamiltonian commutes numerically with `Tx,Ty`; momentum phase convention is written down; basis orbit lengths/stabilizers are tested.
3. **Independent small-cluster test:** direct sparse full-space H vs translation-reduced matrix-free H agree to numerical precision for smaller N.
4. **Lanczos:** output Ritz value `E`, normalized residual `||Hψ-Eψ||`, iteration/restarts, tolerance, convergence flag; no unconverged result described as ground-state energy.
5. **Sanity:** real/Hermitian operator, energy variance, SU(2) multiplet and near-degeneracy assessed cautiously for odd N.
6. **Reproducibility:** exact git commit or code SHA, dependencies, seed, compiler flags, machine specs, wall-clock/RAM, command line.

## 6. 실행 계획과 의사결정 게이트

### Phase A — 검증 및 최저 에너지 (즉시)

**A1.** Codex Cloud가 CPU cores, available RAM, runtime restriction을 보고하고 `H@v` benchmark를 다시 측정.  
**A2.** 3×3에서 실제 triangular K와 Γ momentum 정확히 규정, erroneous M request 제거.  
**A3.** J=J'=1에서 27-spin의 K, Γ lowest eigenvalue 및 잔차를 얻고 independent benchmark와 비교.  
**A4.** 모두 수렴하지 않으면 병목을 객관적으로 기록하고 C++/Numba memory layout, extra symmetries, different cluster, HPC 이전 등을 검토.

**Gate A (통과 기준):** geometry/momentum 검증 + 적어도 핵심 sector의 잔차가 설정 tolerance에 도달한 재현 가능한 eigenpair.

### Phase B — 정적 spin correlation (Gate A 통과 뒤)

**B1.** `J'/J`를 약한 결합부터 1까지 몇 점만 선택해 systematic run.  
**B2.** `⟨S_i·S_j⟩` rim-rim, `⟨T_R·T_R'⟩` hub-hub, hub-rim 구조인자 및 각 sublattice order indicators를 각각 추적.  
**B3.** `S_tot` 판정(가능하면 `⟨S_tot²⟩`) 및 finite-size gap을 별도 기록.  
**B4.** 가능한 서로 다른 cluster/symmetry sector로 교차검증하고 susceptibility maxima를 ordered phase로 과잉해석하지 않기.

**Gate B:** 경향이 3×2의 모멘텀 제한/유한크기 인공효과인지 독립적인 증거가 생길 것.

### Phase C — dynamics / self-energy (Gate B 이후 또는 병행하는 이론팀)

**C1.** Hub dynamic `χ_hub(q,ω)`의 Lehmann/Lanczos continued-fraction 접근을 작은 격자부터 검증.  
**C2.** Rim order/Goldstone modes 계산(LSWT/field theory), rim-hub vertex momentum form factor 재확인.  
**C3.** Double-counting을 피하는 irreducible susceptibility 및 Π 정의를 확정.  
**C4.** `ReΠ`, `ImΠ`를 같은 analytic continuation convention 아래 계산, Kramers–Kronig와 spectral positivity·Goldstone gaplessness를 검사.  
**C5.** `δω(q)`와 `Γ(q)`를 rim Néel order 감소와 비교하여 인과 주장 가능 범위를 판별.

**Gate C:** causality, symmetry, convergence가 모두 충족된 동적 결과와 해석.

## 7. 효과적인 연구 진행 관리 방식

**추천: ‘한 개의 현황 문서 + 버전관리되는 계산 저장소 + 실험 결과 레지스터’**

- 이 문서는 **STATUS.md**로 유지: 날짜, 현재 단계, 다음 행동, 결론의 신뢰등급만 간단히 갱신.
- Codex 코드 저장소: `README.md`, `src/`, `tests/`, `runs/`, `results/`, `notes/`로 관리. 원본 로그와 결과 테이블을 보존하고 수동 복사를 최소화.
- 실험마다 고유 ID 예시: `SH27_N3x3_Jp1p0_k1_2_seed0_20261115`; `config.yaml`, `stdout.log`, `metrics.json` (`E`, residual, sectors, runtime, RAM, git SHA 포함), 간단한 분석 메모를 한 폴더에 저장.
- **증거 등급:** `P`=물리적 예측/가설, `I`=구현·단위시험 통과, `R`=수치적 수렴 및 재현 확인, `F`=finite-size 해석, `T`=thermodynamic limit에서 교차 검증. 문장과 표의 중요한 결론에는 등급 명시.
- 새 결과가 나오면 질문은 세 가지: **① 무엇이 실제로 계산되었나? ② 어떤 대안 설명/유한크기 효과가 남았나? ③ 다음 결과가 어떤 가설을 반증할 수 있나?**
- 정적 계산과 동적 이론 개발을 병렬로 진행하되, 해석에 필요한 전제와 의존성을 기록한다.

### Codex 작업 완료 시 이 문서를 갱신하는 간단한 규칙

1. 이 문서 맨 위 ‘현재 상황’ 표와 기준일을 먼저 수정.
2. `결과 레지스터`에 **코드 버전, cluster, momentum, J'/J, 상태(성공/미수렴), residual, logfile** 추가.
3. 다음 단계의 Gate 통과 여부와 근거를 한두 문장으로 기록.
4. 이전 주장과 모순이 있으면 덮어쓰지 말고 ‘정정 및 철회’ 기록 남기기.

## 8. 결과 레지스터 (초기 상태)

| ID | 항목 | 현재 검증 등급 | 증거/문제 |
|---|---|---|---|
| RUN-18-STATIC | 18 spin 3×2 static χ 및 Δ1 | 이전 대화 보고, 미재실행 | 로그/코드 원본과 재현 필요 |
| DEV-27-REDUCE | 27 spin translation symmetry reduction | 이전 대화에서 I로 보고 | current Codex environment 독립 확인 필요 |
| DEV-27-MATVEC | 27 spin matrix-free H@v, Hermiticity | 이전 대화에서 I로 보고 | benchmark 및 small-system oracle 필요 |
| RUN-27-K | 올바른 K-sector lowest eigenvalue | 미실행/미수렴 | 코드의 기존 ‘K’ 좌표부터 점검 |
| RUN-27-GAMMA | Γ-sector lowest eigenvalue | 미확보 | 추가 계산 필요 |
| RUN-27-M | 3×3에서 exact M | **실행 불가** | PBC 허용 운동량 불일치; 더 적합한 torus 필요 |
| RUN-DYN-CHI | χ_hub(q,ω) | 미실행 | 도구/근사/계산 규모 결정 필요 |
| TH-SELFENERGY | ReΠ/ImΠ | 미실행 | vertex, irreducibility, Goldstone constraints 정립 필요 |

## 9. 선행연구 및 출발 자료

- U. F. P. Seifert and M. Vojta, “Theory of partial quantum disorder in the stuffed honeycomb Heisenberg antiferromagnet,” *Phys. Rev. B* **99**, 155156 (2019), DOI: `10.1103/PhysRevB.99.155156`, arXiv:1812.08168.
- M. G. Gonzalez *et al.*, “Correlated partial disorder in a weakly frustrated quantum antiferromagnet,” *Phys. Rev. Lett.* **122**, 017201 (2019). [세부 모델 및 관측값 비교 필요]
- 코드 초안: `tri27_matrix_free_k_ed.py` (이전 채팅에서 작성, Codex에 전달했다고 사용자 확인).

**향후 새 주장 전 확인:** 위 선행연구들이 이미 어떤 동적 응답/Goldstone self-energy 계산까지 수행했는지 별도 literature matrix로 비교해야 하며, 아직 novelty가 확립됐다고 볼 수 없다.

---

### 변경 이력

- **2026-11-15 v1.0:** 현재 연구 질문, 이전 18-spin 수치, 27-spin 개발 상태, Codex 작업, gate 기반 계획을 통합. 특히 **3×3에서는 exact M 불가**, **코드의 `(1/3,1/3)` K label 검증 필요**를 명시.
