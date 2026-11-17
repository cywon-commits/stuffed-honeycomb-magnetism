# Stuffed honeycomb 양자자성 — 연구 현황 (2026-11-17)

## 모형
스핀 1/2, 3-site A/B/H 단위격자. Rim–rim 교환 J, hub–rim 교환 J', r=J'/J. Hub–hub 직접 교환은 없다.

## 검증된 결과 요약
- N=27: r≈0.2에서 전체 S=1/2, rim Néel 상관이 크고 hub는 저에너지 스펙트럼을 가짐.
- N=27: 총스핀의 1/2 → 7/2 → 5/2 및 더 큰 r에서 3/2 → 1/2 변화를 확인했으나, 각 crossing을 열역학적 상경계로 볼 수는 없음.
- N=18 K18/C18: hub polarization의 급증과 rim-spin compensation이라는 정성적 메커니즘이 재현되지만, crossing 위치와 ground-state momentum에는 큰 크기/형상 의존성이 있음.
- N=27 약결합: 검증된 동적 pole과 rim–hub 교차 spectral residue가 존재. 유한계의 이산 pole 폭은 magnon damping이 아님.
- r→1: ABC 삼각형의 local S_triangle=1/2 weight 증가, dimer/chirality 스칼라 스펙트럼 연구. 이것만으로 별도의 trimer resonance 상은 입증되지 않음.

## 현재의 핵심 목표
18/27-spin에 이어 36-spin 저스핀 sector에서 정확한 대각화를 수행하여 hub polarization과 rim compensation이 큰 계에서도 유지되는지 검증. 제공된 36-spin C++ 파일은 고자화 sector 파일럿이고, 완성된 저스핀 solver가 아니다.

## 시작
- `docs/README_START_HERE.md`
- `docs/HPC_36SPIN_EXECUTION_REQUEST_EN.md`
- `src/workspaces/finite_size_18_27_36pilot/` 및 `src/history/`
- `external_data/README.md`: 큰 체크포인트/원본 아카이브 링크
