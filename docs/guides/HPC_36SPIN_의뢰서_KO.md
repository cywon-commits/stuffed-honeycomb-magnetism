# 36-spin 대규모 ED 이식·실행 의뢰서

**목적:** 기존 18-spin(두 형상)·27-spin에서 관찰한 허브 집단적 편극과 림의 보상 스핀 증가가 36-spin 토러스에서도 재현되는지 확인한다. 전체 결합비 영역의 무차별 스캔보다는 r=J'/J=0.20, 0.30, 0.45, 0.60, 0.68 대표점을 우선한다.

**중요:** 포함한 `kernel36pilot.cpp`/`pilot36.py`는 고자화 sector 파일럿이며, 낮은 총스핀의 36-spin 바닥상태 계산기가 **아니다**. 대규모 Hilbert 공간의 분산 저장/압축 인덱싱, 운동량 및 점군 대칭 분해, 메모리 제약 아래 Krylov eigensolver를 새로 구현해 실제 저스핀 sector에서 검증해야 한다. `guides/HPC_36SPIN_EXECUTION_REQUEST_EN.md`의 상세 요구사항에 따라 계산을 실제로 수행한다.

**기존 재현:** `workspaces/finite_size_18_27_36pilot/README.md`의 `python run_all.py --tests`를 먼저 실행하고 18-spin benchmark와 27-spin 결과를 비교한다. 36-spin 토러스의 셀 좌표는 T1=(2,2), T2=(-2,4) (12개 단위격자, 36개 스핀), 기존 셀 reciprocal momentum을 그대로 사용한다.

**필수 결과:** 각 대표 결합값에서 전체 바닥상태 E, 총스핀 S, 운동량, 결합 상태에서의 <S_H²>, <S_R²>, <S_H·S_R>, rim Néel correlation, hub 구조인자(Γ 및 비영 운동량)를 계산하고 18·27·36 데이터를 동일 정규화로 비교한다. 정확한 총스핀 보존 법칙 및 residual, sector coverage 확인 후에만 물리적으로 해석한다. 재현 코드/환경/수치 데이터/고유벡터/그림/보고서를 제출한다.

**자료 구성:** `workspaces/`는 최신 18-spin·27-spin의 실행 가능한 작업공간, `original_archives/`는 과거 Codex 코드·체크포인트가 들어있는 원본 ZIP 모두, `source_code_all_stages/`는 코드·검증·정적 데이터 모음, `manifests/`는 원본 무결성 확인용 SHA256이다. 파일명과 의도는 `README_START_HERE.md`를 따른다.
