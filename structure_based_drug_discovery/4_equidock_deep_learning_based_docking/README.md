# EQUIDOCK — 딥러닝 기반 단백질-단백질 도킹 실습

> **LAIDD — 단백질 구조 기반 약물탐색**  
> 강좌: 딥러닝을 이용한 단백질 도킹 (이유한, 카카오브레인)  
> 실습 환경: **Google Colab T4 GPU** (CPU 모드 강제 실행)  
> 실습 파일: [`equidock_practice.ipynb`](equidock_practice.ipynb)  
> 논문: [EquiDock: Independent SE(3)-Equivariant Models for End-to-End Rigid Protein Docking](https://arxiv.org/pdf/2111.07786)  
> 소스코드: https://github.com/octavian-ganea/equidock_public

---

## 실습 목표

1. EQUIDOCK 환경 구축 (Python 3.9 가상환경 + PyTorch 1.10.2 + DGL 0.7.2)
2. DIPS 테스트셋 100개 구조에 대한 단백질-단백질 도킹 추론
3. 논문 결과(C-RMSD, I-RMSD) 재현 검증
4. py3Dmol로 예측 결과 시각화 및 Ground Truth 비교

---

## EQUIDOCK 핵심 개념

### 기존 방법 vs EQUIDOCK

| 항목 | 기존 방법 (CLUSPRO 등) | EQUIDOCK |
|------|----------------------|---------|
| 방식 | 수십만 개 후보 샘플링 후 scoring | **단 한 번의 End-to-End 예측** |
| 소요 시간 | 수 시간 (CLUSPRO: 10,475초) | **3~5초 (GPU)** |
| 샘플링 | 필요 | ❌ 불필요 |
| Ranking | 필요 | ❌ 불필요 |
| Fine-tuning | 필요 | ❌ 불필요 |

### 3가지 핵심 아이디어

| Idea | 기법 | 해결하는 조건 |
|------|------|-------------|
| 1 | Residue-Residue Alignment + Kabsch Algorithm | 결합 전/후 구조 일치 가정 우회 |
| 2 | SE(3)-Equivariant GNN | 초기 위치/방향 무관한 일관된 예측 |
| 3 | Graph Matching Network | 수용체/리간드 역할 교환 불변성 |

### Loss 함수

| Loss | 역할 |
|------|------|
| `L_OT` (Optimal Transport) | 예측 keypoint → 실제 결합 포켓 위치 학습 |
| `L_MSE` (Coordinate Loss) | 예측 좌표 → Ground truth 좌표 오차 최소화 |
| `L_NI` (Intersection Loss) | 두 단백질 원자 간 물리적 충돌 방지 |

---

## 실습 환경

### 하드웨어

| 항목 | 값 |
|------|-----|
| GPU | NVIDIA Tesla T4 (15 GB VRAM) |
| CUDA | 12.8 |
| OS | Google Colab (Ubuntu) |

### 소프트웨어 (Python 3.9 가상환경)

| 패키지 | 버전 | 역할 |
|--------|------|------|
| Python | 3.9 | Colab 기본 Python 3.13과 DGL 호환 불가 → 별도 venv 필요 |
| PyTorch | 1.10.2+cu113 | 딥러닝 프레임워크 |
| DGL | 0.7.2 | 그래프 신경망 라이브러리 |
| biopandas | 0.2.8 | PDB 파일 파싱 |
| POT | 0.8.2 | Optimal Transport 계산 |
| rdkit | - | 분자 처리 |
| dgllife | 0.2.8 | 생명과학 그래프 모델 |

---

## 환경 설정

### Step 1 — GPU 런타임 설정

Colab 상단: `런타임 → 런타임 유형 변경 → T4 GPU`

```python
!nvidia-smi
```

### Step 2 — 저장소 Clone

```python
!git clone https://github.com/octavian-ganea/equidock_public.git
%cd equidock_public
```

### Step 3 — Python 3.9 가상환경 구축

```python
# Python 3.9 설치
!apt-get install -y python3.9 python3.9-venv python3.9-dev -q
!python3.9 -m venv /opt/equidock_env
!/opt/equidock_env/bin/pip install --upgrade pip -q

# 패키지 설치
!/opt/equidock_env/bin/pip install torch==1.10.2+cu113 \
    -f https://download.pytorch.org/whl/torch_stable.html -q
!/opt/equidock_env/bin/pip install dgl==0.7.2 \
    -f https://data.dgl.ai/wheels/repo.html -q
!/opt/equidock_env/bin/pip install \
    biopandas==0.2.8 POT==0.8.2 rdkit-pypi \
    dgllife==0.2.8 joblib==1.1.0 numpy==1.22.1 -q

# 확인
!/opt/equidock_env/bin/python -c "import torch, dgl; print(torch.__version__, torch.cuda.is_available(), dgl.__version__)"
!/opt/equidock_env/bin/python -c "import biopandas, ot, rdkit, dgllife; print('All packages OK')"
```

---

## 소스코드 수정 (트러블슈팅)

### 수정 1 — `inference_rigid.py` hardcoded `jean/` 경로 제거

```python
# 125~127번째 줄 주석 처리
with open('src/inference_rigid.py', 'r') as f:
    content = f.read()

content = content.replace(
    "    input_dir = './test_sets_pdb/jean/'\n"
    "    ground_truth_dir = './test_sets_pdb/jean/'\n"
    "    output_dir = './test_sets_pdb/jean_out/'",
    "    #input_dir = './test_sets_pdb/jean/'\n"
    "    #ground_truth_dir = './test_sets_pdb/jean/'\n"
    "    #output_dir = './test_sets_pdb/jean_out/'"
)

with open('src/inference_rigid.py', 'w') as f:
    f.write(content)
```

### 수정 2 — device CPU 강제 설정

```python
# DGL CUDA 버전 호환 문제 → CPU 강제 설정
content = content.replace(
    'args[\'device\'] = torch.device("cuda:0") if torch.cuda.is_available() else torch.device("cpu")',
    'args[\'device\'] = torch.device("cpu")  # forced CPU for DGL compatibility'
)
```

---

## 추론 실행

```python
%cd /content/equidock_public

!/opt/equidock_env/bin/python -m src.inference_rigid -data dips > output.log 2>&1
!tail -10 output.log
```

---

## 실행 환경 및 결과

### 추론 실행 환경

| 항목 | 값 |
|------|-----|
| 실행 모드 | CPU (DGL CUDA 버전 호환 문제로 강제 설정) |
| 총 테스트 구조 수 | 100개 (DIPS 테스트셋) |
| 평균 추론 시간 | **7.10초 / 구조** |
| 표준편차 | 5.76초 |
| 최소 추론 시간 | 1.49초 |
| 최대 추론 시간 | 26.15초 |
| 총 소요 시간 | 약 12분 |

### 논문 대비 속도 비교

| 방법 | 평균 추론 시간 (초) | 비고 |
|------|------------------|------|
| CLUSPRO | 10,475 | 기존 전통적 방법 |
| PATCHDOCK | 7,378 | 기존 전통적 방법 |
| EQUIDOCK (논문, GPU) | **3 ~ 5** | Tesla V100 |
| **EQUIDOCK (이번 실습, CPU)** | **7.10** | Colab T4 CPU 모드 |

> CPU 모드임에도 기존 방법 대비 **1,000배 이상 빠름**

### 출력 파일 구조

| 파일 종류 | 설명 |
|----------|------|
| `*_EQUIDOCK.pdb` | 기본 도킹 결과 (후처리 없음) |
| `*_EQUIDOCK_NO_CLASHES.pdb` | 충돌 제거(Intersection Loss) 후처리 결과 |

---

## RMSD 평가

```python
!/opt/equidock_env/bin/python -m src.test_all_methods.eval_pdb_outputset \
    -data dips -method equidock 2>&1 | tail -5
```

### 주요 지표

| 지표 | Median | Mean | Std |
|------|--------|------|-----|
| Complex RMSD (C-RMSD) | **13.30 Å** | 14.53 Å | ±7.13 |
| Interface RMSD (I-RMSD) | **10.19 Å** | 11.92 Å | ±7.01 |

### 논문 Table 1 대비 검증

| 지표 | 논문 결과 | 이번 실습 | 일치 여부 |
|------|---------|---------|---------|
| C-RMSD Median | 13.29 Å | **13.30 Å** | ✅ 재현 성공 |
| I-RMSD Median | 10.18 Å | **10.19 Å** | ✅ 재현 성공 |

> 논문 결과와 소수점 2자리까지 완벽히 재현됨

### RMSD 해석

| 항목 | 설명 |
|------|------|
| **C-RMSD** | 전체 복합체 CA 원자 좌표 편차 → 전체 도킹 자세 정확도 |
| **I-RMSD** | 결합 계면 근처 CA 원자 편차 → 결합 포켓 예측 정확도 |
| I-RMSD < C-RMSD | 전체 구조보다 **결합 포켓 예측이 더 정확** |
| Median 사용 이유 | 이상치에 robust — 일부 실패 케이스 영향 최소화 |

---

## 실습의 의의

### 입력 → 출력

| 항목 | 내용 |
|------|------|
| **입력** | 두 단백질 각각의 unbound 구조 (결합 전, 무작위 회전/이동 상태) |
| 입력 파일 | `receptor_unbound.pdb` + `ligand_unbound.pdb` |
| **추론** | EQUIDOCK End-to-End 예측 (평균 7초) |
| **출력** | 리간드가 수용체에 결합한 복합체 구조 예측 |
| 출력 파일 | `dock_EQUIDOCK.pdb` |

### 신약개발 관점에서의 의미

| 항목 | 기존 방법 | EQUIDOCK |
|------|----------|---------|
| 방식 | 수십만 개 후보 샘플링 후 scoring | 단 한 번의 End-to-End 예측 |
| 소요 시간 | 수 시간 | **7초** |
| 활용 단계 | Hit 발굴, 타겟 검증 | 동일 |

> **핵심:** 단백질-단백질 결합 구조를 빠르게 예측  
> → 어떤 단백질이 어디에 어떻게 결합하는지 파악  
> → 해당 결합을 **억제하거나 강화하는 약물 설계**의 출발점

---

## py3Dmol 시각화

```python
!pip install py3Dmol -q
import py3Dmol, os

def visualize_docking(case_name, title):
    pred_f = f'test_sets_pdb/dips_equidock_results/{case_name}_l_b_EQUIDOCK.pdb'
    gt_l_f = f'test_sets_pdb/dips_test_random_transformed/complexes/{case_name}_l_b_COMPLEX.pdb'
    gt_r_f = f'test_sets_pdb/dips_test_random_transformed/complexes/{case_name}_r_b_COMPLEX.pdb'

    with open(pred_f) as f: pred_pdb = f.read()
    with open(gt_l_f) as f: gt_pdb = f.read()
    with open(gt_r_f) as f: rec_pdb = f.read()

    print(title)
    view = py3Dmol.view(width=800, height=500)
    view.addModel(rec_pdb, 'pdb')
    view.setStyle({'model': 0}, {'cartoon': {'color': 'gray', 'opacity': 0.7}})
    view.addModel(gt_pdb, 'pdb')
    view.setStyle({'model': 1}, {'cartoon': {'color': 'green'}})
    view.addModel(pred_pdb, 'pdb')
    view.setStyle({'model': 2}, {'cartoon': {'color': 'red'}})
    view.zoomTo()
    view.show()

# 대표 케이스 (median 근처)
visualize_docking('a9_1a95.pdb1_3.dill', '대표 케이스 — C-RMSD ≈ 13.8 Å')

# 최선 케이스
visualize_docking('sk_3sk2.pdb1_0.dill', '✅ 최선 케이스 — C-RMSD = 3.09 Å')

# 최악 케이스
visualize_docking('dm_1dm3.pdb1_4.dill', '❌ 최악 케이스 — C-RMSD = 36.94 Å')
```

### 색상 범례

| 색상 | 구조 | 의미 |
|------|------|------|
| 🔘 회색 | 수용체 (Receptor) | Ground truth 결합 파트너 |
| 🟢 초록 | Ground truth 리간드 | 실험으로 밝혀진 실제 결합 위치 |
| 🔴 빨간 | EQUIDOCK 예측 리간드 | AI가 예측한 결합 위치 |

### 케이스별 RMSD

| 케이스 | PDB ID | C-RMSD | 관찰 |
|--------|--------|--------|------|
| 최선 | sk_3sk2 | **3.09 Å** | 예측(빨강)과 실제(초록) 거의 완벽히 일치 |
| 대표 | a9_1a95 | 13.83 Å | 결합 포켓 위치는 맞으나 말단 루프 편차 |
| 최악 | dm_1dm3 | **36.94 Å** | 예측이 실제 결합 위치와 크게 벗어남 |

---

## 트러블슈팅

| 문제 | 원인 | 해결 |
|------|------|------|
| Colab Python 3.13 + DGL 호환 불가 | `collections.Mapping` Python 3.10에서 제거 | Python 3.9 venv `/opt/equidock_env` 생성 |
| `torch==1.10.2` 설치 불가 | Colab PyPI에서 구버전 제거 | 가상환경에서 `torch_stable.html` 직접 설치 |
| `FileNotFoundError: jean/` | `inference_rigid.py` hardcoded 경로 | 125~127번째 줄 주석 처리 |
| DGL CUDA GPU 에러 | DGL 0.7.2 CPU 버전 + CUDA 불일치 | device를 CPU로 강제 설정 |
| `POT` 빌드 실패 | Python 3.9 + 최신 POT 소스 빌드 오류 | `POT==0.8.2` 버전 지정 |

---

## 명시적 한계

| 한계 | 설명 |
|------|------|
| Rigid-body 가정 | 결합 시 단백질 구조 변화 (Induced fit) 미반영 |
| C-RMSD 13.3 Å | 원자 수준 정밀도는 낮음 — 결합 포켓 탐색 용도 |
| unbound → bound | 실제로는 결합 전/후 구조가 다름 (Conformational flexibility) |
| CPU 모드 | DGL CUDA 호환 문제로 GPU 미사용 → 논문 대비 약 2배 느림 |
| DIPS 테스트셋만 | DB5.5 테스트셋 추론 미수행 |

---

## 참고문헌

1. Ganea OE, et al. **Independent SE(3)-Equivariant Models for End-to-End Rigid Protein Docking.** *arXiv:2111.07786* (2021).
2. Berman HM, et al. **The Protein Data Bank.** *Nucleic Acids Res.* 28(1), 235–242 (2000).

---

## 관련 링크

| 리소스 | URL |
|--------|-----|
| 논문 | https://arxiv.org/pdf/2111.07786 |
| 소스코드 | https://github.com/octavian-ganea/equidock_public |
| DIPS 데이터셋 | https://github.com/drorlab/DIPS |
| DB5.5 데이터셋 | https://zlab.umassmed.edu/benchmark/ |
| Google Colab | https://colab.research.google.com |
