# PDB 통계 분석 실습 — Cα 거리 히스토그램 & 라마찬드란 플롯

> **LAIDD — 단백질 구조 기반 약물탐색**  
> 강좌: 도킹 프로그램 사용 실습  
> 실습 환경: **Google Colab**  
> 실습 파일: [`4_ca_dist.py`](4_ca_dist.py) | [`5_rama.py`](5_rama.py)

---

## 실습 목표

PISCES 비중복 단백질 데이터셋을 활용하여 단백질 뼈대(backbone)의 기하학적 특성을 통계적으로 분석한다:

1. **Exercise 4:** Cα–Cα 거리 히스토그램으로 단백질 2차 구조 패턴 정량화
2. **Exercise 5:** 라마찬드란 플롯(φ–ψ 이면각 분포)으로 아미노산별 입체화학적 자유도 시각화

---

## 데이터셋

| 항목 | 값 |
|------|-----|
| 출처 | PISCES (Protein Culling Server) — 비중복 PDB 구조 큐레이션 |
| 형식 | `.pdb.gz` (gzip 압축 PDB) |
| 선택 기준 | 서열 유사도 < 25%, 해상도 < 2.5Å, R-factor < 0.25 |
| 분석 단위 | PDB ID + Chain ID (list.txt 기반) |

> **비중복 데이터셋을 사용하는 이유:** 서열이 유사한 단백질은 구조도 유사하므로, 중복 포함 시 특정 구조 유형이 과대 표현되어 통계 편향 발생. PISCES 큐레이션으로 이를 방지.

---

## 공통 환경 설정

```python
from google.colab import drive
drive.mount('/content/drive')

import os, numpy as np, gzip, matplotlib.pyplot as plt

data_dir = '/content/drive/MyDrive/.../practice_3/pisces/pdb'
list_file = data_dir + '/list.txt'
pdb_list = [x.strip().split() for x in open(list_file)]
```

---

## 공통 모듈 — PDB 파싱 (`read_pdb`)

두 실습 파일에서 공통으로 사용하는 PDB 파서입니다.

```python
main_atom_list = [' N  ', ' CA ', ' C  ', ' O  ']  # backbone 원자만 추출

def read_pdb(pdb_file):
    """
    gzip PDB 파일 파싱.
    반환: chain_dict, atom_coor_dict, modres_dict

    처리 항목:
      - REMARK 465: 누락 잔기(missing residues) 등록
      - MODRES: 변형 잔기(modified residue) → 표준 잔기 매핑
      - ATOM: 표준 아미노산 backbone 좌표 추출
      - HETATM: MODRES에 등록된 변형 잔기만 처리
      - alt_loc: 'A' 또는 공백만 허용 (alternate conformation 제거)
    """
```

**가정 및 한계:**

| 항목 | 내용 |
|------|------|
| 표준 아미노산 | 20종만 처리 (비표준은 MODRES 통해서만 허용) |
| Alternate conformation | 'A' 또는 공백만 사용 — 다중 구조 무시 |
| 누락 잔기 | REMARK 465 기반 등록하나 좌표 없음 → 거리 계산에서 자동 제외 |
| HETATM | MODRES 매핑된 변형 잔기만 포함 |

---

## Exercise 4 — Cα 거리 히스토그램 (`4_ca_dist.py`)

### 분석 원리

서열상 `i`번째와 `i+j`번째 잔기 (j=1~5)의 Cα 원자 간 유클리디안 거리를 수집하여 히스토그램으로 시각화합니다.

```python
def cal_single(chain_dict, atom_coor_dict, modres_dict, chain_id, ca_dist_dict):
    for i in range(num_keys):
        for j in range(1, 6):           # gap 1~5
            dist = np.linalg.norm(pos_ca_i - pos_ca_j)
            ca_dist_dict[j] += [dist]
```

### 결과 및 해석

#### Gap 1 (i, i+1) — 인접 잔기 거리

```
분포: 3.8 Å 부근 단일 예리한 봉우리
```

**해석:** 펩타이드 결합은 이중 결합 특성으로 회전 불가능한 **평면 구조(trans)** 로 고정됨. 3.0 Å 부근의 미세한 흔적은 극히 드문 **cis 결합** 을 의미.

---

#### Gap 2 (i, i+2) — Bimodal 분포

```
짧은 봉우리: ~5.4 Å  →  α-helix
넓은 봉우리: ~6.5 Å  →  β-strand (extended)
```

**해석:** 뼈대가 꺾이는 방식에 따라 2차 구조가 갈림.
- α-helix: 사슬이 둥글게 말려 가깝게 위치
- β-strand: 지그재그로 뻗어 멀리 위치

> β-sheet는 n=2로 가닥의 형태를 판별하고, 가닥 간 수소결합은 서열상 멀리 떨어진 잔기 간 원거리 상호작용이므로 Distance Matrix 교차 검증 필요.

---

#### Gap 3 (i, i+3) — α-helix 시그니처

```
좁고 높은 봉우리: ~5.0~5.5 Å  →  α-helix (3.6잔기/회전)
넓은 분포: 9~10 Å              →  non-helix
```

**해석:** α-helix는 3.6개 잔기마다 한 바퀴 회전 → i와 i+3이 같은 방향을 향해 가깝게 위치. 더 조밀한 **3₁₀-helix** 판별에도 활용됨.

---

#### Gap 4 (i, i+4) — α-helix 수소결합 주기

```
뾰족한 봉우리: ~6.0~6.5 Å  →  나선 1 Turn 수직 높이 (pitch)
넓은 분포: 12~14 Å           →  β-strand, random coil
```

**해석:** α-helix의 수소결합은 **i번째 C=O ↔ i+4번째 N-H** 사이에 형성됨. 6.0~6.5Å는 나선이 정확히 한 바퀴 완성된 수직 높이(1 Turn pitch).

---

#### Gap 5 (i, i+5) — 구조적 자유도 증가

```
α-helix 잔류 봉우리: ~8.5~9.0 Å
전체적으로 18 Å까지 넓고 평평하게 분산
```

**해석:** 서열 간격이 멀어질수록 뼈대의 기하학적 자유도가 기하급수적으로 증가하여 분포가 확산(diffuse)됨.

---

#### 프롤린(Proline) 특이 분석

```python
def cal_pro(...):
    # 다음 잔기가 PRO일 때만 gap 1 거리 계산
    if res_j_name != 'PRO':
        continue
```

```
결과: 3.0 Å (cis) 봉우리가 일반 gap 1 대비 현저히 증가
```

**해석:** 프롤린은 측쇄가 뼈대 질소(N)와 결합해 **5원환 고리(pyrrolidine ring)** 를 형성 → cis 펩타이드 결합 확률이 다른 아미노산 대비 압도적으로 높음.

---

## Exercise 5 — 라마찬드란 플롯 (`5_rama.py`)

### 분석 원리

각 잔기의 **φ(phi)** 와 **ψ(psi)** 이면각(dihedral angle)을 계산하여 산점도로 시각화합니다.

```python
def cal_torsional_angle(pos_1, pos_2, pos_3, pos_4):
    """
    4개 원자 좌표로 이면각 계산.
    법선 벡터(normal vector) 외적 기반.
    반환: 각도 (도, degree)
    """
    na = cal_normal_vector(pos_1, pos_2, pos_3)
    nb = cal_normal_vector(pos_2, pos_3, pos_4)
    # 부호 결정: cross product의 방향으로 +/-
    z = n32.dot(nd)
    angle2 = z * 180 * angle / np.pi
    return angle2

# φ (phi): C(i-1) - N - Cα - C
phi = cal_torsional_angle(pos_cm, pos_n, pos_ca, pos_c)

# ψ (psi): N - Cα - C - N(i+1)
psi = cal_torsional_angle(pos_n, pos_ca, pos_c, pos_np)
```

**한계:** 말단 잔기(N-terminal, C-terminal)는 φ 또는 ψ 계산 불가 → 자동 제외.

---

### 아미노산별 라마찬드란 플롯 결과 및 해석

#### 알라닌 (ALA) — 일반 아미노산의 표준

```
밀집 영역:
  좌측 상단 (φ < 0, ψ > 0): β-strand 영역
  좌측 하단 (φ < 0, ψ < 0): α-helix 영역
우측 (φ > 0): 거의 비어 있음
```

**해석:** 측쇄와 카보닐 산소 간의 **입체적 충돌(steric clash)** 로 우측 평면이 구조적으로 불안정. 좌상향 나선(left-handed helix) 등 극히 일부만 허용.

---

#### 글라이신 (GLY) — 가장 유연한 아미노산

```
분포: 평면 전체에 넓고 대칭적으로 분산
우측 (φ > 0) 포함: 다른 아미노산에서 금지된 영역까지 접근 가능
```

**해석:** 측쇄(R기)가 수소(H) 원자 하나뿐 → 입체적 충돌 없음. 단백질 뼈대에 **구조적 유연성** 을 부여하는 역할 (예: 루프, 턴 구조에 빈번하게 등장).

---

#### 프롤린 (PRO) — 가장 제약된 아미노산

```
분포: φ ≈ -50° ~ -90° 수직 띠로 극단적으로 제한
ψ: 비교적 자유로움
```

**해석:** 측쇄가 뼈대 N 원자와 공유결합하여 **5원환 고리** 형성 → N-Cα 결합 축 회전을 물리적으로 잠금. φ 각도가 좁은 범위에 갇힘. α-helix를 끊는 **helix breaker** 로 작용.

---

### 아미노산별 입체화학적 특성 요약

| 아미노산 | φ 허용 범위 | ψ 허용 범위 | 특성 |
|----------|------------|------------|------|
| ALA (일반) | 주로 < 0° | 전 범위 (좌측) | 표준적 제약 |
| GLY | 전 범위 | 전 범위 | 최대 유연성 |
| PRO | -50° ~ -90° (고정) | 비교적 자유 | 최대 제약, helix breaker |

---

## 알고리즘 한계

| 한계 | 설명 |
|------|------|
| 정적 구조 기반 | X선 결정 구조 스냅샷 — 동적 거동(MD 시뮬레이션) 미반영 |
| 이면각 계산 정밀도 | 법선 벡터 외적 기반 근사 — 좌표 노이즈에 민감 |
| Alternate conformation | 'A'만 사용 — 다중 구조 정보 손실 |
| PISCES 데이터셋 편향 | 결정화 가능한 단백질만 포함 — IDPs(내재적 무질서 단백질) 과소 표현 |
| Gap 5 이상 | 서열 간격 증가 시 구조 신호 소실 — 3차 구조 분석은 Distance Matrix 필요 |

---

## 파일 구조

```
pdb_statistics/
├── README.md
├── 4_ca_dist.py          ← Cα 거리 히스토그램 (Gap 1~5 + Proline 특이 분석)
├── 5_rama.py             ← 라마찬드란 플롯 (ALA, GLY, PRO)
└── pisces/
    └── pdb/
        ├── list.txt      ← PDB ID + Chain ID 목록
        └── [A-Z]/        ← PDB ID 첫 글자별 하위 폴더
            └── *.pdb.gz  ← gzip 압축 PDB 파일
```

---

## 참고문헌

1. Wang G, Dunbrack RL. **PISCES: a protein sequence culling server.** *Bioinformatics* 19(12), 1589–1591 (2003).
2. Ramachandran GN, Ramakrishnan C, Sasisekharan V. **Stereochemistry of polypeptide chain configurations.** *J Mol Biol* 7(1), 95–99 (1963).
3. Richardson JS. **The anatomy and taxonomy of protein structure.** *Adv Protein Chem* 34, 167–339 (1981).

---

## 관련 링크

| 리소스 | URL |
|--------|-----|
| PISCES | https://dunbrack.fccc.edu/pisces/ |
| RCSB PDB | https://www.rcsb.org |
| Ramachandran Plot 이론 | https://en.wikipedia.org/wiki/Ramachandran_plot |
