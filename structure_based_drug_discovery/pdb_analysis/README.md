# Protein Data Bank (PDB) 구조 분석 실습

> **LAIDD — 단백질 구조 기반 약물탐색**  
> 강좌: 신약개발을 위한 단백질 구조 예측 및 상호작용 예측 (2강)  
> 실습 환경: Google Colab

---

## 실습 목표

1. PDB 파일에서 단백질 서열(SEQRES) 추출 및 구조 불완전성 분석
2. 단백질(ATOM)과 리간드(HETATM) 분리
3. 3D 좌표 → Distance Map / Contact Map 시각화

---

## 분석 대상

| 항목 | 값 |
|------|-----|
| PDB ID | [3HMM](https://www.rcsb.org/structure/3HMM) |
| 단백질 | TGF-β receptor type 1 (TGFR1) kinase domain |
| 실험 방법 | X-ray crystallography |
| 분석 체인 | Chain A |

---

## 환경 설정

```python
# Google Colab 기준
!pip install biopython
!pip install py3Dmol
!pip install biotite
```

**사용 라이브러리:**

| 라이브러리 | 용도 |
|-----------|------|
| `Biopython` | PDB 파싱, 서열 추출, 구조 분리 |
| `NumPy` | CA 좌표 행렬 연산 |
| `SciPy` | 유클리디안 거리 행렬 계산 |
| `Matplotlib` | Distance Map 시각화 |
| `py3Dmol` | Colab 내 인터랙티브 3D 구조 뷰어 |
| `Biotite` | Contact Map 계산 |

---

## 실습 1 — PDB 서열 추출 및 구조 불완전성 분석

### 1-1. SEQRES 서열 추출 (원래 전체 서열)

```python
from Bio import SeqIO
import os

pdb_file = 'pdb/3HMM.pdb'

for record in SeqIO.parse(pdb_file, "pdb-seqres"):
    chain_id = record.annotations['chain']
    sequence = record.seq
    print(f">Chain {chain_id} (길이: {len(sequence)} aa)")
    print(sequence)
```

**`pdb-seqres` 포맷:** PDB 파일 내 `SEQRES` 레코드를 파싱하여 실험에 사용된 단백질의 원래 전체 서열을 반환합니다.

### 1-2. ATOM 기반 서열 추출 (실제 관측된 서열)

```python
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', pdb_file)
ppb = PPBuilder()

for idx, pp in enumerate(ppb.build_peptides(structure)):
    sequence = pp.get_sequence()
    print(f">Fragment {idx+1} (길이: {len(sequence)} aa)")
    print(sequence)
```

**PPBuilder:** 3D 좌표(`ATOM` 레코드)가 실제로 존재하는 잔기만으로 폴리펩타이드 체인을 재구성합니다.

### 1-3. 두 서열 대조 — 누락 잔기 탐색

```python
# SEQRES와 ATOM 서열을 문자열 매칭으로 대조
# 전체 서열 안에서 관측 조각(Fragment)의 시작 위치를 찾아
# 공백 구간을 누락(missing) 잔기로 판정
frag_start_index = full_sequence.find(frag, current_search_index)
if frag_start_index > current_search_index:
    missing_seq = full_sequence[current_search_index:frag_start_index]
```

### 분석 결과 (3HMM Chain A)

| 항목 | 값 | 해석 |
|------|-----|------|
| SEQRES 전체 서열 길이 | 303 aa | 실험에 사용된 원래 단백질 |
| ATOM 관측 서열 | 293 aa (Fragment 2개) | 3D 좌표가 실제로 존재하는 잔기 |
| 누락 잔기 수 | 10 aa | X선에 좌표가 찍히지 않은 유연한 구간 |
| 루프 영역 누락 | 5 aa (`PNHRV`) | Fragment 1-2 사이 연결 루프 |
| C-말단 꼬리 누락 | 5 aa (`EGIKM`) | C-말단 자유 말단부 |

**구조 불완전성 원인:**

- **말단 꼬리 (N/C-말단):** 한쪽 끝이 허공을 향해 열려 있어 구조적 제약 없이 자유롭게 움직임 → X선 회절 시 평균 전자 밀도 불명확
- **루프 영역:** α-helix 및 β-sheet를 연결하는 유연한 관절 구조로, 표적 결합 시 형태 변화(Induced fit)가 필수적이어서 본질적으로 유동적

> **한계:** 누락 잔기 탐색은 문자열 정합(string matching)에 의존하므로, 서열 반복(repeat) 구간에서 오매핑 가능성 존재.

---

## 실습 2 — 단백질 / 리간드 분리

```python
from Bio.PDB import PDBParser, Select, PDBIO

# 단백질(표준 아미노산 ATOM)만 선택
class ProteinSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0] == ' '   # ' ' = 표준 아미노산

# 리간드(HETATM, 물 분자 제외)만 선택
class LigandSelect(Select):
    def accept_residue(self, residue):
        return residue.id[0].startswith('H_')  # 'H_' = 헤테로 원자 리간드
                                                # 'W'  = 물 분자 → 자동 제외
```

**출력 파일:**

| 파일 | 내용 |
|------|------|
| `3HMM_protein.pdb` | 표준 아미노산만 포함 |
| `3HMM_ligand.pdb` | 결합된 리간드(약물)만 포함 |

**PDB 잔기 ID 체계:**

| `residue.id[0]` | 의미 |
|-----------------|------|
| `' '` (공백) | 표준 아미노산 (ATOM) |
| `'H_XXX'` | 헤테로 원자 리간드 (HETATM) |
| `'W'` | 물 분자 (HOH) |

---

## 실습 3 — 3D 좌표 → Distance Map / Contact Map

### 3-1. Distance Map (전체 거리 히트맵)

```python
import numpy as np
from scipy.spatial import distance_matrix
import matplotlib.pyplot as plt

# Chain A의 CA 좌표 추출
ca_coords = []
for residue in structure[0]['A']:
    if residue.id[0] == ' ' and 'CA' in residue:
        ca_coords.append(residue['CA'].get_coord())

coords_array = np.array(ca_coords)  # (N, 3) 행렬

# 유클리디안 거리 행렬 계산: O(N²) 연산
dist_matrix = distance_matrix(coords_array, coords_array)

plt.imshow(dist_matrix, cmap='hot_r', origin='lower')
plt.colorbar(label='Distance (Å)')
```

**알고리즘 복잡도:** `scipy.spatial.distance_matrix`는 O(N²) 메모리와 연산을 요구합니다. N > 1,000 잔기 이상에서는 메모리 병목 발생 가능 → `cdist` 또는 KD-tree 기반 방법 고려 필요.

**Distance Map 해석:**

- **대각선 (거리 0Å):** 자기 자신과의 거리
- **대각선 근방 밝은 덩어리:** 국소적으로 조밀하게 뭉친 구조적 도메인(Domain) 경계
- **원거리 밝은 점:** 서열상 멀지만 3D 공간에서 인접한 잔기 쌍 → 단백질 3차 접힘 정보

### 3-2. 인터랙티브 3D 구조 뷰어 (py3Dmol)

```python
import py3Dmol

with open("3HMM_protein.pdb") as f:
    pdb_data = f.read()

view = py3Dmol.view(width=800, height=500)
view.addModel(pdb_data, 'pdb')
view.setStyle({'model': -1}, {'cartoon': {'color': 'spectrum'}})
view.zoomTo()
view.show()
```

Colab 내에서 마우스로 줌인/아웃 및 회전 가능한 인터랙티브 3D 뷰어를 제공합니다.

### 3-3. Contact Map (8Å 기준 이진 접촉 지도)

```python
import biotite.structure.io as strucio
import biotite.structure as struc

protein = strucio.load_structure("3HMM_protein.pdb")
ca_atoms = protein[protein.atom_name == "CA"]

# CellList: 공간 분할 기반 인접 탐색 (Distance Map보다 효율적)
cell_list = struc.CellList(ca_atoms, cell_size=8.0)
contact_matrix = cell_list.create_adjacency_matrix(8.0)  # 8Å 임계값

plt.imshow(contact_matrix, cmap='Greys', origin='lower')
```

**Distance Map vs Contact Map 비교:**

| 항목 | Distance Map | Contact Map |
|------|-------------|-------------|
| 출력값 | 연속적 거리 (Å) | 이진값 (0/1) |
| 임계값 | 없음 | 8Å (표준) |
| 계산 방법 | `scipy.distance_matrix` O(N²) | `CellList` 공간 분할 (효율적) |
| 정보량 | 풍부 (전체 거리 분포) | 간결 (접촉 여부만) |
| 주요 활용 | 구조 전체 조망 | 2차/3차 구조 패턴 분석 |

**Contact Map 해석:**

| 패턴 | 해석 |
|------|------|
| 대각선 근방 연속 선형 패턴 | α-helix (나선형 인접 잔기 밀착) |
| 대각선과 평행/수직 선형 점열 | β-sheet (가닥 간 수소결합) |
| 대각선 바깥 원거리 점 | 3차 접힘 핵심 골조 — 활성 포켓 후보 |

---

## 실습 결과 요약

| 분석 항목 | 주요 발견 | 신약개발 관련성 |
|-----------|-----------|----------------|
| 서열 불완전성 | 303→293 aa, 루프+C말단 10 aa 누락 | 누락 루프가 결합 포켓일 가능성 |
| 단백질/리간드 분리 | ATOM/HETATM 기준 분리 성공 | 도킹 시뮬레이션 입력 준비 |
| Distance Map | 도메인 경계 시각화 | 구조적 도메인 단위 약물 타겟 탐색 |
| Contact Map | α-helix/β-sheet 패턴, 원거리 접촉 | 활성 포켓 후보 잔기 쌍 식별 |

---

## 명시적 한계

| 한계 | 설명 |
|------|------|
| 단일 구조 분석 | 하나의 PDB 구조만 분석 — 동적 거동 미반영 |
| 정적 스냅샷 | X선 결정 구조는 결정 환경의 정적 상태 → 생리적 조건과 차이 가능 |
| 누락 잔기 미모델링 | 루프 영역은 homology modeling 또는 MD로 별도 보완 필요 |
| CA 원자 기준 | CA만으로 계산한 거리/접촉은 실제 원자 수준 상호작용의 근사값 |

---

## 파일 구조

```
pdb_analysis/
├── README.md                  ← 이 파일
├── pdb_prac1.py               ← 전체 실습 코드
└── pdb/
    └── 3HMM.pdb               ← 실습 대상 PDB 파일
```

---

## 참고문헌

1. Berman HM, et al. **The Protein Data Bank.** *Nucleic Acids Res.* 28(1), 235–242 (2000). https://www.rcsb.org
2. Cock PJ, et al. **Biopython: freely available Python tools for computational molecular biology.** *Bioinformatics* 25(11), 1422–1423 (2009).
3. Cock PJ, et al. **The Sanger FASTQ file format for sequences with quality scores, and the Solexa/Illumina FASTQ variants.** *Nucleic Acids Res.* (2010).

---

## 관련 링크

| 리소스 | URL |
|--------|-----|
| 3HMM 구조 (RCSB PDB) | https://www.rcsb.org/structure/3HMM |
| Biopython 문서 | https://biopython.org/wiki/Documentation |
| py3Dmol | https://github.com/3dmol/3Dmol.js |
| Biotite | https://www.biotite-python.org |
