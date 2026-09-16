# 대규모 가상 탐색 — 멀티프로세싱 도킹 (`mdock_vina.py`)

> **LAIDD — 단백질 구조 기반 약물탐색**  
> 강좌: 도킹 프로그램 사용 실습 (홍승환, 한국제약바이오협회)  
> 실습 환경: **WSL2 Ubuntu** (12코어 LG 그램, Intel i5-1334U)  
> 실습 파일: [`prac4.ipynb`](prac4.ipynb) | [`mdock_vina.py`](mdock_vina.py)

---

## 실습 목표

1. Queue 기반 멀티프로세싱으로 135개 리간드 병렬 도킹 수행
2. 단일 도킹 대비 속도 향상 확인
3. 상위 Hit 화합물 선별

---

## 분석 대상

| 항목 | 값 |
|------|-----|
| 타겟 단백질 | TGFR1_HUMAN (TGF-β receptor type 1 kinase) |
| 수용체 | `pdb/3HMMA_receptor_HOH.pdbqt` |
| 도킹 리간드 수 | 135개 (ChEMBL TGFR1 관련 화합물) |
| 사용 코어 수 | 4 (`-p 4`, 총 12코어 중) |
| 도킹 도구 | Quick Vina 2 (`qvina02`) |

---

## 실습 환경

| 항목 | 사양 |
|------|------|
| CPU | Intel Core i5-1334U (13세대, 12코어) |
| RAM | 16 GB |
| OS | Windows 11 + WSL2 Ubuntu 22.04 |
| conda 환경 | `miniconda310/envs/docking` (Python 3.10) |

---

## 멀티프로세싱 원리

### Queue 기반 병렬화

단순 분할 방식은 화합물마다 계산 시간이 달라 CPU 낭비가 발생합니다.

```
단순 분할 방식 (비효율):
  Core 1: [분자1, 분자2, 분자3]  → 30초 완료 후 대기
  Core 2: [분자4, 분자5, 분자6]  → 60초 (큰 분자)
  Core 3: [분자7, 분자8, 분자9]  → 45초 완료 후 대기
```

```
Queue 방식 (효율적):
  Queue: [분자1, 분자2, ..., 분자135]
  Core 1: 분자1 완료 → Queue에서 분자4 즉시 가져옴
  Core 2: 분자2 완료 → Queue에서 분자5 즉시 가져옴
  Core 3: 분자3 완료 → Queue에서 분자6 즉시 가져옴
  → 항상 모든 코어가 바쁨 = CPU 100% 활용
```

### 코드 구조

```python
from multiprocessing import Queue, Process, Manager

def creator(q, data, num_sub_proc):
    """Master: 작업 목록을 Queue에 순서대로 쌓음"""
    for d in data:
        q.put((idx, d))
    for i in range(num_sub_proc):
        q.put('DONE')  # 종료 신호

def worker(q, param, return_dict):
    """Worker: Queue에서 작업을 하나씩 꺼내 도킹 실행"""
    while True:
        qqq = q.get()
        if qqq == 'DONE':
            break
        # gen_3d → pdb2pdbqt → docking → pdbqt2pdb
        return_dict[idx] = dock_score
```

---

## 도킹 파이프라인 (각 Worker 내부)

```
SMILES
  ↓  add_hydrogen(pH=7.4)       # openbabel
pH 보정 SMILES
  ↓  ligandtools.gen_3d()        # RDKit
3D PDB
  ↓  PDBtools.ligand_to_pdbqt() # obabel (gasteiger 전하)
PDBQT
  ↓  qvina02 --config            # Quick Vina 2 도킹
dock_*.pdbqt
  ↓  pdbqt_to_pdb_ref()         # bond 정보 복원
dock_*.pdb + score
```

---

## config.txt

```
receptor=pdb/3HMMA_receptor_HOH.pdbqt
center_x=17.994
center_y=68.852
center_z=7.446
size_x=26.718
size_y=19.984
size_z=22.033
cpu=1
num_modes=10
exhaustiveness=1
```

---

## 실행 명령어

```bash
# 환경 활성화
wsl
conda activate docking

# 기존 결과 삭제 (재실행 시)
rm -rf dock_mp docking_mp.txt

# 4코어 병렬 도킹 실행
~/miniconda310/envs/docking/bin/python mdock_vina.py \
    -v ~/docking_tools/qvina02 \
    -c config.txt \
    -s ligand_list.smi \
    -d dock_mp \
    -o docking_mp.txt \
    -p 4
```

**인자 설명:**

| 인자 | 값 | 설명 |
|------|-----|------|
| `-v` | `~/docking_tools/qvina02` | 도킹 실행 파일 |
| `-c` | `config.txt` | 도킹 박스 설정 |
| `-s` | `ligand_list.smi` | 리간드 SMILES 목록 |
| `-d` | `dock_mp` | 출력 폴더 |
| `-o` | `docking_mp.txt` | 점수 결과 파일 |
| `-p` | `4` | 병렬 프로세스 수 |

---

## 상위 10개 Hit 화합물 분석

```bash
~/miniconda310/envs/docking/bin/python -c "
results = []
with open('docking_mp.txt') as f:
    next(f)
    for line in f:
        parts = line.strip().split()
        if len(parts) >= 3 and parts[2] != 'None':
            mol_id = parts[0]
            best_score = float(parts[2].split(',')[0])
            results.append((mol_id, best_score))
results.sort(key=lambda x: x[1])
for i, (mol_id, score) in enumerate(results[:10], 1):
    print(f'{i:>4}  {mol_id:<18} {score:.1f}')
"
```

---

## 결과

### 상위 10개 Hit 화합물

| Rank | Mol ID | Best Score (kcal/mol) |
|------|--------|----------------------|
| 1 | **CHEMBL409356** | **-14.5** |
| 2 | CHEMBL260239 | -13.4 |
| 3 | CHEMBL406411 | -13.3 |
| 4 | CHEMBL567287 | -13.1 |
| 5 | CHEMBL436944 | -13.1 |
| 6 | CHEMBL585505 | -12.6 |
| 7 | CHEMBL260015 | -12.6 |
| 8 | CHEMBL204211 | -12.5 |
| 9 | CHEMBL206233 | -12.5 |
| 10 | CHEMBL566234 | -12.3 |

### 평가

- **CHEMBL409356 (-14.5 kcal/mol)** → 135개 중 가장 강한 결합 친화도
- 단일 도킹 Best Hit인 CHEMBL517068 (-10.5 kcal/mol)보다 **4.0 kcal/mol 더 강한 결합**
- 4코어 병렬 처리(`-p 4`)로 135개 리간드를 효율적으로 완료

> **Score 해석:** -14.5 kcal/mol은 매우 강한 결합 친화도  
> 일반적으로 -8 kcal/mol 이하면 strong binder로 간주

---

## 결과 파일 구조

```
3_4_multi_processing_for_virtual_screening/
├── prac4.py                     ← 전체 실습 코드
├── mdock_vina.py                ← 멀티프로세싱 도킹 스크립트
├── config.txt                   ← 도킹 박스 설정
├── ligand_list.smi              ← 135개 리간드 SMILES 목록
├── docking_mp.txt               ← 전체 도킹 점수 결과
├── pdb/                         ← 수용체 파일
│   └── 3HMMA_receptor_HOH.pdbqt
└── dock_mp/                     ← 도킹 결과 (mol_id 앞 7자로 분류)
    ├── CHEMBL3/
    │   └── CHEMBL341280.pdb
    ├── CHEMBL4/
    │   └── CHEMBL480742.pdb
    └── CHEMBL5/
        ├── CHEMBL517068.pdb
        └── CHEMBL519948.pdb
```

---

## 트러블슈팅

| 문제 | 원인 | 해결 |
|------|------|------|
| `dock_*.pdbqt` 파일 없음 | `PDBtools.ligand_to_pdbqt()` 실패 | obabel로 교체 |
| `prepare_ligand4` PDB 파싱 오류 | `gen_3d` 생성 PDB의 좌표 컬럼 포맷 불일치 | `pdbtools.py`에서 obabel 직접 호출로 변경 |
| `pdbqt_to_pdb_ref` KeyError | obabel pdbqt와 원본 pdb의 원자 인덱스 불일치 | try/except로 감싸서 경고 출력 후 계속 진행 |
| `miniconda3` 환경 우선 호출 | PATH 우선순위 문제 | `pdbtools.py`에 절대경로 (`/home/kucz11/miniconda310/...`) 하드코딩 |

### 핵심 수정 사항

**1. `pdbtools.py` — `ligand_to_pdbqt` 함수 교체**

```python
# 기존 (prepare_ligand4 → PDB 파싱 오류)
command = ['prepare_ligand4', '-l', pdb_file, '-o', pdbqt_file, '-U', 'nphs_lps']

# 수정 (obabel — PDB 포맷 무관하게 동작)
command = [
    '/home/kucz11/miniconda310/envs/docking/bin/obabel',
    pdb_file, '-o', 'pdbqt', '-O', pdbqt_file,
    '--partialcharge', 'gasteiger',
]
```

**2. `mdock_vina.py` — `pdbqt_to_pdb_ref` 예외 처리**

```python
# 기존
ligandtools.pdbqt_to_pdb_ref(dock_pdbqt_file, dock_pdb_file, pdb_file)

# 수정
try:
    ligandtools.pdbqt_to_pdb_ref(dock_pdbqt_file, dock_pdb_file, pdb_file)
except Exception as e:
    print(f'pdbqt_to_pdb_ref warning: {e}', flush=True)
```

---

## 명시적 한계

| 한계 | 설명 |
|------|------|
| `exhaustiveness=1` | 탐색 철저도 최소값 → 속도 우선, 정확도 낮음 |
| `pdbqt_to_pdb_ref` 실패 | obabel 변환 pdbqt와 원본 pdb 원자 인덱스 불일치로 일부 PDB 복원 실패 |
| 강체 도킹 | 수용체 유연성 미반영 |
| 단일 PC 병렬화 | NFS 기반 클러스터 분산 처리 미수행 (PVSdock) |
| Score 해석 | 도킹 점수는 결합 자유에너지 근사값 — 실험적 검증 필요 |

---

## 참고문헌

1. Alhossary A, et al. **Fast, Accurate, and Reliable Molecular Docking with QuickVina 2.** *Bioinformatics* 31(13), 2214–2216 (2015).
2. Python Software Foundation. **multiprocessing — Process-based parallelism.** https://docs.python.org/3/library/multiprocessing.html

---

## 관련 링크

| 리소스 | URL |
|--------|-----|
| Quick Vina 2 | https://github.com/QVina/qvina |
| PBI Toolkit | https://github.com/gicsaw/PBI |
| ChEMBL | https://www.ebi.ac.uk/chembl/ |
