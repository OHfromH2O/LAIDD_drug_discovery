# 도킹 자동화 실습 — Bash / Python / Pharmacophore 분석

> **LAIDD — 단백질 구조 기반 약물탐색**  
> 강좌: 도킹 프로그램 사용 실습 (홍승환, 한국제약바이오협회)  
> 실습 환경: **WSL2 Ubuntu** + **PyMOL Windows**  
> 실습 파일: [`prac3.ipynb`](prac3.ipynb)

---

## 실습 목표

1. Bash 스크립트로 다수 리간드 도킹 자동화 (`dock_vina.sh`)
2. Python 스크립트로 동일한 파이프라인 자동화 (`dock_vina.py`)
3. SMARTS 기반 Pharmacophore 추출 (`find_pharmacophore_pdb.py`)
4. 도킹 결과의 Pharmacophore 충족 여부 판정 (`pharmacophore.py`)

---

## 도킹 자동화가 필요한 이유

신약개발에서 도킹은 단일 분자가 아닌 **수천~수백만 개의 후보 분자**를 대상으로 수행됩니다.

| 단계 | 규모 | 방법 |
|------|------|------|
| 초기 화합물 라이브러리 | 100만 ~ 1억 개 | 고속 가상 탐색 (HTVS) |
| 1차 필터링 | 수만 개 | 표준 도킹 |
| 2차 필터링 | 수백 개 | 정밀 도킹 + MD 시뮬레이션 |
| 실험 검증 | 수십 개 | 실제 합성 + 활성 측정 |

> ⚠️ 수동으로 하나씩 도킹하면 **100개만 해도 수일 소요** → 자동화 필수

---

## 분석 대상

| 항목 | 값 |
|------|-----|
| 타겟 단백질 | TGFR1_HUMAN (TGF-β receptor type 1 kinase) |
| 수용체 | `pdb/3HMMA_receptor_HOH.pdb` |
| Reference 리간드 | `pdb/3HMMA_855.pdb` |
| 도킹 리간드 | CHEMBL519948, CHEMBL480742, CHEMBL341280, CHEMBL517068 |
| Docking Box | center=(17.994, 68.852, 7.446), size=(26.718, 19.984, 22.033) |

---

## 전체 자동화 파이프라인

```
화합물 라이브러리 (SMILES)
        ↓  dock_vina.sh / dock_vina.py
  대규모 도킹 자동화
        ↓  find_pharmacophore_pdb.py
  Pharmacophore 패턴 추출 (SMARTS)
        ↓  pharmacophore.py
  공통 결합 특징 분석
        ↓
  Hit 화합물 선별 → 실험 검증
```

---

## 환경 설정

```bash
wsl
conda activate docking
cd "/mnt/c/Users/PC/.../3_3_docking_automation/TGFR1"
```

---

## 실습 1 — Bash 자동화 (`dock_vina.sh`)

### 스크립트 구조

```bash
#!/bin/bash
# 사용법: ./dock_vina.sh ligand_list.smi

list_file=$1
dock_dir=dock
mkdir -p $dock_dir

while read line
do
    array=($line)
    mol_id=${array[0]}
    smi=${array[1]}

    # 1. pH 7.4 protonation state 계산
    smi_p=$(obabel -:"$smi" -p7.4 -osmi 2>/dev/null | head -1)

    # 2. 3D conformer 생성
    gen3d.py -i "$smi_p" -o $dock_dir/$mol_id.pdb

    # 3. pdb → pdbqt 변환
    pdb2pdbqt.py -i $dock_dir/$mol_id.pdb -o $dock_dir/$mol_id.pdbqt -l

    # 4. 도킹 실행 (Quick Vina 2)
    qvina02 --config config.txt \
        --ligand $dock_dir/$mol_id.pdbqt \
        --out $dock_dir/dock_$mol_id.pdbqt

    # 5. pdbqt → pdb 변환 (bond 정보 복원)
    pdbqt2pdb_ref.py \
        -i $dock_dir/dock_$mol_id.pdbqt \
        -o $dock_dir/dock_$mol_id.pdb \
        -r $dock_dir/$mol_id.pdb

done < $list_file
```

> **트러블슈팅:** WSL에서 `obabel`, `qvina02` 등이 PATH에 없는 경우,  
> `sed`로 절대경로 (`~/miniconda310/envs/docking/bin/`) 로 일괄 치환 필요.  
> OneDrive 경로에서는 `sed -i`가 권한 문제로 실패 → Windows 메모장으로 직접 수정.

### ligand_list.smi

```
CHEMBL519948 Cc4cccc(c1nc(N)sc1c3ccc2nccnc2c3)n4
CHEMBL480742 Cc4nc(c1nc(N)[nH]c1c3ccc2nccnc2c3)ccc4F
CHEMBL341280 Brc4cccc(c1[nH]ncc1c2ccnc3ccccc23)n4
CHEMBL517068 CC(=O)Nc4nc(c1ccc(F)c(C)n1)c(c3ccc2ncnn2c3)[nH]4
```

### 실행

```bash
chmod +x dock_vina.sh
./dock_vina.sh ligand_list.smi
ls dock/
```

---

## 실습 2 — Python 자동화 (`dock_vina.py`)

### 주요 함수

```python
def add_hydrogen(smi0, pH=7.4):
    """SMILES에 pH 7.4 기준 수소 추가."""
    m = pybel.readstring("smi", smi0)
    m.OBMol.AddHydrogens(False, True, pH)
    return m.write("smi").strip()

def docking(vina, config_file, ligand_file, output_file):
    """Quick Vina 2 도킹 실행."""
    run_line = f'{vina} --config {config_file} --ligand {ligand_file} --out {output_file}'
    subprocess.check_output(run_line.split(), ...)
```

### 실행

```bash
~/miniconda310/envs/docking/bin/python dock_vina.py \
    ~/docking_tools/qvina02 \
    config.txt \
    ligand_list.smi \
    dock_py

ls dock_py/
```

> **트러블슈팅:** `ligandtools.gen_3d()` 에서 `timeout` 인자 미지원 오류 발생  
> → `dock_vina.py` 66번째 줄에서 `timeout=20` 인자 제거 필요 (Windows 메모장으로 수정)

### Bash vs Python 비교

| 항목 | dock_vina.sh (Bash) | dock_vina.py (Python) |
|------|--------------------|-----------------------|
| 언어 | Bash | Python |
| 가독성 | 단순, 선형 | 모듈화, 예외 처리 |
| 확장성 | 제한적 | 높음 (함수 추가 용이) |
| 오류 처리 | 제한적 | try/except로 개별 처리 |
| 적합한 용도 | 빠른 프로토타입 | 대규모 자동화 파이프라인 |

---

## 실습 3 — Pharmacophore 분석

### SMARTS란?

SMILES의 확장 언어로, **분자 내 특정 패턴**을 검색하는 쿼리 언어입니다.

| 표현 | 예시 | 의미 |
|------|------|------|
| SMILES | `Cc1cccc(n1)` | 분자 구조 표현 |
| SMARTS | `[nH]` | 방향족 N-H 패턴 검색 |
| SMARTS | `a1:a:a:a:a:a:1` | 6원환 방향족 고리 검색 |

📎 참고: https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

### Pharmacophore란?

**Pharmacophore** = 단백질과의 결합에 **필수적인** 리간드의 3D 화학적 특징 집합

| 특징 유형 | 해당 작용기 |
|-----------|------------|
| H-bond donor (수소결합 공여체) | N-H, O-H |
| H-bond acceptor (수소결합 수용체) | N, O |
| Hydrophobic (소수성 접촉) | 방향족, 지방족 |
| Cation (양이온성) | 양전하 질소 |

---

### Step 6.1 — find_pharmacophore_pdb.py

Reference 리간드 855의 **소수성 방향족 고리 중심 좌표** 추출

```bash
~/miniconda310/envs/docking/bin/python find_pharmacophore_pdb.py
```

**작동 원리:**

1. SMARTS 패턴 (`a1:a:a:a:a:a:1` 등)으로 방향족 고리 탐색
2. 고리를 구성하는 원자들의 평균 좌표 계산 → **pseudo-atom** (가상 원자) 생성
3. 방향족 고리는 개별 원자가 아닌 **π 전자구름 전체**로 단백질과 결합 (Pi-Pi stacking)
   → 고리 중심 1개 좌표로 단순화하여 연산 효율 극대화

**실행 결과:**

| # | X (Å) | Y (Å) | Z (Å) | 의미 |
|---|-------|-------|-------|------|
| 1 | 12.978 | 64.211 | 4.400 | 소수성 고리 중심 1 |
| 2 | 15.542 | 66.681 | 6.614 | 소수성 고리 중심 2 |
| 3 | 17.387 | 68.980 | 2.449 | 소수성 고리 중심 3 |
| 4 | 16.222 | 67.623 | 8.757 | 소수성 고리 중심 4 |

> 이 4개 좌표가 `pharmacophore.py`의 `pharmacophore_coor_ref_list`로 사용됨

---

### Step 6.2 — pharmacophore.py

도킹된 리간드가 Pharmacophore 조건을 충족하는지 판정

```bash
~/miniconda310/envs/docking/bin/python pharmacophore.py
```

**판정 기준:**
- cutoff = **2.5 Å**
- 새 리간드의 방향족 고리 중심이 reference 포인트 3개 이상과 2.5Å 이내 → prediction = **1** (Hit)

**실행 결과:**

| 리간드 | Best Score (kcal/mol) | Pharmacophore 충족 (Mode 1) | 최종 평가 |
|--------|----------------------|----------------------------|-----------|
| CHEMBL519948 | -9.0 | ✅ 1 | Hit |
| CHEMBL480742 | -9.6 | ✅ 1 | Hit |
| CHEMBL341280 | -10.3 | ✅ 1 | Hit |
| CHEMBL517068 | **-11.2** | ✅ 1 | **Best Hit** |

---

### Mode란?

도킹 프로그램은 리간드의 결합 자세를 **여러 개** 탐색합니다. 각 자세를 **Mode** 라고 합니다.

| Mode | 의미 |
|------|------|
| Mode 1 | 가장 낮은 score (최적 결합 자세) |
| Mode 2 | 두 번째로 낮은 score |
| Mode N | N번째 결합 자세 |

**CHEMBL517068 예시:**

| Mode | Score (kcal/mol) | Pharmacophore 충족 |
|------|-----------------|-------------------|
| 1 | -11.2 | ✅ 1 |
| 2 | -11.0 | ✅ 1 |
| 3 | -10.4 | ✅ 1 |
| 4 | -10.0 | ✅ 1 |
| 5~10 | -9.3 ~ -8.9 | ❌ 0 |

> Mode 1~4: score 높음 + pharmacophore 충족 → 결합 포켓 내 올바른 자세  
> Mode 5~: pharmacophore 미충족 → 비특이적 자세

---

## 결과 시각화 (PyMOL)

```bash
# WSL에서 (X server 필요)
pymol pdb/3HMMA_receptor_HOH.pdb pdb/3HMMA_855.pdb dock/dock_*.pdb
```

```python
# PyMOL Windows 콘솔에서
cd C:/.../3_3_docking_automation/TGFR1
load pdb/3HMMA_receptor_HOH.pdb
load pdb/3HMMA_855.pdb
load dock_py/dock_CHEMBL519948.pdb
load dock_py/dock_CHEMBL480742.pdb
load dock_py/dock_CHEMBL341280.pdb
load dock_py/dock_CHEMBL517068.pdb

hide everything
show cartoon, 3HMMA_receptor_HOH
show sticks, 3HMMA_855
show sticks, dock_CHEMBL519948
show sticks, dock_CHEMBL480742
show sticks, dock_CHEMBL341280
show sticks, dock_CHEMBL517068
color cyan, 3HMMA_receptor_HOH
color yellow, 3HMMA_855
zoom
```

---

## 결과 파일 구조

```
3_3_docking_automation/
└── TGFR1/
    ├── prac3.py                      ← 전체 실습 코드
    ├── dock_vina.sh                  ← Bash 자동화 스크립트
    ├── dock_vina.py                  ← Python 자동화 스크립트
    ├── find_pharmacophore_pdb.py     ← Pharmacophore 좌표 추출
    ├── pharmacophore.py              ← Pharmacophore 충족 판정
    ├── ligand_list.smi               ← 도킹 대상 리간드 목록
    ├── config.txt                    ← 도킹 박스 설정
    ├── pdb/                          ← 수용체 + reference 리간드
    ├── dock/                         ← Bash 자동화 도킹 결과
    │   ├── CHEMBL*.pdb
    │   ├── CHEMBL*.pdbqt
    │   ├── dock_CHEMBL*.pdb
    │   └── dock_CHEMBL*.pdbqt
    └── dock_py/                      ← Python 자동화 도킹 결과
        ├── CHEMBL*.pdb
        └── dock_CHEMBL*.pdb
```

---

## 명시적 한계

| 한계 | 설명 |
|------|------|
| 강체 도킹 | 수용체 유연성 미반영 (Induced fit 무시) |
| Pharmacophore 단순화 | 소수성 고리만 고려 — 수소결합, 전하 상호작용 미포함 |
| 소규모 검증 | 4개 리간드만 실습 — 실제 가상 탐색은 수만~수백만 개 |
| cutoff 고정 | 2.5Å cutoff는 경험적 값 — 타겟에 따라 최적화 필요 |
| WSL PATH 충돌 | miniconda3 (구 환경)과 miniconda310 (신 환경) 혼재 — 절대경로 명시 필요 |

---

## 트러블슈팅

| 문제 | 원인 | 해결 |
|------|------|------|
| `obabel: command not found` | conda 환경 PATH 미등록 | `/usr/local/bin`에 심볼릭 링크 |
| `gen_3d() unexpected keyword argument 'timeout'` | PBI 버전 불일치 | `dock_vina.py`에서 `timeout=20` 제거 |
| `sed -i` 권한 오류 | OneDrive 경로 제한 | Windows 메모장으로 직접 수정 |
| miniconda3 환경 호출 | PATH 우선순위 문제 | `~/miniconda310/envs/docking/bin/` 절대경로 명시 |

---

## 참고문헌

1. Trott O, Olson AJ. **AutoDock Vina.** *J Comput Chem* 31(2), 455–461 (2010).
2. Alhossary A, et al. **QuickVina 2.** *Bioinformatics* 31(13), 2214–2216 (2015).
3. Daylight Chemical Information Systems. **SMARTS Theory.** https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html

---

## 관련 링크

| 리소스 | URL |
|--------|-----|
| AutoDock Vina | https://github.com/ccsb-scripps/AutoDock-Vina |
| Quick Vina 2 | https://github.com/QVina/qvina |
| PBI Toolkit | https://github.com/gicsaw/PBI |
| SMARTS 문서 | https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html |
