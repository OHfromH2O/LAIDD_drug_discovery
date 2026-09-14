# PBI Toolkit을 활용한 단백질-리간드 도킹 전처리 실습

> **LAIDD — 단백질 구조 기반 약물탐색**  
> 강좌: 도킹 프로그램 사용 실습  
> 실습 환경: **Google Colab** (전체 실행) + **PyMOL Windows** (구조 시각화)  
> 실습 파일: [`pdb_prac2.ipynb`](pdb_prac2.ipynb)

---

## 실습 목표

타겟 단백질 **TGFR1_HUMAN** (TGF-β receptor type 1)에 대해 다음을 수행한다:

1. UniProt DB로부터 타겟 단백질의 서열 및 관련 PDB 목록 자동 수집
2. 도메인 기준 PDB 필터링 및 다운로드
3. Mutation 확인 및 구조 정렬
4. Chain 분리 및 단백질-리간드 분리
5. 포켓 근처 리간드 식별 및 Pharmacophore 탐색
6. Docking을 위한 단백질/리간드 전처리 및 Docking Box 계산

---

## 분석 대상

| 항목 | 값 |
|------|-----|
| 타겟 단백질 | TGFR1_HUMAN (TGF-β receptor type 1 kinase) |
| UniProt ID | [P36897](https://www.uniprot.org/uniprot/P36897) |
| 분석 도메인 | Kinase domain (잔기 205-495) |
| Reference 구조 | [3HMM](https://www.rcsb.org/structure/3HMM) Chain A |
| Template 리간드 | 3HMMA_855 |

---

## 실습 환경 구성

| 환경 | 용도 | 이유 |
|------|------|------|
| **Google Colab** | 전체 분석 (Python + bash 스크립트) | Ubuntu 기반 → `.sh` 직접 실행 가능 |
| **PyMOL (Windows)** | 구조 시각화 및 정렬 | GUI 필요 |

> **WSL2 conda를 사용하지 않는 이유:** `openbabel`, `pdbfixer`, `rdkit`이 WSL2 conda solver와 버전 충돌 → `Solving environment` 무한 대기 발생. Colab에서 동일 패키지 정상 설치됨.

---

## 워크플로우 전체 개요

```
[Colab]  패키지 설치 및 PBI toolkit 설정
      ↓
[Colab]  UniProt DB 파싱 (uniprot_dict.py make_pickle)
      ↓
[Colab]  타겟 단백질 서열 + PDB 목록 추출 (uniprot_dict.py pdb)
      ↓
[Colab]  Kinase domain 포함 PDB 필터링 (filter_pdb_list.py)
      ↓
[Colab]  PDB 파일 다운로드 (dw_pdb.py)
      ↓
[Colab]  Mutation 확인 + 리간드 종류 확인 (find_mutation_pdb.py)
      ↓
[PyMOL] 구조 시각화 및 정렬 (fetch + alignto 3HMM)
      ↓
[Colab]  Chain 분리 → TGFR1 해당 chain 선택
      ↓
[Colab]  TMalign 기반 3D 정렬 (rot.sh)
      ↓
[Colab]  단백질 / 리간드 분리 (split_ligand.sh)
      ↓
[Colab]  포켓 근처 리간드 필터링 (dist.sh)
      ↓
[PyMOL] 단백질-리간드 상호작용 분석
      ↓
[Colab]  HOH 추출 / Ligand fix / Protein fix
      ↓
[Colab]  Docking Box 계산 (auto_box.sh)
```

---

## Step 0 — Colab 환경 설정

### 0-1. Google Drive 마운트

```python
from google.colab import drive
drive.mount('/content/drive')
```

### 0-2. 패키지 설치

```python
# apt 패키지
!apt-get install -y openbabel libopenbabel-dev cmake -q

# pip 패키지
!pip install rdkit pdbfixer biopython openbabel-wheel -q

# 설치 확인
from rdkit import Chem
from openbabel import openbabel
print("rdkit OK")
print("openbabel OK")
```

### 0-3. PBI toolkit 설치 및 작업 폴더 설정

```python
import os

# 작업 폴더로 이동
os.chdir("/content/drive/MyDrive/.../practice_2")

# PBI.zip 압축 해제 및 script 복사
!unzip PBI.zip -d .
!cp -rf PBI/script TGFR1

# PBI 패키지 설치 (setup.py 기반)
!pip install -e PBI/ -q
!which uniprot_dict.py   # /usr/local/bin/uniprot_dict.py 확인

# 작업 폴더로 이동
os.chdir(".../TGFR1")
```

> **`pip install -e PBI/`가 필요한 이유:** `uniprot_dict.py`, `filter_pdb_list.py` 등 실행 스크립트가 `PBI/bin/`에 있으므로, `setup.py` 기반 설치로 PATH에 등록해야 Colab 어디서든 호출 가능.

---

## Step 1 — UniProt DB 구축

```python
# UniProt DB 저장 폴더 생성 (Drive에 영구 보존)
!mkdir -p "/content/drive/MyDrive/Database/uniprot"

# Swiss-Prot DB 다운로드 (~600MB, 최초 1회)
!wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.dat.gz \
    -P "/content/drive/MyDrive/Database/uniprot/" --show-progress

# DB 경로 설정 및 파싱
!echo "db_dir=/content/drive/MyDrive/Database/uniprot" > db_dir.txt
!uniprot_dict.py make_pickle
```

**생성 파일:**

| 파일 | 내용 |
|------|------|
| `uniprot.pkl*` | 전체 Swiss-Prot DB pickle |
| `uniprot_human.pkl` | 인간 단백질만 필터링된 pickle |

> **한계:** Swiss-Prot은 수동 검증 단백질만 포함 (TrEMBL 미포함).

---

## Step 2 — 타겟 단백질 정보 추출

```python
!uniprot_dict.py pdb TGFR1_HUMAN
!cat P36897.fasta
!cat P36897_pdb_list.txt
```

**결과:** P36897.fasta (522 aa canonical 서열), 46개 PDB ID 목록 생성

---

## Step 3 — Kinase Domain 포함 PDB 필터링

```python
!filter_pdb_list.py P36897_pdb_list.txt -r 205-495 > P36897_pdb_list_filter.txt
!cat P36897_pdb_list_filter.txt
```

**결과:** 46개 → **38개** (X-ray 구조, kinase domain 포함)

- NMR (`2L5S`), EM (`9FK5`) 구조 제외
- 세포외 도메인만 포함 구조 제외

> **한계:** 해당 영역이 일부 누락된 경우도 통과 가능.

---

## Step 4 — PDB 파일 다운로드

```python
!dw_pdb.py P36897_pdb_list_filter.txt pdb
!ls pdb
```

**결과:** 38개 중 **35개 다운로드 완료**

- 실패 (HTTP 404): `8YHF`, `8YHL`, `9J9D`, `9F6X` — 최신 PDB ID로 URL 형식 변경됨 → 무시

---

## Step 5 — Mutation 확인 및 리간드 종류 확인

```python
!find_mutation_pdb.py P36897_pdb_list_filter.txt P36897.fasta pdb > list_mutation.txt
!cat list_mutation.txt
```

**결과 분류:**

| 분류 | PDB 목록 | 수 |
|------|----------|-----|
| Wild type + 리간드 있음 | 1RW8, 1VJY, 2WOT, 2WOU, 3GXL, 3HMM, 3TZM, 4X0M, 4X2F, 4X2G, 4X2J, 4X2K, 5FRI, 5QIM, 5USQ | 15개 |
| Wild type + 리간드 없음 | 1B6C, 1IAS, 4X2N, 5E8S (SO4만 존재) | 4개 |
| Mutation 구조 | 5E8T, 5E8U, 5E8W, 5E8X, 5E8Z, 5E90, 5QIK, 5QIL, 5QTZ, 5QU0, 6B8Y (R199-, T204D 포함) | 11개 |
| 다운로드 실패 | 8YHF, 8YHL, 9F6X, 9J9D | 4개 |

> **주의:** `R199-`, `T204D`는 결정학 목적 engineered mutation으로, 임상적 내성 돌연변이와 구분 필요.

---

## Step 6 — PyMOL 구조 정렬 (Windows)

Wild type + 리간드 있는 15개 구조를 PyMOL에서 직접 fetch 후 정렬:

```python
# PyMOL 콘솔에서 실행
fetch 3HMM 1RW8 1VJY 2WOT 2WOU 3GXL 3TZM 4X0M 4X2F 4X2G 4X2J 4X2K 5FRI 5QIM 5USQ, async=0
alignto 3HMM
```

**정렬 결과 (RMSD vs 3HMM):**

| 구조 | RMSD (Å) | 구조 | RMSD (Å) |
|------|----------|------|----------|
| 3TZM | 0.316 | 5QIM | 0.457 |
| 5FRI | 0.337 | 2WOT | 0.513 |
| 2WOU | 0.510 | 4X0M | 0.550 |
| 3GXL | **1.393** | 나머지 | 0.5~0.7 |

전체 RMSD < 1.5Å → 동일한 결합 포켓 공유 확인.  
`3GXL`이 가장 큰 편차 → 루프 구조 차이 가능성.

---

## Step 7 — Chain 분리 및 TGFR1 Chain 선택

```python
# split_chain.py로 각 PDB를 chain별로 분리 (-e: water 제외)
import os
os.makedirs("chain", exist_ok=True)
with open("list_mutation.txt") as f:
    for line in f:
        if "does not exist" in line or not line.strip():
            continue
        pdb_id = line.split(";")[0]
        pdb_file = f"pdb/{pdb_id}.pdb"
        if os.path.exists(pdb_file):
            !split_chain.py -i {pdb_file} -d chain -e

# TGFR1 해당 chain 선택
!python gen_chain_list.py -i list_mutation.txt -o list_chain.txt
!cat list_chain.txt
```

**결과:** 34개 PDB의 TGFR1 chain 선택 완료  
(`1B6C`→B, 나머지→A)

---

## Step 8 — Reference 구조 기준 3D 정렬

```python
# TMalign 설치 (Colab에 없으므로 소스에서 빌드)
!wget https://zhanggroup.org/TM-align/TMalign.cpp -O TMalign.cpp -q
!g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp
!mv TMalign /usr/local/bin/TMalign

# 3HMMA.pdb 기준 정렬
!bash rot.sh list_chain.txt chain/3HMMA.pdb
!ls chain/*rotate.pdb | head -5
```

**결과:** 34개 `*_rotate.pdb` 생성

**Step 6과의 차이:**

| 단계 | 도구 | 목적 |
|------|------|------|
| Step 6 | PyMOL `alignto` | 시각적 확인 |
| Step 8 | TMalign (`rot.sh`) | 도킹용 좌표계 수치 통일 |

---

## Step 9 — 단백질 / 리간드 분리

```python
!bash split_ligand.sh list_chain.txt

# 리간드 파일 목록 생성 (HOH 제외)
!ls chain/?????_???.pdb | grep -v HOH > list_ligand.txt
```

**결과:** receptor `*_rotate.pdb` + 리간드 `*_XXX.pdb` 분리 완료

> **주의:** `split_chain.py -e` 옵션으로 water가 제거됨 → HOH 파일은 Step 12에서 원본 PDB에서 별도 추출 필요.

---

## Step 10 — 포켓 근처 리간드 식별

```python
!bash dist.sh list_ligand.txt chain/3HMMA_855.pdb > list_ligand_new.txt
!cp list_ligand_new.txt list_ligand_select.txt
!cat list_ligand_new.txt
```

**결과:** **28개 리간드** 선별 (3HMMA_855와 동일 결합 포켓)

| 컬럼 | 의미 |
|------|------|
| 원자수 | 리간드 원자 개수 |
| RMSD | template 중심 대비 이동 거리 (Å) |
| 최대거리 | 리간드 최대 반경 (Å, 범위: 10-18Å) |

---

## Step 11 — 단백질-리간드 상호작용 분석 (PyMOL)

```python
# PyMOL Windows에서 실행
load 3HMMA_rotate.pdb
load 3HMMA_855.pdb
# Action → Preset → Ligand Sites → Cartoon
```

**확인 항목:**

| 항목 | 설명 |
|------|------|
| Binding pose | 리간드(855)의 3D 결합 자세 |
| H-bond | 단백질-리간드 수소결합 잔기 |
| Hydrophobic contact | 소수성 접촉 잔기 |
| Bridging water | 단백질-리간드 양쪽과 수소결합하는 HOH |

---

## Step 12 — HOH 추출 / Ligand fix / Protein fix

### 12-1. HOH 파일 추출 (원본 PDB에서)

```python
# split_chain.py -e로 water가 제거됐으므로 원본 pdb에서 직접 추출
import os

with open("list_final.txt") as f:
    entries = [line.strip().split() for line in f if line.strip()]

for pdb_id, lig_id in entries:
    raw_pdb_id = pdb_id[:4]   # 3HMMA → 3HMM
    pdb_file = f"pdb/{raw_pdb_id}.pdb"
    hoh_file = f"chain/{pdb_id}_HOH.pdb"

    hoh_lines = [line for line in open(pdb_file)
                 if "HOH" in line and
                 (line.startswith("HETATM") or line.startswith("ATOM"))]

    with open(hoh_file, "w") as f:
        f.writelines(hoh_lines)
    print(f"{pdb_id}: {len(hoh_lines)}개 water")
```

### 12-2. Bridging water 식별

```python
!bash check_water.sh list_ligand_select.txt
!ls select/*_HOH.pdb > water_file_list.txt
```

**결과:** 28개 구조 모두 `select/*_HOH.pdb` 생성

### 12-3. Ligand fix (mmCIF 기준)

```python
!bash dw_lig_ref.sh list_final.txt   # RCSB에서 CIF 다운로드
!bash fix.sh list_final.txt           # fix/*_p.pdb 생성
```

### 12-4. Protein fix (pdbfixer)

```python
import os
os.makedirs("fix", exist_ok=True)

with open("list_final.txt") as f:
    entries = [line.strip().split() for line in f if line.strip()]

for pdb_id, lig_id in entries:
    receptor = f"select/{pdb_id}_receptor.pdb"
    hoh      = f"select/{pdb_id}_HOH.pdb"
    combined = f"select/{pdb_id}_receptor_HOH.pdb"
    fixed    = f"fix/{pdb_id}_receptor_HOH.pdb"

    # receptor + HOH 합치기
    with open(combined, "w") as out:
        out.write(open(receptor).read())
        if os.path.exists(hoh):
            out.write(open(hoh).read())

    !fix_protein.py -i {combined} -o {fixed}
```

**결과:** 28개 구조 `fix/*_receptor_HOH.pdb` 생성 완료

| 처리 | 목적 |
|------|------|
| Bridging water 식별 | 결합에 기여하는 물 분자 보존 |
| Ligand fix (mmCIF) | 수소, 결합 차수, 형식 전하 보정 |
| Protein fix (pdbfixer) | 누락 원자 보완, 수소 추가, 양성자화 |

> **한계:** pdbfixer 수소화는 pH 7 기준 근사값.

---

## Step 13 — Docking Box 계산

```python
!bash auto_box.sh list_final.txt > config.txt
!cat config.txt
```

**출력 결과:**

```
center_x=17.994
center_y=68.852
center_z=7.446
size_x=26.718
size_y=19.984
size_z=22.033
```

| 파라미터 | 값 (Å) | 설명 |
|----------|---------|------|
| `center_x/y/z` | 17.994 / 68.852 / 7.446 | 결합 포켓 중심 좌표 |
| `size_x` | 26.718 | Box X축 크기 |
| `size_y` | 19.984 | Box Y축 크기 |
| `size_z` | 22.033 | Box Z축 크기 |

> **한계:** Box 크기가 너무 작으면 결합 자세가 탐색 범위 밖으로 벗어남. 너무 크면 계산 시간 증가.

---

## 전체 실습 결과 요약

| Step | 작업 | 결과 |
|------|------|------|
| 0 | 환경 설정 | Colab + PBI toolkit 설치 완료 |
| 1 | UniProt DB 파싱 | `uniprot_human.pkl` 생성 |
| 2 | 타겟 단백질 추출 | `P36897.fasta`, 46개 PDB 목록 |
| 3 | Kinase domain 필터링 | 46개 → 38개 (X-ray 구조) |
| 4 | PDB 다운로드 | 35개 완료 (4개 404 오류) |
| 5 | Mutation 확인 | 15개 WT+리간드, 11개 mutation |
| 6 | PyMOL 정렬 | 모든 구조 RMSD < 1.5Å 확인 |
| 7 | Chain 분리 | 34개 TGFR1 chain 선택 |
| 8 | TMalign 3D 정렬 | 34개 `*_rotate.pdb` 생성 |
| 9 | 리간드 분리 | receptor + ligand 분리 완료 |
| 10 | 포켓 필터링 | 28개 리간드 선별 |
| 11 | 상호작용 분석 | 3HMMA_855 결합 포켓 확인 |
| 12 | Water/Ligand/Protein fix | 28개 구조 전처리 완료 |
| 13 | Docking Box | center=(17.994, 68.852, 7.446) |

---

## 전체 생성 파일 구조

```
pbi_toolkit_analysis/
├── README.md
├── pdb_prac2.ipynb                       ← 전체 실습 코드 (Colab .ipynb 변환)
├── PBI.zip                            ← PBI toolkit 원본

```

---

## 명시적 한계 및 트러블슈팅

| 이슈 | 원인 | 해결 |
|------|------|------|
| conda solver 무한 대기 | WSL2 Python 3.7 + 최신 패키지 충돌 | Colab으로 전환 |
| `TMalign` not found | Colab에 미설치 | 소스 빌드 후 `/usr/local/bin`에 배치 |
| HOH 파일 없음 | `split_chain.py -e`로 water 제거됨 | 원본 pdb에서 직접 추출 |
| PDB 404 오류 | 최신 PDB ID URL 형식 변경 | 해당 4개 구조 무시 |
| Colab 세션 초기화 | 세션 종료 시 파일 삭제 | Google Drive 마운트 필수 |
| PyMOL wildcard load 오류 | PyMOL 콘솔 단일 명령어 제한 | `fetch` 명령어로 직접 다운로드 |

---

## 참고문헌

1. Trott O, Olson AJ. **AutoDock Vina: improving the speed and accuracy of docking.** *J Comput Chem* 31(2), 455–461 (2010).
2. Berman HM, et al. **The Protein Data Bank.** *Nucleic Acids Res.* 28(1), 235–242 (2000).
3. Eastman P, et al. **PDBFixer.** https://github.com/openmm/pdbfixer (2013).
4. Zhang Y, Skolnick J. **TM-align: a protein structure alignment algorithm based on the TM-score.** *Nucleic Acids Res.* 33(7), 2302–2309 (2005).

---

## 관련 링크

| 리소스 | URL |
|--------|-----|
| TGFR1_HUMAN (UniProt) | https://www.uniprot.org/uniprot/P36897 |
| 3HMM (RCSB PDB) | https://www.rcsb.org/structure/3HMM |
| UniProt FTP | https://ftp.uniprot.org/pub/databases/uniprot/ |
| PDBFixer | https://github.com/openmm/pdbfixer |
| TMalign | https://zhanggroup.org/TM-align/ |
| PyMOL | https://pymol.org |
| Google Colab | https://colab.research.google.com |
