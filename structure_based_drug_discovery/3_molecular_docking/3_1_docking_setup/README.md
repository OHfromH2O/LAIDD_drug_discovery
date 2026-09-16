# 도킹 프로그램 사용 실습 — 환경 설정 및 도킹 실행

> **LAIDD — 단백질 구조 기반 약물탐색**  
> 강좌: 도킹 프로그램 사용 실습 (홍승환, 한국제약바이오협회)  
> 실습 환경: **WSL2 Ubuntu** (CLI) + **UCSF Chimera Windows** (GUI)  
> 실습 파일: [`prac1.ipynb`](prac1.ipynb)

---

## 실습 목표

1. 도킹 전용 conda 환경 구성 및 필수 패키지 설치
2. 수용체(receptor)와 리간드(ligand)를 PDBQT 형식으로 변환
3. UCSF Chimera GUI로 AutoDock Vina 도킹 실행 및 결과 분석
4. Linux CLI로 Quick Vina 2 도킹 실행 및 결과 비교
5. PyMOL로 도킹 결과 시각화

---

## 분석 대상

| 항목 | 값 |
|------|-----|
| 타겟 단백질 | TGFR1_HUMAN (TGF-β receptor type 1 kinase) |
| 수용체 파일 | `3HMMA_receptor_HOH.pdb` |
| 리간드 | 3HMM 공결정 리간드 855 (`3HMMA_855.pdb`) |
| Docking Box | center=(17.994, 68.852, 7.446), size=(26.718, 19.984, 22.033) |

---

## 환경 구성

### WSL2 conda 환경 문제 및 해결

```bash
# 문제: WSL2 base 환경 (Python 3.7) + conda 23.1.0
# → vina, smina 설치 시 'Solving environment' 무한 대기
# → conda 26.x 업그레이드 불가 (Python 3.7 호환 불가)

# 해결: Python 3.10 기반 Miniconda 별도 설치
wget https://repo.anaconda.com/miniconda/Miniconda3-py310_23.11.0-2-Linux-x86_64.sh
bash Miniconda3-py310_23.11.0-2-Linux-x86_64.sh -b -p ~/miniconda310
source ~/miniconda310/etc/profile.d/conda.sh

# 도킹 전용 환경 생성
conda create -n docking python=3.10 -y
conda activate docking
conda install -c conda-forge vina smina openbabel -y
```

### 설치 패키지 목록

| 패키지 | 설치 방법 | 용도 |
|--------|-----------|------|
| AutoDock Vina 1.2.3 | `conda install -c conda-forge vina` | 표준 도킹 |
| Quick Vina 2 | 바이너리 직접 다운로드 | 고속 도킹 |
| Smina | `conda install -c conda-forge smina` | 커스텀 scoring |
| AutoDockTools (py3) | `git clone` + `setup.py install` | PDBQT 변환 |
| PBI Toolkit | `git clone` + `setup.py install` | 전처리 자동화 |
| numpy | `pip install numpy` | ADT 의존성 |
| rdkit | `pip install rdkit` | PBI 의존성 |
| openmm | `pip install openmm` | PBI 의존성 |
| pdbfixer | `pip install pdbfixer` | PBI 의존성 |

> **기존 miniconda3 (Python 3.7):** trDesign 실습용으로 그대로 유지 — 충돌 없음

---

## 리간드 전처리 워크플로우

### SMILES 변환 및 3D Conformer 생성

```bash
# SMILES 변환
obabel TGFR1/pdb/3HMMA_855.pdb -osmi
# → Cc1cccc(n1)c1nc(c2ccccc2n1)Nc1ccncc1

# pH 7.4 protonation state 확인
obabel -:"Cc1cccc(n1)c1nc(c2ccccc2n1)Nc1ccncc1" -p7.4 -osmi
# → 동일 SMILES (중성 분자)

# 3D conformer 생성 (obabel — 원자 이름 중복)
obabel -:"Cc1cccc(n1)c1nc(c2ccccc2n1)Nc1ccncc1" --gen3d -O M855_0.pdb

# 3D conformer 생성 (gen3d.py — 원자 이름 고유 번호)
gen3d.py -i "Cc1cccc(n1)c1nc(c2ccccc2n1)Nc1ccncc1" -o M855.pdb
```

### obabel vs gen3d.py 비교

| 항목 | M855_0.pdb (obabel) | M855.pdb (gen3d.py) |
|------|---------------------|---------------------|
| 원자 이름 | `C`, `N` (중복) | `C1`, `C2`, `N1` (고유) |
| 최적화 방법 | MMFF94 force field | RDKit 기반 |
| 용도 | 참고용 | 도킹 전처리 입력용 ✅ |

> **gen3d.py를 사용하는 이유:** AutoDock의 `prepare_ligand4`가 원자를 개별 식별하기 위해 고유한 원자 이름 필요.

### 전체 흐름

```
PDB 파일 (3D 결정 구조)
      ↓  obabel -osmi
SMILES (1D 문자열)
      ↓  obabel -p7.4
pH 7.4 protonation state 확인
      ↓  gen3d.py
3D Conformer (원자 이름 고유 번호)
      ↓  prepare_ligand4 / pdb2pdbqt.py
PDBQT (도킹 입력 파일)
```

---

## PDBQT 변환

```bash
# 방법 1: AutoDockTools
prepare_receptor4 -U nphs_lps -r 3HMMA_receptor_HOH.pdb -o 3HMMA_receptor_HOH.pdbqt
prepare_ligand4 -U nphs_lps -l 3HMMA_855.pdb -o 3HMMA_855.pdbqt

# 방법 2: PBI Toolkit
pdb2pdbqt.py -i 3HMMA_receptor_HOH.pdb -o 3HMMA_receptor_HOH.pdbqt -r
pdb2pdbqt.py -i 3HMMA_855.pdb -o 3HMMA_855.pdbqt -l
```

> **PDBQT:** AutoDock 계열 전용 파일 형식.  
> PDB에 **원자 타입(q)** 과 **부분 전하(t)** 정보가 추가된 형식.  
> 분자의 fragment를 강체로 취급하며 bond 정보가 없음.

---

## 실습 A — UCSF Chimera GUI 도킹

### 사전 준비

| 소프트웨어 | 다운로드 |
|-----------|---------|
| UCSF Chimera | https://www.cgl.ucsf.edu/chimera/download.html |
| AutoDock Vina 1.2.3 (Windows) | https://github.com/ccsb-scripps/AutoDock-Vina/releases/tag/v1.2.3 |

> ⚠️ **주의:** 파일 경로에 **한글 포함 금지**

### 실행 순서

```
File → Open → 3HMMA_receptor_HOH.pdb
File → Open → 3HMMA_855.pdb
Tools → Surface/Binding Analysis → AutoDock Vina
```

### 파라미터 설정

| 항목 | 값 |
|------|-----|
| Output file | `C:/Users/PC/docking_GUI/out.pdbqt` |
| Receptor | `3HMMA_receptor_HOH.pdb` |
| Ligand | `3HMMA_855.pdb` |
| Center X | 17.994 |
| Center Y | 68.852 |
| Center Z | 7.446 |
| Size X | 26.718 |
| Size Y | 19.984 |
| Size Z | 22.033 |
| Vina 실행파일 | `vina_1.2.3_windows_x86_64.exe` 전체 경로 |

### 마우스 조작법

| 동작 | 마우스 |
|------|--------|
| 회전 | 좌클릭 드래그 |
| 줌인/아웃 | 우클릭 드래그 |
| 수평 이동 | 휠 클릭 드래그 |

### 결과 파일 구조

| 파일 | 내용 |
|------|------|
| `out.conf` | Vina 파라미터 설정 파일 |
| `out.receptor.pdb` | 도킹용 수용체 원본 |
| `out.receptor.pdbqt` | 도킹용 수용체 (PDBQT) |
| `out.ligand.pdb` | 도킹용 리간드 초기 구조 |
| `out.ligand.pdbqt` | 도킹용 리간드 초기 구조 (PDBQT) |
| `out.pdbqt` | **도킹 결과** — 모든 결합 자세 + 스코어 |

### Chimera GUI 도킹 결과

| Rank | Score (kcal/mol) | RMSD lb (Å) | RMSD ub (Å) |
|------|-----------------|-------------|-------------|
| **1** | **-10.588** | 0.000 | 0.000 |
| 2 | -9.668 | 1.185 | 5.917 |
| 3 | -9.030 | 1.727 | 6.757 |
| 4 | -9.027 | 1.715 | 5.207 |
| 5 | -8.971 | 1.765 | 6.945 |

### Rank 1 에너지 분해

| 항목 | 값 (kcal/mol) | 설명 |
|------|--------------|------|
| **VINA RESULT** | **-10.588** | 최종 결합 친화도 |
| INTER + INTRA | -13.247 | 상호작용 + 내부 에너지 합 |
| INTER | -12.973 | 단백질-리간드 상호작용 에너지 |
| INTRA | -0.273 | 리간드 내부 에너지 |
| UNBOUND | -0.802 | 결합 전 리간드 자유 에너지 |

**활성 비틀림 (Active Torsions): 3개**

| # | 결합 원자 쌍 |
|---|-------------|
| 1 | C6 — C8 |
| 2 | C10 — N5 |
| 3 | C12 — N5 |

---

## 실습 B — Linux CLI 도킹 (Quick Vina 2)

### config.txt 설정

```bash
cat > config.txt << 'EOF'
receptor=3HMMA_receptor_HOH.pdbqt
center_x=17.994
center_y=68.852
center_z=7.446
size_x=26.718
size_y=19.984
size_z=22.033
cpu=10
num_modes=10
exhaustiveness=10
EOF
```

### 도킹 실행

```bash
# Quick Vina 2로 도킹
~/docking_tools/qvina02 --config config.txt \
    --ligand 3HMMA_855.pdbqt \
    --out dock_3HMMA_855.pdbqt

# pdbqt → pdb 변환 (bond 정보 복원)
~/miniconda310/envs/docking/bin/pdbqt2pdb_ref.py \
    -i dock_3HMMA_855.pdbqt \
    -o dock_3HMMA_855.pdb \
    -r 3HMMA_855.pdb
```

> **pdbqt2pdb_ref.py 실행 시 주의:**  
> `miniconda3` (구 환경)이 아닌 `miniconda310` 경로 명시 필요.  
> `obabel`을 `/usr/local/bin`에 심볼릭 링크 등록 필요:  
> `sudo ln -sf ~/miniconda310/envs/docking/bin/obabel /usr/local/bin/obabel`

### Quick Vina 2 도킹 결과

| Mode | Score (kcal/mol) | RMSD lb (Å) | RMSD ub (Å) |
|------|-----------------|-------------|-------------|
| **1** | **-11.4** | 0.000 | 0.000 |
| 2 | -10.3 | 1.089 | 5.849 |
| 3 | -10.2 | 1.695 | 2.119 |
| 4 | -9.5 | 1.795 | 6.835 |
| 5 | -9.3 | 1.779 | 6.877 |
| 6~10 | -9.1 ~ -8.8 | 2.3 ~ 3.2 | 6.3 ~ 8.0 |

- **계산 시간: 1.955초** (Quick Vina 2의 속도 장점)
- **Mode 1 RMSD = 0.000** → 기준 자세
- **Mode 2-3 RMSD lb < 2Å** → 결정 구조와 매우 유사

---

## Chimera GUI vs Quick Vina 2 비교

| 항목 | Chimera GUI (Vina 1.2.3) | Quick Vina 2 |
|------|--------------------------|--------------|
| Rank 1 Score | -10.588 kcal/mol | **-11.4 kcal/mol** |
| 계산 속도 | 상대적으로 느림 | **1.955초** |
| 사용 편의성 | GUI (직관적) | CLI (자동화 가능) |
| 적합한 용도 | 단일 구조 확인 | 대규모 가상 탐색 |

---

## PyMOL 시각화

```python
# PyMOL 콘솔에서
cd C:/Users/PC/.../3_1_docking_setup

load 3HMMA_receptor_HOH.pdb
load 3HMMA_855.pdb
load dock_3HMMA_855.pdb

# 스타일 설정
hide everything
show cartoon, 3HMMA_receptor_HOH
show sticks, 3HMMA_855
show sticks, dock_3HMMA_855

# 색상 구분
color cyan, 3HMMA_receptor_HOH     # 단백질
color yellow, 3HMMA_855             # 결정 구조 리간드
color magenta, dock_3HMMA_855       # 도킹 결과 리간드

zoom
```

> 결정 구조 리간드(노랑)와 도킹 결과(마젠타)가 겹쳐 보이면 재현 성공 ✅

---

## 결과 파일 구조

```
3_1_docking_setup/
├── prac1.py                              ← 전체 실습 코드
├── 3HMMA_receptor_HOH.pdb               ← 수용체 원본
├── 3HMMA_receptor_HOH.pdbqt             ← 수용체 PDBQT
├── 3HMMA_855.pdb                        ← 리간드 원본
├── 3HMMA_855.pdbqt                      ← 리간드 PDBQT
├── config.txt                           ← 도킹 박스 설정
├── dock_3HMMA_855.pdbqt                 ← 도킹 결과 (PDBQT)
├── dock_3HMMA_855.pdb                   ← 도킹 결과 (PDB)
├── docking_practice_with_UCSF_chimera/  ← Chimera GUI 결과
├── docking_practice_with_linux_and_pymol/ ← CLI + PyMOL 결과
├── TGFR1/                               ← PBI toolkit 전처리 파일
└── TGFR1_test/                          ← 테스트용 파일
```

---

## 명시적 한계

| 한계 | 설명 |
|------|------|
| 강체 도킹 | 수용체 유연성 미반영 (Induced fit 무시) |
| 단일 리간드 | 1개 리간드만 실습 — 가상 탐색 미수행 |
| scoring function | Vina scoring은 결합 자유에너지 근사값 |
| 물 분자 | 일부 bridging water 포함했으나 완전하지 않음 |
| WSL 환경 | miniconda3/miniconda310 이중 환경으로 PATH 충돌 주의 |



---

## 참고문헌

1. Trott O, Olson AJ. **AutoDock Vina: improving the speed and accuracy of docking.** *J Comput Chem* 31(2), 455–461 (2010).
2. Alhossary A, et al. **Fast, Accurate, and Reliable Molecular Docking with QuickVina 2.** *Bioinformatics* 31(13), 2214–2216 (2015).
3. Eberhardt J, et al. **AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings.** *J Chem Inf Model* 61(8), 3891–3898 (2021).

---

## 관련 링크

| 리소스 | URL |
|--------|-----|
| AutoDock Vina | https://github.com/ccsb-scripps/AutoDock-Vina |
| Quick Vina 2 | https://github.com/QVina/qvina |
| AutoDockTools py3 | https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3 |
| UCSF Chimera | https://www.cgl.ucsf.edu/chimera/download.html |
| PBI Toolkit | https://github.com/gicsaw/PBI |
