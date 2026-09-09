# LAIDD Drug Discovery Portfolio

> **LAIDD (Lectures on AI-driven Drug Discovery)**  
> 한국제약바이오협회 AI 신약융합연구원 온라인 교육 실습 포트폴리오  
> https://www.laidd.org

---

## 개요

본 리포지토리는 LAIDD에서 제공하는 AI 기반 신약개발 강의 중 **실습이 포함된 과정**의 코드, 결과, 분석을 체계적으로 정리한 포트폴리오입니다.

강의는 두 개의 대주제로 구성됩니다:

| 대주제 | 폴더 | 핵심 방법론 |
|--------|------|------------|
| 단백질 구조 기반 약물탐색 | `structure-based/` | 구조 예측, 도킹, 단백질 디자인, MD 시뮬레이션 |
| 리간드 기반 약물탐색 | `ligand-based/` | QSAR, GNN, 분자 생성 모델, 독성 예측 |

---

## 리포지토리 구조

```
LAIDD-drug-discovery/
├── README.md                          
│
├── structure-based/
│   └── trdesign-hallucination/        
│       ├── README.md
│       ├── analyze_trdesign.py
│       ├── output_L50.fa
│       ├── analysis_result.png
│       └── esm_structure.png
│
└── ligand-based/
    └── (실습 추가 예정)
```

---

## 단백질 구조 기반 약물탐색 (Structure-Based)

**강좌명:** 신약개발을 위한 단백질 구조 예측 및 상호작용 예측  
**교수자:** 석차옥 (서울대학교)  
**총 강의:** 11강

| # | 강의 제목 | 핵심 내용 | 실습 | 상태 |
|---|-----------|-----------|------|------|
| 1 | 단백질 구조예측 및 상호작용 예측 개괄 | 구조예측·도킹·디자인 발전사, 신약개발 활용 분야 | trDesign hallucination | ✅ [`trdesign-hallucination/`]([structure-based/trdesign-hallucination/](https://github.com/OHfromH2O/LAIDD_drug_discovery/tree/main/structure_based_drug_discovery/trdesign_hallucination)) |
| 2 | 단백질 구조의 기초 | 4단계 계층 구조, 구조 뷰어(PyMOL/Chimera), SCOP/CATH 분류 | 구조 뷰어 실습 | — |
| 3 | 단백질 접힘 원리 | 소수성 효과, 수소결합, 샤페론 매개 접힘 | — | — |
| 4 | 호몰로지 모델링 | 주형 기반 구조 예측, MODELLER, Swiss-Model | 웹서버 실습 | — |
| 5 | AI 기반 구조 예측 | AlphaFold/trRosetta (2D CNN), AlphaFold2/RoseTTAFold (Attention) | 웹서버 실습 | — |
| 6 | 단백질-단백질 도킹 (원리 기반) | PPI 입문, Rosetta, GalaxyTongDock | GalaxyTongDock 실습 | — |
| 7 | 단백질-단백질 도킹 (정보 기반) | GalaxyHomomer/Heteromer, CASP/CAPRI 사례 | 실습 | — |
| 8 | 분자동력학 시뮬레이션 | MD 원리, 신약개발 적용 사례 | — | — |
| 9 | Rosetta 기반 단백질 디자인 | 서열·뼈대 디자인, 항체/항원 디자인 | Rosetta 실습 | — |
| 10 | AI 기반 단백질 디자인 | 생물정보학 활용 디자인, trDesign 논문 리뷰 | trDesign 실습 | — |
| 11 | 단백질-단백질 결합 친화도 | 기계학습+실험 결합, PRODIGY/SSIPe 실습 | 웹서버 실습 | — |

**추가 강좌 (Structure-Based 확장):**

| 강좌명 | 교수자 | 실습 |
|--------|--------|------|
| 도킹 프로그램 사용 실습 | 홍승환 (한국제약바이오협회) | AutoDock Vina, rDock, Python 스크립트 |
| AI in Predicting Protein-Ligand Interaction (structure-based) | 김우연 (KAIST) | 이론 강의 (3D CNN, GNN 기반 모델) |
| 구조 기반 가상 탐색을 활용한 유효물질 발굴과 최적화 | 이세한 (㈜히츠) | 가상 탐색, Hit-to-Lead |
| 딥러닝을 이용한 단백질 도킹 | 이유한 (카카오브레인) | SE(3)-equivariant, End-to-end docking |
| Molecular design with deep generative models | 임재창 (HITS) | SMILES/그래프 기반 생성모델 |
| 분자생성모델 연구동향 리뷰 | 임재창 (HITS) | 이론 강의 |
| AI 기반 protein-ligand interaction 예측 최신동향 | 황상연 (HITS) | 이론 강의 |
| 신약후보물질 탐색 및 최적화를 위한 딥러닝 모델 | 김동섭 (KAIST) | GCN, 강화학습, 도킹 |
| 단백질-리간드 상호작용 계산을 위한 분자동역학 시뮬레이션 | 최정모 (부산대학교) | MD 시뮬레이션, 자유에너지 계산 |
| 생물정보학을 활용한 단백질 간 상호작용 및 복합체 모델링 | 이윤지 (중앙대학교) | 이론 강의 |
| 단백질 언어 모델을 활용한 컨텍트 예측 | 김재훈 (카카오브레인) | HuggingFace, ESM, Contact prediction |
| 면역정보학과 단백질 재설계 | 최윤주 (전남대학교) | 항체 재설계, MHC-펩타이드 결합 예측 |

---

## 리간드 기반 약물탐색 (Ligand-Based)

| # | 강좌명 | 교수자 | 핵심 내용 | 실습 | 상태 |
|---|--------|--------|-----------|------|------|
| 1 | QSAR | 김동섭 (KAIST) | Descriptor, ML 기반 활성 예측, Proteochemometrics | Python | — |
| 2 | 독성예측 인공지능 모델 활용 | 신현길 (안전성평가연구소) | pandas, pubchempy, MOPAC, scikit-learn | Python | — |
| 3 | AI in Predicting Drug-protein Interaction (sequence-based) | 남호정 (GIST) | 서열 기반 화합물-단백질 상호작용 예측 | 이론 강의 | — |
| 4 | Graph Neural Networks for Molecular Property Prediction | 류성옥 (Galux) | GCN, GIN, GAT, GGNN, Bayesian Learning | Python | — |
| 5 | 그래프 트랜스포머를 활용한 분자물성 예측 | 이유한 (카카오브레인) | Attention, Graph Transformer | Python | — |
| 6 | Deep Learning Based Molecular Generation | 이일구 (팜캐드) | RNN, ChemicalVAE, de novo 분자 생성 | Python (PyTorch) | — |
| 7 | Disease-Target-Drug relationship analysis from multi-dimensional data | 김현욱 (KAIST) | 약물반응 예측, 다차원 데이터 분석 | 이론 강의 | — |

---

## 실습 완료 목록

| 폴더 | 강의 | 주요 결과 | 사용 도구 |
|------|------|-----------|-----------|
| [`structure-based/trdesign-hallucination/`](structure-based/trdesign-hallucination/) | 1강 — 단백질 구조예측 개괄 | 50잔기 de novo 서열 디자인, ESMFold 구조 예측, 접힘 원리 분석 | trDesign, TF 1.14, ESMFold |

---

## 환경 요약

| 구분 | 도구 | 버전 |
|------|------|------|
| OS | Windows 11 + WSL2 (Ubuntu 22.04) | — |
| 패키지 관리 | Miniconda | 23.1.0 |
| Python (structure-based) | conda env `trdesign` | 3.7 |
| Python (ligand-based) | conda env (예정) | 3.9+ |
| 딥러닝 프레임워크 | TensorFlow (structure) / PyTorch (ligand) | 1.14 / 최신 |

---

## 참고문헌 (전체 과정)

- Anishchenko et al. **De novo protein design by deep network hallucination.** *Nature* 600 (2021)
- Yang et al. **Improved protein structure prediction using predicted interresidue orientations.** *PNAS* 117 (2020)
- Jumper et al. **Highly accurate protein structure prediction with AlphaFold.** *Nature* 596 (2021)
- Lin et al. **Evolutionary-scale prediction of atomic-level protein structure with a language model.** *Science* 379 (2023)
- De Vivo et al. **Role of Molecular Dynamics and Related Methods in Drug Discovery.** *J. Med. Chem.* (2016)
- Norn et al. **Protein sequence design by explicit energy landscape optimization.** *PNAS* (2021)

---

## 관련 링크

| 리소스 | URL |
|--------|-----|
| LAIDD 공식 사이트 | https://www.laidd.org |
| Protein Data Bank | https://www.rcsb.org |
| ESMFold | https://esmatlas.com/resources?action=fold |
| AlphaFold2 Colab | https://colab.research.google.com/github/sokrypton/ColabFold |
| trDesign GitHub | https://github.com/gjoni/trDesign |
| PRODIGY (결합 친화도) | https://wenmr.science.uu.nl/prodigy/ |
