"""
trDesign 결과 분석 스크립트
=============================
목표: 디자인된 서열이 단백질 접힘 원리에 얼마나 부합하는지 정량적으로 평가

분석 항목:
  1. MCMC trajectory score 수렴 곡선
  2. 아미노산 조성 vs 자연 단백질 평균 (UniProtKB/Swiss-Prot)
  3. Kyte-Doolittle 소수성 프로파일 (window=9)
  4. 잔기 분류 및 등전점(pI) 추정
  5. 단백질 접힘 원리 부합성 종합 평가

가정 및 한계:
  - 자연 단백질 평균 조성: Doolittle (1989), UniProtKB 기준
  - KD 소수성: 1차원 서열 기반 근사 (3D 구조 미반영)
  - pI 추정: Henderson-Hasselbalch 단순화 모델 (이온 강도 무시)
  - 2차 구조 예측은 외부 도구(PSIPRED, ESMFold) 필요 → 미포함

사용법:
  python analyze_trdesign.py

출력:
  - 터미널: 정량 분석 보고서
  - analysis_result.png: 3패널 시각화 그래프

Author: Jay (실습)
Date: 2026-09-09
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')   # GUI 없는 WSL 환경용 백엔드
import matplotlib.pyplot as plt
import csv
from collections import Counter


# ── 상수 정의 ─────────────────────────────────────────────────────────────────

FINAL_SEQ = "LNYSHLKKIARDYARHYGLQDVEWESTPEDGIFTVKGRFNGRDAQVMFPI"
CSV_FILE  = "output_L50.csv"

# 자연 단백질 평균 아미노산 조성 (%, UniProtKB/Swiss-Prot)
NATURAL_COMPOSITION = {
    'A': 8.25, 'R': 5.53, 'N': 4.06, 'D': 5.45, 'C': 1.37,
    'E': 6.32, 'Q': 3.93, 'G': 7.07, 'H': 2.27, 'I': 5.96,
    'L': 9.66, 'K': 5.84, 'M': 2.42, 'F': 3.86, 'P': 4.70,
    'S': 6.56, 'T': 5.34, 'W': 1.08, 'Y': 2.92, 'V': 6.87
}

# Kyte-Doolittle 소수성 척도 (Kyte & Doolittle, 1982)
KD_SCALE = {
    'A': 1.8,  'R':-4.5, 'N':-3.5, 'D':-3.5, 'C': 2.5,
    'E':-3.5,  'Q':-3.5, 'G':-0.4, 'H':-3.2, 'I': 4.5,
    'L': 3.8,  'K':-3.9, 'M': 1.9, 'F': 2.8, 'P':-1.6,
    'S':-0.8,  'T':-0.7, 'W':-0.9, 'Y':-1.3, 'V': 4.2
}

# 단순화 pKa (Henderson-Hasselbalch 근사용)
PKA = {
    'D': 3.9, 'E': 4.1, 'H': 6.0,
    'C': 8.3, 'Y':10.1, 'K':10.5,
    'R':12.5, 'N_term': 8.0, 'C_term': 3.1
}

# 잔기 분류 집합
HYDROPHOBIC = set('AILMFWYVP')
CHARGED_POS = set('KRH')
CHARGED_NEG = set('DE')
POLAR       = set('STNQ')


# ── 함수 정의 ─────────────────────────────────────────────────────────────────

def load_trajectory(csv_file: str):
    """
    CSV trajectory 파일 로드.

    Parameters
    ----------
    csv_file : str
        hallucinate.py --ocsv 옵션으로 생성된 CSV 경로

    Returns
    -------
    steps : list[int]
    seqs  : list[str]
    scores: list[float]
    """
    steps, seqs, scores = [], [], []
    with open(csv_file, 'r') as fh:
        for row in csv.DictReader(fh):
            steps.append(int(row['step']))
            seqs.append(row['sequence'])
            scores.append(float(row['score']))
    return steps, seqs, scores


def aa_composition(seq: str) -> dict:
    """
    아미노산 조성 계산 (%).

    가정: 20종 표준 아미노산만 존재.
    비표준 문자는 Counter에 포함되나 NATURAL_COMPOSITION 키에만 출력.
    """
    c = Counter(seq)
    t = len(seq)
    return {aa: round(c.get(aa, 0) / t * 100, 2) for aa in NATURAL_COMPOSITION}


def kd_profile(seq: str, window: int = 9) -> list:
    """
    Kyte-Doolittle sliding window 소수성 프로파일.

    Parameters
    ----------
    seq    : str   아미노산 서열
    window : int   윈도우 크기 (문헌 표준: 9)

    한계: 서열 경계부는 실제 window보다 작은 구간을 평균하므로
          경계 효과(boundary effect) 존재.
    """
    sc   = [KD_SCALE.get(aa, 0.0) for aa in seq]
    half = window // 2
    return [
        np.mean(sc[max(0, i - half):min(len(seq), i + half + 1)])
        for i in range(len(seq))
    ]


def estimate_pI(seq: str) -> float:
    """
    Henderson-Hasselbalch 근사로 등전점(pI) 추정.

    가정:
      - 잔기 pKa는 독립적 (인접 잔기 영향 무시)
      - 이온 강도 0 가정
      - 비표준 pKa 사용 (실제 단백질 pI와 ±0.5 내외 오차 가능)
    """
    def charge_at_pH(pH: float) -> float:
        q  = 1.0 / (1.0 + 10 ** (pH - PKA['N_term']))
        q -= 1.0 / (1.0 + 10 ** (PKA['C_term'] - pH))
        handlers = {
            'D': lambda p: -1.0 / (1.0 + 10 ** (PKA['D'] - p)),
            'E': lambda p: -1.0 / (1.0 + 10 ** (PKA['E'] - p)),
            'H': lambda p:  1.0 / (1.0 + 10 ** (p - PKA['H'])),
            'C': lambda p: -1.0 / (1.0 + 10 ** (PKA['C'] - p)),
            'Y': lambda p: -1.0 / (1.0 + 10 ** (PKA['Y'] - p)),
            'K': lambda p:  1.0 / (1.0 + 10 ** (p - PKA['K'])),
            'R': lambda p:  1.0 / (1.0 + 10 ** (p - PKA['R'])),
        }
        for aa in seq:
            if aa in handlers:
                q += handlers[aa](pH)
        return q

    lo, hi = 0.0, 14.0
    for _ in range(100):          # 이분법 100회 → 수렴 보장
        mid = (lo + hi) / 2.0
        if charge_at_pH(mid) > 0:
            lo = mid
        else:
            hi = mid
    return round((lo + hi) / 2.0, 2)


def evaluate_folding(seq: str, kd: list) -> tuple:
    """
    단백질 접힘 원리 부합성 규칙 기반 평가.

    Returns
    -------
    passes : list[str]   통과 항목
    issues : list[str]   주의 항목
    """
    n_pos   = sum(seq.count(a) for a in 'KRH')
    n_neg   = sum(seq.count(a) for a in 'DE')
    n_core  = sum(1 for v in kd if v > 1.6)
    n_pro   = seq.count('P')
    n_trp   = seq.count('W')

    passes, issues = [], []

    # 소수성 core 후보 잔기
    if n_core >= 8:
        passes.append(f"소수성 core 후보 잔기 충분 ({n_core}개, KD>1.6)")
    else:
        issues.append(f"소수성 core 후보 잔기 부족 ({n_core}개 < 8)")

    # Cys
    if 'C' not in seq:
        passes.append("Cys 없음 → 비정상 disulfide 위험 없음 (--rm_aa=C 효과)")
    else:
        issues.append("Cys 존재 → disulfide 환경 확인 필요")

    # Pro 과다
    if n_pro > 5:
        issues.append(f"Pro {n_pro}개 → 2차 구조(helix/sheet) 형성 방해 가능")
    else:
        passes.append(f"Pro {n_pro}개 → 적정 수준")

    # Trp 과다
    if n_trp > 3:
        issues.append(f"Trp {n_trp}개 → 입체 충돌 위험")
    else:
        passes.append(f"Trp {n_trp}개 → 적정 수준")

    # 전하 균형
    if abs(n_pos - n_neg) > 8:
        issues.append(f"전하 불균형 심함 (양전하 {n_pos}개 / 음전하 {n_neg}개)")
    else:
        passes.append(f"전하 균형 양호 (양전하 {n_pos}개 / 음전하 {n_neg}개)")

    # 서열 길이
    if len(seq) < 40:
        issues.append(f"서열 길이 {len(seq)}잔기 → 독립 구조 형성 어려울 수 있음")
    else:
        passes.append(f"서열 길이 {len(seq)}잔기 → 구조 형성 가능 범위")

    return passes, issues


# ── 메인 실행 ─────────────────────────────────────────────────────────────────

def main():
    seq = FINAL_SEQ
    print("=" * 65)
    print("  trDesign 결과 분석 보고서")
    print("=" * 65)
    print(f"\n  분석 서열 (L={len(seq)}):")
    print(f"  {seq}\n")

    # 데이터 로드
    steps, seqs, scores = load_trajectory(CSV_FILE)
    comp = aa_composition(seq)
    kd   = kd_profile(seq, window=9)
    pI   = estimate_pI(seq)
    n_pos = sum(seq.count(a) for a in 'KRH')
    n_neg = sum(seq.count(a) for a in 'DE')

    # ── 1. 아미노산 조성 ──────────────────────────────────────────────────────
    print("[ 1. 아미노산 조성 분석 ]")
    print(f"  {'AA':<4} {'Designed(%)':>12} {'Natural(%)':>12} {'Diff':>8}  판정")
    print(f"  {'-' * 52}")
    for aa in sorted(NATURAL_COMPOSITION):
        d, n = comp[aa], NATURAL_COMPOSITION[aa]
        flag = "◀ 과다" if d - n > 3 else ("◀ 부족" if n - d > 3 else "")
        print(f"  {aa:<4} {d:>12.1f} {n:>12.1f} {d-n:>+8.1f}  {flag}")

    # ── 2. 잔기 분류 ──────────────────────────────────────────────────────────
    cls = {
        '소수성 (AILMFWYVP)': sum(a in HYDROPHOBIC for a in seq),
        '양전하 (KRH)':       sum(a in CHARGED_POS for a in seq),
        '음전하 (DE)':        sum(a in CHARGED_NEG for a in seq),
        '극성   (STNQ)':      sum(a in POLAR       for a in seq),
    }
    print(f"\n[ 2. 잔기 분류 (총 {len(seq)}잔기) ]")
    for k, v in cls.items():
        bar = '█' * v
        print(f"  {k:<22}: {v:>2}개 ({v/len(seq)*100:>4.1f}%) {bar}")

    # ── 3. 전하 / pI ──────────────────────────────────────────────────────────
    print(f"\n[ 3. 전하 분석 ]")
    print(f"  양전하 잔기 (K+R+H) : {n_pos}개")
    print(f"  음전하 잔기 (D+E)   : {n_neg}개")
    print(f"  추정 등전점 (pI)    : {pI}")
    if pI < 5.0:
        note = "산성 단백질 (생리적 pH 7.4에서 음전하 우세)"
    elif pI > 9.0:
        note = "염기성 단백질 (생리적 pH 7.4에서 양전하 우세)"
    else:
        note = "중성 범위"
    print(f"  → {note}")

    # ── 4. 소수성 프로파일 ────────────────────────────────────────────────────
    n_core = sum(1 for v in kd if v > 1.6)
    print(f"\n[ 4. Kyte-Doolittle 소수성 프로파일 (window=9) ]")
    print(f"  평균 KD          : {np.mean(kd):+.3f}")
    print(f"  최대 KD          : {np.max(kd):+.3f}")
    print(f"  최소 KD          : {np.min(kd):+.3f}")
    print(f"  Core 잔기 (>1.6) : {n_core}개 / {len(seq)}개")

    # ── 5. 종합 평가 ──────────────────────────────────────────────────────────
    passes, issues = evaluate_folding(seq, kd)
    print(f"\n[ 5. 단백질 접힘 원리 부합성 종합 평가 ]")
    print("\n  ✅ 통과 항목:")
    for p in passes:
        print(f"     · {p}")
    print("\n  ⚠️  주의 항목:")
    for i in issues:
        print(f"     · {i}")
    print(f"\n  ※ 한계: 1차 서열 기반 근사.")
    print(f"         실제 접힘 검증은 ESMFold/AlphaFold2 구조 예측 필요.")

    # ── 6. 시각화 ─────────────────────────────────────────────────────────────
    fig, axes = plt.subplots(3, 1, figsize=(13, 15))
    fig.suptitle(
        f'trDesign Sequence Analysis\n{seq}',
        fontsize=11, fontweight='bold'
    )

    # 패널 1: Score 수렴 곡선
    ax = axes[0]
    ax.plot(steps, scores, color='steelblue', lw=1.3, zorder=2)
    ax.axhline(scores[-1], color='red', ls='--', lw=1.0, alpha=0.8,
               label=f'Final score = {scores[-1]:.3f}')
    ax.fill_between(steps, scores, scores[-1], alpha=0.08, color='steelblue')
    ax.set_xlabel('MCMC Step', fontsize=10)
    ax.set_ylabel('Score (Loss)', fontsize=10)
    ax.set_title('Score Convergence (MCMC Trajectory)', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # 패널 2: AA 조성 비교
    ax = axes[1]
    aa_list  = sorted(NATURAL_COMPOSITION.keys())
    designed = [comp[aa] for aa in aa_list]
    natural  = [NATURAL_COMPOSITION[aa] for aa in aa_list]
    x, w     = np.arange(len(aa_list)), 0.35
    ax.bar(x - w/2, designed, w, label='Designed', color='steelblue', alpha=0.85)
    ax.bar(x + w/2, natural,  w, label='Natural avg (Swiss-Prot)',
           color='coral', alpha=0.85)
    ax.set_xticks(x)
    ax.set_xticklabels(aa_list, fontsize=9)
    ax.set_ylabel('Frequency (%)', fontsize=10)
    ax.set_title('Amino Acid Composition vs Natural Protein Average', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')

    # 패널 3: KD 소수성 프로파일
    ax = axes[2]
    pos    = list(range(1, len(seq) + 1))
    colors = ['steelblue' if v >= 0 else 'coral' for v in kd]
    ax.bar(pos, kd, color=colors, alpha=0.85)
    ax.axhline(0,   color='black', lw=0.8)
    ax.axhline(1.6, color='green', lw=1.2, ls='--',
               label='Hydrophobic core threshold (1.6)')
    for i, aa in enumerate(seq):
        ax.text(i + 1, min(kd) - 0.22, aa, ha='center', va='top',
                fontsize=6.5, color='black')
    ax.set_xlabel('Residue Position', fontsize=10)
    ax.set_ylabel('KD Hydrophobicity', fontsize=10)
    ax.set_title('Kyte-Doolittle Hydrophobicity Profile (window=9)', fontsize=11)
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_xticks(pos[::5])

    plt.tight_layout()
    out = 'analysis_result.png'
    plt.savefig(out, dpi=150, bbox_inches='tight')
    print(f"\n  그래프 저장 완료: {out}")

    print("\n" + "=" * 65)
    print("  다음 단계: ESMFold 구조 예측")
    print("  URL: https://esmatlas.com/resources?action=fold")
    print(f"  서열: {seq}")
    print("=" * 65)


if __name__ == '__main__':
    main()
