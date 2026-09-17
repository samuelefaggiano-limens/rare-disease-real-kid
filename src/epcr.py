"""EPCR scoring for Track 1, Step 5.

``epcr`` (Estimated Probability of Causal Relationship) is the hackathon's
custom confidence score in (0, 1]. The leaderboard computes rank points and
F-max (see the Space's ``evaluation.py``); F-max sweeps every distinct epcr
value as a threshold and keeps the best F1, counting each *variant* as one
prediction. Consequences that shape this module:

  * Separation matters more than absolute calibration. A single confident row
    followed by a clearly lower tail can reach F-max = 1.0; ten near-tied rows
    cap precision at 1/10 even when the answer is among them.
  * A compound-het pair must be one row (two variants). If only one of the
    two true variants is in a row, the row earns half rank points.

Every component below is a number in [0, 1]. They are combined with a
weighted *geometric* mean so weak evidence on any axis pulls the score down
multiplicatively, which is what creates separation between a candidate that
is strong everywhere and the many that are strong on one axis only.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

# Default weights; sum to 1. Phenotype fit and pathogenicity dominate because
# after Step 2/3 every candidate is already rare and protein-altering, so the
# frequency and impact axes have little remaining spread.
DEFAULT_WEIGHTS = {
    "phenotype": 0.35,
    "pathogenicity": 0.30,
    "inheritance": 0.20,
    "frequency": 0.15,
}
COMPONENT_FLOOR = 0.02  # keeps ln() finite; a zero axis still costs ~x0.25 at w=0.35

IMPACT_BASE = {"HIGH": 0.90, "MODERATE": 0.45, "LOW": 0.15, "MODIFIER": 0.05}


def pathogenicity_score(impact: str | None, cadd_phred: float | None,
                        revel: float | None, clinvar_sig: str | None,
                        avi_phred: float | None = None) -> float:
    """Combine snpEff impact with CADD/REVEL/AVI/ClinVar into [0, 1].

    HIGH impact (stop-gain, frameshift, canonical splice) is close to
    maximal on its own. MODERATE (missense, in-frame indel) starts at 0.45
    and is lifted by in-silico scores when available: CADD_PHRED >= 20 is the
    conventional "top 1%" line, REVEL >= 0.5 the common pathogenic cut, and
    AlphaGenome AVI_PHRED shares CADD's Phred scale (>= 20 top 1 %, >= 30 top
    0.1 % genome-wide). A ClinVar Pathogenic/Likely pathogenic assertion
    overrides to >= 0.95; Benign/Likely benign caps at 0.10.
    """
    base = IMPACT_BASE.get((impact or "").upper(), 0.30)
    evidence = []
    if cadd_phred is not None:
        evidence.append(min(max(cadd_phred, 0.0) / 30.0, 1.0))   # 30+ ~ top 0.1 %
    if avi_phred is not None:
        evidence.append(min(max(avi_phred, 0.0) / 30.0, 1.0))    # same Phred scale
    if revel is not None:
        evidence.append(min(max(revel, 0.0), 1.0))
    score = base
    if evidence:
        ins = max(evidence)
        # in-silico evidence can lift MODERATE towards HIGH, or temper it.
        score = 0.5 * base + 0.5 * max(ins, 0.15 if base >= 0.45 else 0.0)
        if base >= 0.90:                        # HIGH stays high unless scores are poor
            score = max(score, 0.75)
    sig = (clinvar_sig or "").lower()
    if "pathogenic" in sig and "benign" not in sig and "conflicting" not in sig:
        score = max(score, 0.95)
    elif "benign" in sig and "pathogenic" not in sig and "conflicting" not in sig:
        score = min(score, 0.10)
    return round(min(max(score, 0.0), 1.0), 4)


def frequency_score(gnomad_af: float | None) -> float:
    """Rarer is better. Missing gnomAD AF means absent from gnomAD => 1.0.

    Step 2 already removed dbSNP Common (>= ~1 %), so anything left is rare;
    this axis only orders "never seen" ahead of "seen a few times". Kept
    linear in -log10(AF) between 1e-2 (0.0) and 1e-5 (1.0). Population-
    specific caveat (docs/00-references.md section 4): a founder variant can
    reach several percent in one subpopulation, which is why this is one of
    the lighter weights and not a hard filter.
    """
    if gnomad_af is None or gnomad_af <= 0:
        return 1.0
    lo, hi = 2.0, 5.0  # -log10 bounds
    v = (-math.log10(gnomad_af) - lo) / (hi - lo)
    return round(min(max(v, 0.0), 1.0), 4)


def inheritance_score(zygosity: str, n_het_in_gene: int, mva_recessive: bool = True) -> float:
    """How well the genotype pattern fits the expected model.

    MVA is autosomal recessive in every described gene, so biallelic
    evidence (homozygous alt, or >= 2 heterozygous variants in the gene that
    can be paired as compound het) scores 1.0. A lone heterozygous variant
    still fits a de novo dominant or a het-plus-undetected-CNV scenario, so
    it is penalised, not excluded.
    """
    z = zygosity.upper()
    if z == "HOM_ALT":
        return 1.0
    if z == "HET":
        if n_het_in_gene >= 2:
            return 1.0
        return 0.45 if mva_recessive else 0.85
    return 0.2


def combine(components: dict[str, float], weights: dict[str, float] | None = None,
            bonus: float = 1.0) -> float:
    """Weighted geometric mean of components, times a bonus, clipped to (0, 1]."""
    w = weights or DEFAULT_WEIGHTS
    total_w = sum(w.values())
    log_sum = 0.0
    for k, wk in w.items():
        c = max(components.get(k, 0.0), COMPONENT_FLOOR)
        log_sum += wk * math.log(c)
    score = math.exp(log_sum / total_w) * bonus
    return round(min(max(score, 1e-4), 1.0), 4)


@dataclass
class CandidateRow:
    """One submission row: 1 variant, or 2 for a compound-het pair."""
    variants: tuple[tuple[str, int, str, str], ...]   # (chrom, pos, ref, alt)
    raw_score: float
    gene: str
    notes: str = ""
    finding_type: str = "primary"
    epcr: float = field(default=0.0)


def assign_epcr(rows: list[CandidateRow], max_rows: int = 5, tail_ratio: float = 0.6,
                top_cap: float = 0.95) -> list[CandidateRow]:
    """Turn raw scores into a ranked, separated epcr list.

    * Sort by raw score, keep at most ``max_rows`` (leaderboard max is 10;
      fewer rows means fewer false positives at every threshold).
    * Row 1 gets ``min(raw, top_cap)``.
    * Each following row is capped at ``tail_ratio`` x the previous epcr, so
      the list always decays geometrically even when raw scores are near-tied.
      With tail_ratio 0.6 the 5th row is <= 13 % of the first.

    Ties in raw scores are what F-max punishes; a forced decay converts
    "these two look equally good" into an explicit ranking decision, which
    is what the metric rewards. If two candidates are genuinely equal the
    order is decided by the stable sort on the input list.
    """
    ranked = sorted(rows, key=lambda r: -r.raw_score)[:max_rows]
    prev = None
    for r in ranked:
        e = min(r.raw_score, top_cap)
        if prev is not None:
            e = min(e, prev * tail_ratio)
        e = max(round(e, 4), 1e-4)
        r.epcr = e
        prev = e
    return ranked


# ---- Local re-implementation of the leaderboard metrics (for what-if reports) ----

RANK_POINT_TIERS = [(1, 100), (3, 50), (5, 25), (10, 10)]


def _rank_points(rank: int) -> int:
    for max_rank, pts in RANK_POINT_TIERS:
        if rank <= max_rank:
            return pts
    return 0


def score_against(rows: list[CandidateRow], true_variants: frozenset) -> dict:
    """Rank points and F-max exactly as the Space's evaluation.py computes them."""
    ordered = sorted(rows, key=lambda r: -r.epcr)
    full = partial = None
    for i, r in enumerate(ordered, 1):
        if frozenset(r.variants) == true_variants:
            full = i
            break
    if len(true_variants) == 2 and full is None:
        for i, r in enumerate(ordered, 1):
            if frozenset(r.variants) & true_variants:
                partial = i
                break
    if full is not None:
        rp = float(_rank_points(full))
    elif partial is not None:
        rp = 0.5 * _rank_points(partial)
    else:
        rp = 0.0

    best_f, best_t = 0.0, None
    for t in sorted({r.epcr for r in ordered}, reverse=True):
        pred: set = set()
        for r in ordered:
            if r.epcr >= t:
                pred |= set(r.variants)
        tp = len(pred & true_variants)
        fp = len(pred - true_variants)
        fn = len(true_variants - pred)
        p = tp / (tp + fp) if tp + fp else 0.0
        rcl = tp / (tp + fn) if tp + fn else 0.0
        f = 2 * p * rcl / (p + rcl) if p + rcl else 0.0
        if f > best_f:
            best_f, best_t = f, t
    return {"rank_points": rp, "f_max": round(best_f, 3), "f_max_threshold": best_t,
            "full_match_rank": full, "partial_match_rank": partial}


def what_if_table(rows: list[CandidateRow]) -> list[dict]:
    """For each row k: the leaderboard result if row k were the true answer."""
    out = []
    for k, r in enumerate(sorted(rows, key=lambda x: -x.epcr), 1):
        res = score_against(rows, frozenset(r.variants))
        out.append({"if_true_row": k, "epcr": r.epcr, **res})
    return out
