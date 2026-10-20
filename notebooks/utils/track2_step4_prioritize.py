"""Track 2, Step 4: integrate evidence and prioritise the final hypotheses.

This is a research prioritisation of repurposing *hypotheses* for follow-up
experiments, as specified in `docs/01-procedure.md` Step 4. It is explicitly
non-clinical: nothing here is a treatment recommendation, and the ranking is
of what is worth testing in a laboratory, not of what should be given to
anyone.

Input: `notebooks/utils/track2_step4_candidates.tsv`, where each row is one
hypothesis scored 0-3 on the five axes the procedure lists:

  mechanistic_relevance   proximity of the target to the BUB1B defect
  directionality          evidence the drug class acts in the desired direction
  experimental_evidence   BUB1B-deficient / SAC / aneuploidy model evidence
  human_exposure_evidence prior human exposure and characterised pharmacology
  tractability            an approved or otherwise obtainable compound exists

The weights encode the procedure's own rule: a drug is not prioritised merely
because it interacts with a protein somewhere in the BUB1B pathway, so
mechanistic relevance and directionality dominate, and drug availability
alone cannot carry a candidate (H5 is the worked example of that - best drug
availability, lowest mechanistic relevance, lowest rank).

Scores are a transparent weighted sum, not a trained model: there is no
labelled outcome data for MVA repurposing, so any fitted model would be
unvalidatable. The weights are stated here so a reader can disagree with them
explicitly and recompute.

Output (output/track2/):
  track2_step4_prioritised_hypotheses.csv
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

WEIGHTS = {
    "mechanistic_relevance": 0.35,
    "directionality": 0.20,
    "experimental_evidence": 0.20,
    "human_exposure_evidence": 0.15,
    "tractability": 0.10,
}
MAX_SCORE = 3.0

# Hypotheses that are not candidates but are reported for completeness:
# H3 is a mechanism-derived research question rather than a repurposing
# proposal. Every other hypothesis, H4 included, is ranked as a candidate.
NON_CANDIDATE_CONCEPTS = {"pharmacogenomic_caution"}


def score(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for axis in WEIGHTS:
        df[axis] = pd.to_numeric(df[axis], errors="coerce")
    # An axis that does not apply (n/a) contributes nothing and its weight is
    # dropped from the denominator, so H3 is not penalised for having no
    # "desired direction" to be concordant with.
    num = sum(df[a].fillna(0) * w for a, w in WEIGHTS.items())
    den = sum(df[a].notna() * w for a, w in WEIGHTS.items())
    df["priority_score"] = (num / (den * MAX_SCORE)).round(3)
    df["is_candidate"] = ~df["therapeutic_concept"].isin(NON_CANDIDATE_CONCEPTS)
    df = df.sort_values(["is_candidate", "priority_score"], ascending=[False, False])
    df["rank"] = range(1, len(df) + 1)
    return df


def build(candidates_path: Path, out_dir: Path, verbose: bool = True) -> pd.DataFrame:
    df = score(pd.read_csv(candidates_path, sep="\t"))
    out_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_dir / "track2_step4_prioritised_hypotheses.csv", index=False)
    if verbose:
        cols = ["rank", "hypothesis_id", "therapeutic_concept", "priority_score",
                *WEIGHTS, "is_candidate"]
        print(df[cols].to_string(index=False))
    return df


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--candidates", type=Path, default=here / "track2_step4_candidates.tsv")
    ap.add_argument("--out", type=Path, default=here.parents[1] / "output" / "track2")
    a = ap.parse_args()
    build(a.candidates, a.out)
