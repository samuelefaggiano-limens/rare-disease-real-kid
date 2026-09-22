"""Track 1, Step 5: epcr score and submission CSV.

Input: data/track1_step4_candidates.tsv (Step 4). Output: the CSV in the exact
format the Space validates (proband PROBAND01, chromosomes with the chr
prefix, <= 10 rows, epcr in (0,1], finding_type).

How each row is built:
  * Per gene. If a gene has >= 2 heterozygous variants, the two best go
    together in ONE row as a compound-heterozygous pair (columns _1 and _2),
    which is the expected model in MVA. A homozygous-alt variant goes alone.
    A gene with a single het goes alone, with an inheritance penalty.
  * raw epcr = weighted geometric mean of phenotype, pathogenicity,
    inheritance and frequency (src/epcr.py), with a x1.2 bonus for genes in
    the MVA panel or annotated to the MVA disease in HPO.
  * Separation: --max-rows rows are kept and each one is capped at
    --tail-ratio times the previous one. F-max sweeps every distinct epcr as
    a threshold, so near-tied rows sink precision; a geometric decay forces a
    ranking decision.

At the end it prints an "if row k were the true one" table: rank points and
F-max the leaderboard would return, with the same logic as its evaluation.py.
It shows what every extra row costs.

Usage:
    python track1_step5_submission.py --username <hf-user> --model variant-filter-v1
    python track1_step5_submission.py --max-rows 3 --tail-ratio 0.5
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from src.epcr import (CandidateRow, assign_epcr, combine, frequency_score,  # noqa: E402
                      inheritance_score, pathogenicity_score, what_if_table)

DEFAULT_INPUT = REPO_ROOT / "data" / "track1_step4_candidates.tsv"
OUTPUT_DIR = REPO_ROOT / "output" / "track1"
PROBAND_ID = "PROBAND01"          # the only id the Space accepts (tabs/submit_track1.py)
MVA_BONUS = 1.2


def fnum(v):
    try:
        return float(v) if v not in ("", None, "None") else None
    except ValueError:
        return None


def fbool(v):
    return str(v).lower() in ("true", "1", "yes")


def score_variant(v: dict) -> dict:
    comps = {
        "phenotype": fnum(v.get("pheno_score")) or 0.0,
        "pathogenicity": pathogenicity_score(v.get("PUTATIVE_IMPACT"), fnum(v.get("CADD_PHRED")),
                                             fnum(v.get("REVEL")), v.get("CLINVAR_SIG"),
                                             avi_phred=fnum(v.get("AVI_PHRED"))),
        "inheritance": inheritance_score(v.get("zygosity", ""), int(fnum(v.get("n_het_in_gene")) or 0)),
        "frequency": frequency_score(fnum(v.get("GNOMAD_AF"))),
    }
    bonus = MVA_BONUS if (fbool(v.get("mva_panel")) or fbool(v.get("mva_disease_annotated"))) else 1.0
    return {"components": comps, "bonus": bonus, "raw": combine(comps, bonus=bonus)}


def variant_key(v: dict, chr_prefix: bool):
    c = v["chrom"]
    c = f"chr{c}" if chr_prefix and not c.startswith("chr") else c
    return (c, int(v["pos"]), v["ref"].upper(), v["alt"].upper())


def note_for(v: dict, s: dict) -> str:
    c = s["components"]
    avi = v.get("AVI_PHRED")
    return (f"{v['gene']} {v.get('CONSEQUENCE', '')} {v.get('PUTATIVE_IMPACT', '')} {v['zygosity']}"
            f" pheno={c['phenotype']:.2f} patho={c['pathogenicity']:.2f} inh={c['inheritance']:.2f}"
            f" freq={c['frequency']:.2f}" + (f" AVI_PHRED={avi}" if avi not in ("", None) else "")
            + (" MVA-panel" if s["bonus"] > 1 else ""))


def build_rows(variants: list[dict], chr_prefix: bool, max_af: float = 0.01) -> list[CandidateRow]:
    """Groups by gene and builds rows. Discards variants with gnomAD AF > max_af.

    Step 2's dbSNP Common subtraction misses many indels (dbSNP represents
    them with a different normalisation, so the exact chrom:pos:ref:alt key
    does not match): Step 3A found "rare" indels with gnomAD AF 0.34-0.54.
    The threshold is the procedure's 1 %; it is not lowered further because
    of the founder-effect caution in 00-references.md section 4.
    """
    by_gene: dict[str, list[tuple[dict, dict]]] = defaultdict(list)
    for v in variants:
        if not v.get("gene"):
            continue
        af = max(fnum(v.get("GNOMAD_AF")) or 0.0, fnum(v.get("GNOMAD_EXOME_AF")) or 0.0)
        if af > max_af:
            continue
        by_gene[v["gene"]].append((v, score_variant(v)))

    rows: list[CandidateRow] = []
    for gene, items in by_gene.items():
        items.sort(key=lambda x: -x[1]["raw"])
        hets = [(v, s) for v, s in items if v["zygosity"] == "HET"]
        homs = [(v, s) for v, s in items if v["zygosity"] == "HOM_ALT"]
        used = set()
        for v, s in homs:                                   # homozygous: its own row
            rows.append(CandidateRow((variant_key(v, chr_prefix),), s["raw"], gene, note_for(v, s)))
            used.add(id(v))
        if len(hets) >= 2:                                  # compound-heterozygous pair
            (v1, s1), (v2, s2) = hets[0], hets[1]
            raw = (s1["raw"] * s2["raw"]) ** 0.5             # the pair is dragged down by its weaker member
            rows.append(CandidateRow((variant_key(v1, chr_prefix), variant_key(v2, chr_prefix)),
                                     round(raw, 4), gene, f"compound-het | {note_for(v1, s1)} | {note_for(v2, s2)}"))
            used.update((id(v1), id(v2)))
        for v, s in hets:
            if id(v) in used:
                continue
            rows.append(CandidateRow((variant_key(v, chr_prefix),), s["raw"], gene, note_for(v, s)))
    return rows


def main():
    ap = argparse.ArgumentParser(description="Step 5: epcr and submission CSV")
    ap.add_argument("--input", default=str(DEFAULT_INPUT))
    ap.add_argument("--username", default="team", help="HF user or team name for the CSV file name")
    ap.add_argument("--model", default="hpo-epcr-v1", help="model name for the CSV file name")
    ap.add_argument("--max-rows", type=int, default=5, help="rows to submit (maximum 10)")
    ap.add_argument("--tail-ratio", type=float, default=0.6, help="each row <= ratio x the previous one")
    ap.add_argument("--top-cap", type=float, default=0.95)
    ap.add_argument("--no-chr-prefix", action="store_true", help="emit '1' instead of 'chr1'")
    ap.add_argument("--secondary-from", type=int, default=None,
                    help="rows from this position onwards are marked finding_type=secondary")
    ap.add_argument("--max-af", type=float, default=0.01,
                    help="discard variants with gnomAD AF (genomes or exomes) above this")
    args = ap.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"[-] {args.input} not found. Run Step 4 first.")
    if not 1 <= args.max_rows <= 10:
        sys.exit("[-] --max-rows must be between 1 and 10 (leaderboard limit)")

    with open(args.input, newline="") as fh:
        variants = list(csv.DictReader(fh, delimiter="\t"))
    print(f"[+] {len(variants):,} candidate variants read from {args.input}")

    rows = build_rows(variants, chr_prefix=not args.no_chr_prefix, max_af=args.max_af)
    print(f"    {len(rows):,} candidate rows built (genes with >=2 het go as a pair; "
          f"gnomAD AF > {args.max_af} discarded)")
    ranked = assign_epcr(rows, max_rows=args.max_rows, tail_ratio=args.tail_ratio, top_cap=args.top_cap)
    if args.secondary_from:
        for r in ranked[args.secondary_from - 1:]:
            r.finding_type = "secondary"

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out_path = OUTPUT_DIR / f"{args.username}_{args.model}.csv"
    fields = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2", "ref_2", "alt_2",
              "epcr", "finding_type", "notes"]
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in ranked:
            v1 = r.variants[0]
            v2 = r.variants[1] if len(r.variants) > 1 else ("", "", "", "")
            w.writerow({"proband_id": PROBAND_ID,
                        "chrom_1": v1[0], "pos_1": v1[1], "ref_1": v1[2], "alt_1": v1[3],
                        "chrom_2": v2[0], "pos_2": v2[1], "ref_2": v2[2], "alt_2": v2[3],
                        "epcr": f"{r.epcr:.4f}", "finding_type": r.finding_type, "notes": r.notes})
    print(f"[OK] Submission -> {out_path}")

    print("\n--- SUBMITTED ROWS ---")
    for i, r in enumerate(ranked, 1):
        vs = " + ".join(f"{c}:{p} {a}>{b}" for c, p, a, b in r.variants)
        print(f"#{i} epcr={r.epcr:.4f} raw={r.raw_score:.4f} {r.gene:<10} {vs}")
        print(f"    {r.notes}")

    print("\n--- IF ROW k WERE THE TRUE ONE (same logic as the Space's evaluation.py) ---")
    print(f"{'k':>2} {'epcr':>7} {'rank_pts':>8} {'F-max':>6} {'thresh':>7}")
    for w in what_if_table(ranked):
        print(f"{w['if_true_row']:>2} {w['epcr']:>7.4f} {w['rank_points']:>8.1f} {w['f_max']:>6.3f} "
              f"{w['f_max_threshold'] if w['f_max_threshold'] is not None else '-':>7}")


if __name__ == "__main__":
    main()
