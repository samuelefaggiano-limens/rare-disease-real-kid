"""Track 1, Step 3B: regulatory enrichment with AlphaGenome Atlas.

Adds the AVI (AlphaGenome Variant Impact) scores to every variant coming out
of Step 3A. Pure enrichment step: no variant is discarded.

INFO tags added:
  AVI_SCORE           raw model score
  AVI_QUANTILE        calibrated quantile (1 - CDF) of the variant against all
                      genome-wide SNVs
  AVI_PHRED           -10*log10(1 - quantile). >= 20 is the top 1 % of the
                      genome, >= 30 the top 0.1 %
  REGULATORY_FEATURE  modality with the largest attribution (MERGED_SPLICING,
                      ALPHAMISSENSE, MAX_ABS_RNA_SEQ, ...) from the
                      AVI_SCORE_FEATURE_IMPORTANCE scorer
  AVI_STATUS          'ok' or the reason the variant could not be scored
                      (e.g. ref_mismatch when the REF allele does not match the
                      server's GRCh38 reference)

TARGET_GENE, planned in docs/01-procedure.md, is not returned by the AVI
scorers; the gene remains snpEff's GENE_SYMBOL from Step 3A.

How the API is used (alphagenome 0.9):
  - The client is created with atlas.create(api_key); a bare AtlasClient()
    does not work (it needs the gRPC stub). The key is read from
    ALPHAGENOME_API_KEY in ~/.env (AlphaGenome skills convention) or from
    the environment.
  - query_variants() returns {scorer: AnnData}; obs['variant'] identifies the
    row, X[:,0] is the score and layers['quantiles'] the quantile.
  - The server validates the REF allele: one mis-referenced variant fails the
    WHOLE batch. So when a batch fails with ValueError it is retried variant
    by variant, and only the offenders are marked in AVI_STATUS.
  - Infrastructure errors (network, quota) are retried with backoff and, if
    they persist, the script ABORTS instead of writing a half-annotated VCF.

Usage:
    python step3b_alphagenome.py
    python step3b_alphagenome.py --batch-size 100 --workers 8
"""

from __future__ import annotations

import argparse
import gzip
import math
import os
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_INPUT = REPO_ROOT / "track1_step3a_snpeff_filter.vcf.gz"
DEFAULT_OUTPUT = REPO_ROOT / "track1_step3b_regulatory_annotated.vcf.gz"
SCORERS = ["AVI_SCORE", "AVI_SCORE_FEATURE_IMPORTANCE"]
MAX_RETRIES = 5

NEW_HEADERS = [
    '##INFO=<ID=AVI_SCORE,Number=1,Type=Float,Description="AlphaGenome Variant Impact raw score">',
    '##INFO=<ID=AVI_QUANTILE,Number=1,Type=Float,Description="AlphaGenome Variant Impact calibrated tail quantile (1-CDF) vs genome-wide SNVs">',
    '##INFO=<ID=AVI_PHRED,Number=1,Type=Float,Description="-10*log10(1-AVI_QUANTILE); >=20 is top 1 percent genome-wide">',
    '##INFO=<ID=REGULATORY_FEATURE,Number=1,Type=String,Description="Top contributing AlphaGenome modality (feature importance)">',
    '##INFO=<ID=AVI_STATUS,Number=1,Type=String,Description="ok, or reason the variant could not be scored">',
]


def log(msg):
    print(msg, flush=True)


def load_api_key() -> str:
    key = os.environ.get("ALPHAGENOME_API_KEY")
    if not key:
        env_path = os.path.expanduser("~/.env")
        if os.path.exists(env_path):
            with open(env_path) as fh:
                for line in fh:
                    if line.startswith("ALPHAGENOME_API_KEY="):
                        key = line.split("=", 1)[1].strip().strip('"').strip("'")
                        break
    if not key:
        sys.exit("[-] ALPHAGENOME_API_KEY is missing (environment or ~/.env).\n"
                 "    Save it with: printf 'Enter ALPHAGENOME_API_KEY: ' && read -s v && echo \"ALPHAGENOME_API_KEY=$v\" >> ~/.env")
    return key


def build_client():
    try:
        from alphagenome.atlas import atlas
    except ImportError as exc:
        sys.exit(f"[-] Could not import alphagenome ({exc}). Install it with: pip install alphagenome")
    try:
        return atlas.create(load_api_key())
    except Exception as exc:  # noqa: BLE001
        sys.exit(f"[-] Could not create the Atlas client ({exc}).")


def format_variant_id(chrom, pos, ref, alt):
    c = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return f"{c}:{pos}:{ref}:{alt}"


def to_atlas_str(vid: str) -> str:
    """chr15:100:C:T -> chr15:100:C>T (genome.Variant.from_str format)."""
    c, p, r, a = vid.split(":")
    return f"{c}:{p}:{r}>{a}"


def parse_results(res) -> dict[str, dict]:
    """{variant_str: {AVI_SCORE, AVI_QUANTILE, AVI_PHRED, REGULATORY_FEATURE}}."""
    out: dict[str, dict] = {}
    score_ad = res.get("AVI_SCORE")
    if score_ad is not None:
        X = score_ad.X.toarray() if hasattr(score_ad.X, "toarray") else score_ad.X
        q_layer = score_ad.layers.get("quantiles") if hasattr(score_ad, "layers") else None
        if q_layer is not None and hasattr(q_layer, "toarray"):
            q_layer = q_layer.toarray()
        for i, v in enumerate(score_ad.obs["variant"]):
            d = out.setdefault(str(v), {})
            d["AVI_SCORE"] = f"{float(X[i, 0]):.4f}"
            if q_layer is not None:
                q = float(q_layer[i, 0])
                d["AVI_QUANTILE"] = f"{q:.6f}"
                d["AVI_PHRED"] = f"{-10 * math.log10(max(1 - q, 1e-12)):.2f}"
    fi_ad = res.get("AVI_SCORE_FEATURE_IMPORTANCE")
    if fi_ad is not None:
        X = fi_ad.X.toarray() if hasattr(fi_ad.X, "toarray") else fi_ad.X
        names = list(fi_ad.var["name"])
        for i, v in enumerate(fi_ad.obs["variant"]):
            row = X[i]
            j = max(range(len(names)), key=lambda k: abs(float(row[k])))
            out.setdefault(str(v), {})["REGULATORY_FEATURE"] = names[j]
    return out


def query_batch(client, variant_ids: list[str], workers: int) -> tuple[dict[str, dict], dict[str, str]]:
    """Returns (annotations, per-variant errors). Aborts if the API does not respond."""
    from alphagenome.data import genome

    def run(ids):
        vs = [genome.Variant.from_str(to_atlas_str(v)) for v in ids]
        return client.query_variants(vs, requested_scorers=SCORERS, progress_bar=False, max_workers=workers)

    errors: dict[str, str] = {}
    for attempt in range(MAX_RETRIES):
        try:
            return parse_results(run(variant_ids)), errors
        except (ValueError, NotImplementedError) as exc:
            # Server-side validation (REF mismatch, unknown contig, unsupported
            # variant type): it concerns one variant but brings the whole
            # batch down. Resolve variant by variant.
            if len(variant_ids) == 1:
                errors[variant_ids[0]] = classify_error(exc)
                return {}, errors
            annotations: dict[str, dict] = {}
            for vid in variant_ids:
                a, e = query_batch(client, [vid], workers)
                annotations.update(a)
                errors.update(e)
            return annotations, errors
        except Exception as exc:  # noqa: BLE001  network / gRPC / quota errors
            wait = 2 ** attempt + 1
            log(f"    warning: batch failed ({type(exc).__name__}: {str(exc)[:120]}); retrying in {wait}s")
            time.sleep(wait)
    sys.exit(f"[-] The Atlas API did not respond after {MAX_RETRIES} attempts. "
             f"The output VCF is incomplete; do not use it.")


def classify_error(exc: Exception) -> str:
    msg = str(exc)
    if isinstance(exc, NotImplementedError):
        return "unsupported_by_atlas"
    if "reference base" in msg:
        return "ref_mismatch"
    if "chromosome" in msg.lower() or "contig" in msg.lower():
        return "unknown_contig"
    return "invalid:" + msg[:60].replace(";", ",").replace(" ", "_")


def is_snv(vid: str) -> bool:
    _, _, ref, alt = vid.split(":")
    return len(ref) == 1 and len(alt) == 1


def flush_batch(client, batch, outfile, workers, stats):
    # Atlas (alphagenome 0.9, 2026-09) answers NotImplementedError to every
    # indel: they are marked without spending a call. Only SNVs are queried.
    ids = [vid for _, vid in batch if is_snv(vid)]
    annotations, errors = query_batch(client, ids, workers) if ids else ({}, {})
    for parts, vid in batch:
        if not is_snv(vid):
            errors[vid] = "indel_unsupported"
        data = annotations.get(to_atlas_str(vid)) or annotations.get(vid)
        additions = []
        if data:
            additions.extend(f"{k}={v}" for k, v in data.items())
            additions.append("AVI_STATUS=ok")
            stats["scored"] += 1
        else:
            additions.append(f"AVI_STATUS={errors.get(vid, 'no_result')}")
            stats["unscored"] += 1
        parts[7] = parts[7] + ";" + ";".join(additions)
        outfile.write("\t".join(parts) + "\n")
    return len(batch)


def run_step3b(input_vcf, output_vcf, batch_size, workers):
    client = build_client()
    log(f"[+] Enriching {input_vcf} -> {output_vcf}")

    batch, total, header_done = [], 0, False
    stats = {"scored": 0, "unscored": 0}
    start = time.time()
    open_in = gzip.open if str(input_vcf).endswith(".gz") else open

    with open_in(input_vcf, "rt", encoding="utf-8") as infile, \
         gzip.open(output_vcf, "wt", encoding="utf-8") as outfile:
        for line in infile:
            if line.startswith("##"):
                outfile.write(line)
                continue
            if line.startswith("#CHROM"):
                for nh in NEW_HEADERS:
                    outfile.write(nh + "\n")
                outfile.write(line)
                header_done = True
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, _, ref, alt = parts[:5]
            batch.append((parts, format_variant_id(chrom, pos, ref, alt.split(",", 1)[0])))
            if len(batch) >= batch_size:
                total += flush_batch(client, batch, outfile, workers, stats)
                batch = []
                rate = total / max(time.time() - start, 1e-9)
                log(f"    {total:,} variants ({rate:.1f}/s) | scored {stats['scored']:,} | unscored {stats['unscored']:,}")
        if batch:
            total += flush_batch(client, batch, outfile, workers, stats)

    if not header_done:
        sys.exit(f"[-] {input_vcf} has no #CHROM line; the input VCF is empty or corrupt.")

    log(f"\n[OK] Step 3B: {total:,} variants written, {stats['scored']:,} with AVI, "
        f"{stats['unscored']:,} unscored (see AVI_STATUS) in {(time.time() - start) / 60:.1f} min -> {output_vcf}")
    return total


def main():
    ap = argparse.ArgumentParser(description="Step 3B: enrichment with AlphaGenome Atlas (AVI)")
    ap.add_argument("--input", default=str(DEFAULT_INPUT))
    ap.add_argument("-o", "--output", default=str(DEFAULT_OUTPUT))
    ap.add_argument("--batch-size", type=int, default=100, help="variants per call")
    ap.add_argument("--workers", type=int, default=8, help="in-flight gRPC requests within a batch")
    args = ap.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"[-] Input VCF not found: {args.input}\n    Run Step 3A first.")
    run_step3b(args.input, args.output, args.batch_size, args.workers)


if __name__ == "__main__":
    main()
