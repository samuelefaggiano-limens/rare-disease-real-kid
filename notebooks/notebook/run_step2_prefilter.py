"""Track 1, Step 2: pre-filtering of the patient VCF.

Intrinsic QC (FILTER/GQ/DP/GT/AB) + common-variant subtraction against dbSNP
Common. Same thresholds as track1_step2_exp01-vcf_prefiltering.ipynb.

Strategy: three SEQUENTIAL passes, no random lookups, bounded memory
(~250 MB), instead of one tabix query per variant.

  Pass 1  patient VCF -> temporary VCF with the intrinsic-QC survivors.
          Keeps the key (chrom,pos,ref,alt) of every survivor in memory.
  Pass 2  Scans dbSNP once, in binary mode, and discards from the set every
          key found there. What remains is the set of rare variants.
  Pass 3  Rewrites the temporary VCF keeping only the rare variants.

Why this design:
  - Loading all of dbSNP into a set of tuples costs ~287 bytes per variant,
    about 11 GB for its ~38M records: not viable.
  - One tabix query per variant means building ~2.6M iterators and parsing
    every dbSNP record as a full object. Here dbSNP is read as raw bytes,
    splitting only the first 5 fields and never touching its INFO.
  - Keys are stored exact, not hashed: a collision could delete the causal
    variant, and that risk is not worth saving ~100 MB.

Usage:
    python run_step2_prefilter.py
    python run_step2_prefilter.py --input path.vcf.gz --dbsnp common.vcf.gz -o out.vcf.gz
"""

import argparse
import gzip
import os
import sys
import time
from pathlib import Path

import cyvcf2

REPO_ROOT = Path(__file__).resolve().parents[2]

MIN_GQ = 20
MIN_DP = 10
HET_AB_MIN, HET_AB_MAX = 0.30, 0.70
HOM_ALT_AB_MIN = 0.90


def log(msg):
    print(msg, flush=True)


def variant_keys(chrom, pos, ref, alts):
    """Exact keys of a variant, one per alternate allele.

    Normalises the 'chr' prefix: the patient VCF and dbSNP both use numeric
    contigs, but it is better not to depend on that.
    """
    c = chrom[3:] if chrom.startswith("chr") else chrom
    return [f"{c}:{pos}:{ref}:{a}" for a in alts]


def pass1_intrinsic_qc(input_vcf, temp_vcf):
    """Intrinsic QC. Returns the set of surviving keys."""
    vcf_in = cyvcf2.VCF(input_vcf)
    vcf_in.add_info_to_header({
        "ID": "AB_FILTER",
        "Description": "Allelic Balance for Heterozygous calls",
        "Type": "Float",
        "Number": "1",
    })
    vcf_out = cyvcf2.Writer(temp_vcf, vcf_in)

    # Local counters: 3.4x faster than stats['key'] += 1 inside a 5M-iteration loop.
    total = dropped_filter = dropped_missing = 0
    dropped_gq = dropped_dp = dropped_gt = dropped_ab = 0
    survivors = set()
    start = time.time()

    for variant in vcf_in:
        total += 1

        # Filters run cheapest-first. variant.format() allocates a numpy
        # array per call, so AD is not fetched until GQ, DP and GT have
        # rejected what they can.
        if variant.FILTER is not None and variant.FILTER != "PASS":
            dropped_filter += 1
            continue

        gq_array = variant.format("GQ")
        if gq_array is None:
            dropped_missing += 1
            continue
        if gq_array[0][0] < MIN_GQ:
            dropped_gq += 1
            continue

        dp_array = variant.format("DP")
        if dp_array is None:
            dropped_missing += 1
            continue
        if dp_array[0][0] < MIN_DP:
            dropped_dp += 1
            continue

        gt_type = variant.gt_types[0]
        if gt_type == 0 or gt_type == 2:  # HOM_REF and UNKNOWN
            dropped_gt += 1
            continue

        ad_array = variant.format("AD")
        if ad_array is None:
            dropped_missing += 1
            continue

        ref_depth = ad_array[0][0]
        alt_depth = ad_array[0][1] if len(ad_array[0]) > 1 else 0
        calc_dp = ref_depth + alt_depth
        if calc_dp == 0:
            dropped_ab += 1
            continue

        ab = alt_depth / calc_dp
        if gt_type == 1:
            if ab < HET_AB_MIN or ab > HET_AB_MAX:
                dropped_ab += 1
                continue
        elif gt_type == 3:
            if ab < HOM_ALT_AB_MIN:
                dropped_ab += 1
                continue

        variant.INFO["AB_FILTER"] = round(float(ab), 3)
        vcf_out.write_record(variant)
        survivors.update(variant_keys(variant.CHROM, variant.POS, variant.REF, variant.ALT))

        if total % 500_000 == 0:
            log(f"    {total:,} variants read | {len(survivors):,} survivors "
                f"| {time.time() - start:.0f}s")

    vcf_in.close()
    vcf_out.close()

    stats = {
        "total": total,
        "dropped_filter": dropped_filter,
        "dropped_missing_metrics": dropped_missing,
        "dropped_gq": dropped_gq,
        "dropped_dp": dropped_dp,
        "dropped_gt": dropped_gt,
        "dropped_ab": dropped_ab,
    }
    return survivors, stats


def pass2_subtract_common(dbsnp_vcf, survivors):
    """Scans dbSNP once and removes from the set whatever is common.

    Reads in binary mode and splits only the first 5 fields: dbSNP's INFO
    is huge and not needed for anything.
    """
    before = len(survivors)
    discard = survivors.discard  # resolved once, not 38M times
    lines = 0
    start = time.time()

    with gzip.open(dbsnp_vcf, "rb") as f:
        for line in f:
            lines += 1
            if line[0] == 35:  # '#'
                continue

            parts = line.split(b"\t", 5)
            if len(parts) < 5:
                continue

            chrom = parts[0].decode("ascii", "ignore")
            if chrom.startswith("chr"):
                chrom = chrom[3:]
            pos = parts[1].decode("ascii", "ignore")
            ref = parts[3].decode("ascii", "ignore")
            alts = parts[4].decode("ascii", "ignore")

            prefix = f"{chrom}:{pos}:{ref}:"
            if "," in alts:
                for alt in alts.split(","):
                    discard(prefix + alt)
            else:
                discard(prefix + alts)

            if lines % 5_000_000 == 0:
                log(f"    {lines:,} dbSNP lines | {len(survivors):,} still rare "
                    f"| {time.time() - start:.0f}s")

    return before - len(survivors)


def pass3_write_rare(temp_vcf, output_vcf, survivors):
    """Rewrites the temporary VCF keeping only the rare variants."""
    vcf_in = cyvcf2.VCF(temp_vcf)
    vcf_out = cyvcf2.Writer(output_vcf, vcf_in)
    kept = 0

    for variant in vcf_in:
        keys = variant_keys(variant.CHROM, variant.POS, variant.REF, variant.ALT)
        # All of its keys entered the set in pass 1. If any is missing, pass 2
        # removed it as common.
        if any(k not in survivors for k in keys):
            continue
        vcf_out.write_record(variant)
        kept += 1

    vcf_in.close()
    vcf_out.close()
    return kept


def main():
    ap = argparse.ArgumentParser(description="Step 2: VCF pre-filtering (QC + rare variants)")
    ap.add_argument("--input", default=str(REPO_ROOT / "data" / "WGS_EX2312012_HGWCNDSX7.vcf.gz"))
    ap.add_argument("--dbsnp", default=str(REPO_ROOT / "data" / "common_all_20180418.vcf.gz"))
    ap.add_argument("-o", "--output", default=str(REPO_ROOT / "data" / "track1_step2_rare_variants.vcf.gz"))
    ap.add_argument("--keep-temp", action="store_true", help="do not delete the intermediate VCF")
    args = ap.parse_args()

    for path, label in ((args.input, "patient VCF"), (args.dbsnp, "dbSNP Common")):
        if not os.path.exists(path):
            sys.exit(f"[-] {label} not found: {path}")

    temp_vcf = args.output.replace(".vcf.gz", "") + ".intrinsic_qc.tmp.vcf.gz"
    t0 = time.time()

    log(f"[1/3] Intrinsic QC on {args.input}")
    survivors, stats = pass1_intrinsic_qc(args.input, temp_vcf)
    passed_intrinsic = stats["total"] - sum(v for k, v in stats.items() if k.startswith("dropped"))
    log(f"      {passed_intrinsic:,} of {stats['total']:,} pass QC "
        f"({time.time() - t0:.0f}s)")

    log(f"[2/3] Common-variant subtraction against {args.dbsnp}")
    dropped_common = pass2_subtract_common(args.dbsnp, survivors)
    log(f"      {dropped_common:,} keys discarded as common ({time.time() - t0:.0f}s)")

    log(f"[3/3] Writing rare variants to {args.output}")
    kept = pass3_write_rare(temp_vcf, args.output, survivors)

    if not args.keep_temp:
        for suffix in ("", ".tbi", ".csi"):
            try:
                os.remove(temp_vcf + suffix)
            except OSError:
                pass

    elapsed = time.time() - t0
    log("\n--- FILTERING STATISTICS ---")
    log(f"Total variants:                {stats['total']:,}")
    log(f"Dropped by FILTER:             {stats['dropped_filter']:,}")
    log(f"Dropped, missing metrics:      {stats['dropped_missing_metrics']:,}")
    log(f"Dropped by GQ < {MIN_GQ}:           {stats['dropped_gq']:,}")
    log(f"Dropped by DP < {MIN_DP}:           {stats['dropped_dp']:,}")
    log(f"Dropped by GT (ref/unk):       {stats['dropped_gt']:,}")
    log(f"Dropped by allelic balance:    {stats['dropped_ab']:,}")
    log("--------------------------------")
    log(f"Pass intrinsic QC:             {passed_intrinsic:,}")
    log(f"Dropped as common:             {dropped_common:,}")
    log(f"Final rare variants:           {kept:,} "
        f"({kept / stats['total'] * 100:.1f}% of total)")
    log(f"\nCompleted in {elapsed / 60:.1f} min -> {args.output}")


if __name__ == "__main__":
    main()
