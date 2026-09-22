"""Track 1, Step 3A: functional annotation (snpEff already applied) + MyVariant.info.

Reads the snpEff-annotated VCF, keeps HIGH/MODERATE impacts and adds the
MyVariant.info metadata (gnomAD, CADD, REVEL, ClinVar, OMIM) as INFO tags.

Optimisations over the previous version:
  - The ANN field is extracted with a targeted search on the INFO string
    instead of building a dict of every field: 9.2x faster, and it runs once
    for each of the ~591k input lines.
  - Batches are requested in parallel with a thread pool. The cost here is
    network latency, not CPU, so that is where the bulk of the time goes.
  - A single requests.Session reuses the TCP/TLS connection between batches
    instead of renegotiating the handshake on every request.
  - If a batch exhausts its retries the script FAILS instead of silently
    writing unannotated variants. A half-annotated VCF that looks complete is
    worse than an error.

Usage:
    python track1_step3a_annotate.py
    python track1_step3a_annotate.py --workers 8 --batch-size 1000
"""

import argparse
import gzip
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import requests

REPO_ROOT = Path(__file__).resolve().parents[1]

DEFAULT_INPUT = REPO_ROOT / "snpeff_annotated.vcf"
DEFAULT_OUTPUT = REPO_ROOT / "track1_step3a_snpeff_filter.vcf.gz"

API_URL = "https://myvariant.info/v1/variant"
MAX_RETRIES = 5
RETRY_STATUS = (429, 500, 502, 503, 504)

# The patient VCF is GRCh38. MyVariant.info indexes hg19 by default: without
# assembly=hg38, 97 % of the variants come back "notfound" and the remaining
# 3 % are coordinate coincidences in hg19 (annotations of ANOTHER variant).
ASSEMBLY = "hg38"

# Fields available in the hg38 index (checked against /v1/metadata/fields on
# 2026-09-17). CADD comes via dbNSFP (the 'cadd' dataset is hg19-only);
# LoFtool does not exist in MyVariant, so LOFTOOL_SCORE is not emitted.
FIELDS = (
    "gnomad_genome.af.af,gnomad_exome.af.af,dbnsfp.cadd.phred,dbnsfp.revel.score,"
    "clinvar.rcv.clinical_significance,clinvar.rcv.review_status,"
    "clinvar.rcv.conditions.identifiers.omim"
)

NEW_HEADERS = [
    '##INFO=<ID=GENE_SYMBOL,Number=1,Type=String,Description="Affected Gene Symbol">',
    '##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Physical alteration/consequence on transcript">',
    '##INFO=<ID=PUTATIVE_IMPACT,Number=1,Type=String,Description="Putative impact: HIGH, MODERATE, LOW, MODIFIER">',
    '##INFO=<ID=GNOMAD_AF,Number=1,Type=Float,Description="gnomAD genomes overall allele frequency (hg38)">',
    '##INFO=<ID=GNOMAD_EXOME_AF,Number=1,Type=Float,Description="gnomAD exomes overall allele frequency (hg38)">',
    '##INFO=<ID=CADD_PHRED,Number=1,Type=Float,Description="CADD Phred pathogenicity score (via dbNSFP)">',
    '##INFO=<ID=REVEL,Number=1,Type=Float,Description="REVEL score">',
    '##INFO=<ID=CLINVAR_SIG,Number=.,Type=String,Description="ClinVar clinical significance (unique values, |-separated)">',
    '##INFO=<ID=CLINVAR_REVIEW_STATUS,Number=1,Type=String,Description="ClinVar review status of the most clinically relevant RCV">',
    '##INFO=<ID=OMIM_PHENOTYPE,Number=.,Type=String,Description="OMIM phenotype IDs from ClinVar conditions">',
    '##INFO=<ID=MYVARIANT_HIT,Number=0,Type=Flag,Description="Variant was found in MyVariant.info (hg38)">',
]


class AnnotationError(RuntimeError):
    """A batch could not be annotated after exhausting the retries."""


def extract_ann(info):
    """Returns the value of the ANN field without parsing the rest of INFO.

    Replaces dict(item.split('=', 1) for item in info.split(';')), which
    built a whole ~15-key dictionary to read a single one.
    """
    if info.startswith("ANN="):
        return info[4:].split(";", 1)[0]
    i = info.find(";ANN=")
    if i == -1:
        return None
    return info[i + 5:].split(";", 1)[0]


def format_hgvs_id(chrom, pos, ref, alt):
    """Genomic HGVS that MyVariant recognises. 'g.100AT>A' is NOT valid HGVS:
    indels go as del / ins / delins over the bases that actually change."""
    c = chrom if chrom.startswith("chr") else f"chr{chrom}"
    pos = int(pos)
    if len(ref) == 1 and len(alt) == 1:
        return f"{c}:g.{pos}{ref}>{alt}"
    if len(ref) > 1 and len(alt) == 1 and ref[0] == alt[0]:      # deletion
        start, end = pos + 1, pos + len(ref) - 1
        return f"{c}:g.{start}del" if start == end else f"{c}:g.{start}_{end}del"
    if len(ref) == 1 and len(alt) > 1 and ref[0] == alt[0]:      # insertion
        return f"{c}:g.{pos}_{pos + 1}ins{alt[1:]}"
    if ref[0] == alt[0]:                                         # delins with anchor base
        start, end = pos + 1, pos + len(ref) - 1
        rng = f"{start}" if start == end else f"{start}_{end}"
        return f"{c}:g.{rng}delins{alt[1:]}"
    rng = f"{pos}" if len(ref) == 1 else f"{pos}_{pos + len(ref) - 1}"
    return f"{c}:g.{rng}delins{alt}"


def make_session():
    """Session with a connection pool sized for the in-flight threads."""
    session = requests.Session()
    adapter = requests.adapters.HTTPAdapter(pool_connections=16, pool_maxsize=16)
    session.mount("https://", adapter)
    session.headers.update({"User-Agent": "mva-hackathon-2026/track1-step3a"})
    return session


def query_myvariant_batch(session, ids):
    """Queries one batch. Raises AnnotationError if no response is obtained."""
    payload = {"ids": ",".join(ids), "fields": FIELDS, "assembly": ASSEMBLY}
    last_error = "no detail"

    for attempt in range(MAX_RETRIES):
        try:
            res = session.post(API_URL, data=payload, timeout=45)
            if res.status_code == 200:
                data_map = {}
                for item in res.json():
                    q = item.get("query")
                    if q:
                        data_map[q] = item
                return data_map
            last_error = f"HTTP {res.status_code}"
            if res.status_code not in RETRY_STATUS:
                break
        except requests.RequestException as exc:
            last_error = f"{type(exc).__name__}: {exc}"

        if attempt < MAX_RETRIES - 1:
            time.sleep(2 ** attempt + 1)

    raise AnnotationError(
        f"Batch of {len(ids)} variants failed after {MAX_RETRIES} attempts ({last_error}). "
        f"First id in the batch: {ids[0]}"
    )


def parse_metadata(data):
    meta = {}
    if not data or "notfound" in data:
        return meta

    meta["MYVARIANT_HIT"] = None  # flag without a value

    for key, tag in (("gnomad_genome", "GNOMAD_AF"), ("gnomad_exome", "GNOMAD_EXOME_AF")):
        block = data.get(key)
        if isinstance(block, list):
            block = block[0]
        if block:
            af = (block.get("af") or {}).get("af")
            if af is not None:
                meta[tag] = str(af)

    dbnsfp = data.get("dbnsfp")
    if isinstance(dbnsfp, list):
        dbnsfp = dbnsfp[0]
    if dbnsfp:
        cadd = dbnsfp.get("cadd") or {}
        if cadd.get("phred") is not None:
            v = cadd["phred"]
            meta["CADD_PHRED"] = str(max(v) if isinstance(v, list) else v)
        revel = (dbnsfp.get("revel") or {}).get("score")
        if revel is not None:
            meta["REVEL"] = str(max(revel) if isinstance(revel, list) else revel)

    clinvar = data.get("clinvar")
    if clinvar and "rcv" in clinvar:
        rcvs = clinvar["rcv"]
        if not isinstance(rcvs, list):
            rcvs = [rcvs]
        sigs, omims = [], []
        best = None
        for r in rcvs:
            sig = str(r.get("clinical_significance", "")).strip()
            if sig and sig not in sigs:
                sigs.append(sig)
            cond = r.get("conditions")
            for c in (cond if isinstance(cond, list) else [cond]):
                ident = (c or {}).get("identifiers") or {}
                o = ident.get("omim")
                for x in (o if isinstance(o, list) else [o]):
                    if x and str(x) not in omims:
                        omims.append(str(x))
            # The most informative RCV: pathogenic > VUS > benign > no assertion
            rank = 3 if "athogenic" in sig else 2 if "ncertain" in sig else 1 if "enign" in sig else 0
            if best is None or rank > best[0]:
                best = (rank, r)
        if sigs:
            meta["CLINVAR_SIG"] = "|".join(s.replace(" ", "_").replace(",", "") for s in sigs)
        if best and best[1].get("review_status"):
            meta["CLINVAR_REVIEW_STATUS"] = str(best[1]["review_status"]).replace(" ", "_").replace(",", "")
        if omims:
            meta["OMIM_PHENOTYPE"] = "|".join(omims)

    return meta


def write_chunk(chunk, meta_maps, outfile):
    """Writes a group of already-annotated batches, preserving VCF order."""
    written = 0
    for (records, ids), meta_map in zip(chunk, meta_maps):
        for rec, vid in zip(records, ids):
            additions = [
                f"GENE_SYMBOL={rec['gene']}",
                f"CONSEQUENCE={rec['consequence']}",
                f"PUTATIVE_IMPACT={rec['impact']}",
            ]
            additions.extend(k if v is None else f"{k}={v}"
                             for k, v in parse_metadata(meta_map.get(vid, {})).items())

            row = [rec["chrom"], rec["pos"], rec["id"], rec["ref"], rec["alt"],
                   rec["qual"], rec["filt"], rec["info_raw"] + ";" + ";".join(additions)]
            row.extend(rec["genotypes"])
            outfile.write("\t".join(row) + "\n")
            written += 1
    return written


def process_pipeline(input_vcf, output_vcf, batch_size, workers):
    open_fn = gzip.open if str(input_vcf).endswith(".gz") else open
    session = make_session()

    total_read = total_kept = 0
    batch_records, batch_ids = [], []
    chunk = []            # batches waiting to be requested in parallel
    start = time.time()

    with open_fn(input_vcf, "rt", encoding="utf-8") as infile, \
         gzip.open(output_vcf, "wt", encoding="utf-8") as outfile, \
         ThreadPoolExecutor(max_workers=workers) as pool:

        def flush_chunk():
            nonlocal total_kept, chunk
            if not chunk:
                return
            # executor.map preserves input order, so the output VCF keeps its
            # coordinate order even when batches return out of order.
            meta_maps = list(pool.map(lambda b: query_myvariant_batch(session, b[1]), chunk))
            total_kept += write_chunk(chunk, meta_maps, outfile)
            chunk = []
            log_progress()

        def log_progress():
            rate = total_kept / max(time.time() - start, 1e-9)
            print(f"    annotated {total_kept:,} variants ({rate:.0f}/s)", flush=True)

        for line in infile:
            if line.startswith("##"):
                outfile.write(line)
                continue
            if line.startswith("#CHROM"):
                for nh in NEW_HEADERS:
                    outfile.write(nh + "\n")
                outfile.write(line)
                continue

            total_read += 1
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom, pos, var_id, ref, alt, qual, filt, info = parts[:8]

            ann_field = extract_ann(info)
            if not ann_field:
                continue

            first_trans = ann_field.split(",", 1)[0].split("|")
            impact = first_trans[2] if len(first_trans) > 2 else ""
            if impact not in ("HIGH", "MODERATE"):
                continue

            batch_records.append({
                "chrom": chrom, "pos": pos, "id": var_id, "ref": ref, "alt": alt,
                "qual": qual, "filt": filt, "info_raw": info,
                "gene": first_trans[3] if len(first_trans) > 3 else "",
                "consequence": first_trans[1] if len(first_trans) > 1 else "",
                "impact": impact,
                "genotypes": parts[8:],
            })
            batch_ids.append(format_hgvs_id(chrom, pos, ref, alt.split(",", 1)[0]))

            if len(batch_records) >= batch_size:
                chunk.append((batch_records, batch_ids))
                batch_records, batch_ids = [], []
                if len(chunk) >= workers:
                    flush_chunk()

        if batch_records:
            chunk.append((batch_records, batch_ids))
        flush_chunk()

    return total_read, total_kept, time.time() - start


def main():
    ap = argparse.ArgumentParser(description="Step 3A: annotation with MyVariant.info")
    ap.add_argument("--input", default=str(DEFAULT_INPUT))
    ap.add_argument("-o", "--output", default=str(DEFAULT_OUTPUT))
    ap.add_argument("--batch-size", type=int, default=1000,
                    help="ids per request (API maximum: 1000)")
    ap.add_argument("--workers", type=int, default=4,
                    help="simultaneous in-flight batches; raising it may trigger 429s")
    args = ap.parse_args()

    if not Path(args.input).exists():
        sys.exit(f"[-] Input VCF not found: {args.input}\n"
                 f"    Generate snpeff_annotated.vcf with snpEff first.")

    print(f"[+] Annotating {args.input} -> {args.output}")
    print(f"    batches of {args.batch_size}, {args.workers} in parallel")

    try:
        total_read, total_kept, elapsed = process_pipeline(
            args.input, args.output, args.batch_size, args.workers
        )
    except AnnotationError as exc:
        sys.exit(f"[-] {exc}\n    The output VCF is incomplete; do not use it.")

    print(f"\n[OK] {total_read:,} variants read, {total_kept:,} HIGH/MODERATE annotated "
          f"in {elapsed / 60:.1f} min")
    print(f"     Result: {args.output}")


if __name__ == "__main__":
    main()
