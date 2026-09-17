# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A submission for the SageBio "Rare Disease, Real Kid" MVA Hackathon 2026.

- **Track 1 — Variant Prediction:** isolate the 1–2 causal variants for Mosaic Variegated Aneuploidy (MVA) out of ~5M variants in a real child's WGS VCF, ranked by a custom `epcr` confidence score.
- **Track 2 — Drug Repurposing:** from the Track 1 gene, characterise the mechanism (LoF/GoF, pathway, downstream effect) and propose FDA-approved repurposing candidates. Judge-scored, qualitative.

`README.md` is the authoritative spec for submission formats and scoring. Read it before touching anything that produces a deliverable.

## Setup (macOS)

```bash
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

brew install openjdk@17                      # snpEff needs a JVM
unzip -q snpEff.zip                          # snpEff.zip is versioned; snpEff/ is not
```

`cyvcf2` builds from source against htslib. If the wheel fails, `brew install htslib` first.

The `java17/` directory and `openjdk17.zip` that may appear in a working copy are a **Windows** JDK, gitignored and useless here — use the Homebrew JVM.

### Data
`data/` is gitignored and never committed: it is a real child's genome plus clinical record. Place the hackathon files at `data/` inside the repo — that is where every notebook and script expects them:

```
data/
├── WGS_EX2312012_HGWCNDSX7.vcf.gz(.tbi)   # patient VCF, GRCh38, ~301 MB
├── common_all_20180418.vcf.gz             # dbSNP Common, ~1.5 GB
└── Challenge_Clinical_Phenotype_1.docx
```

dbSNP Common, if not already present:
```bash
curl -O ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/common_all_20180418.vcf.gz
```

**Never commit VCFs, FASTQs, the clinical `.docx`, or extracted patient identifiers.** `.gitignore` blocks `*.vcf*` as a backstop.

## Commands

No build system or CI. Scripts take `--input`/`--output` and default to the right paths under `data/`.

```bash
# Step 2 — QC + rare-variant prefilter (5M -> ~380k). The long one.
python notebooks/notebook/run_step2_prefilter.py

# Step 3A — SnpEff (below) then MyVariant.info annotation
python notebooks/notebook/track1_step3a_annotate.py --workers 4

# Step 3B — AlphaGenome Atlas regulatory enrichment
python notebooks/notebook/step3b_alphagenome.py

# Notebooks assume notebooks/ as cwd (paths are ../data/...)
cd notebooks && jupyter lab

python -m pytest tests/                      # pytest not in requirements yet; tests/ is empty
```

### snpEff (Step 3A, run before the annotate script)
```bash
java -jar snpEff/snpEff.jar download GRCh38.99          # ~5 GB, once
java -Xmx8g -jar snpEff/snpEff.jar -canon GRCh38.99 \
    data/track1_step2_rare_variants.vcf.gz > snpeff_annotated.vcf
```

`-canon` is required by the procedure: one row per variant, canonical transcript, most severe consequence. It is what keeps the output free of transcript duplicates.

`bcftools`/`tabix` are not used anywhere — all VCF I/O goes through `cyvcf2` or raw `gzip`. `brew install bcftools` if you want them for inspection.

## Pipeline architecture (Track 1)

Each step consumes and emits a **bgzipped VCF**, never a CSV. This is a deliberate contract: new annotations are injected as `INFO` tags so the file stays valid input for the next tool in the chain. Only Step 5 emits CSV.

| Step | Input | Output | Adds |
|---|---|---|---|
| 0 EDA | raw VCF + `.docx` | — | `track1_step0_exp01-data_exploration.ipynb` |
| 1 HPO extraction | `Challenge_Clinical_Phenotype_1.docx` | HPO term list | terms come pre-tabulated in the docx; parse `doc.tables` |
| 2 Pre-filtering | patient VCF (5,012,204 variants) | `data/track1_step2_rare_variants.vcf.gz` | `AB_FILTER` INFO tag |
| 3A snpEff + MyVariant.info | step 2 output | `track1_step3a_snpeff_filter.vcf.gz` | `GENE_SYMBOL`, `CONSEQUENCE`, `PUTATIVE_IMPACT`, `GNOMAD_AF`, `CADD_PHRED`, `REVEL`, `CLINVAR_SIG`, `CLINVAR_REVIEW_STATUS`, `OMIM_PHENOTYPE`, `LOFTOOL_SCORE` |
| 3B AlphaGenome Atlas | step 3A output | `track1_step3b_regulatory_annotated.vcf.gz` | `AVI_SCORE`, `TARGET_GENE`, `REGULATORY_FEATURE` |
| 4 Gene prioritisation | 3A/3B + HPO terms | ranked genes | HPO ↔ gene-phenotype overlap |
| 5 Ranking | candidates | `output/track1/track1_submission_template.csv` | `epcr` |

Division of labour in Step 3A: **snpEff (local) owns `GENE_SYMBOL` / `CONSEQUENCE` / `PUTATIVE_IMPACT`** — computed dynamically so novel variants are covered; **MyVariant.info (REST) owns everything else**, batched at 1000 ids per POST. Steps 3A and 3B are deliberately parallel arms, comparing traditional functional annotation against deep-learning regulatory scoring.

Filtering thresholds are justified in `docs/01-procedure.md` — read the rationale before changing one. In particular `GQ >= 20`, `DP >= 10`, `FILTER` in {`PASS`, missing}, allelic balance 0.30–0.70 for hets, and dbSNP Common subtraction (~1% MAF). The 1% cutoff is a deliberate compromise: `docs/00-references.md` §4 documents a founder-effect variant reaching 3.6% in a subpopulation, so do not tighten frequency filtering without re-reading that section.

### How Step 2 is implemented, and why
`run_step2_prefilter.py` does three **sequential** passes with bounded memory (~250 MB measured), rather than one tabix lookup per variant:

1. patient VCF → temp VCF of intrinsic-QC survivors, collecting `chrom:pos:ref:alt` keys
2. dbSNP read once in binary mode, splitting only the first 5 fields, discarding common keys
3. temp VCF rewritten keeping only the rare variants

Two constraints drove this. Loading dbSNP into a `set` of tuples costs ~287 bytes/variant — about **11 GB** for its ~38M records. And a per-variant tabix query means ~2.6M iterator constructions plus full parsing of dbSNP records whose `INFO` field is never used. Keys are stored **exact, not hashed**: a hash collision would silently drop a variant, and that variant could be the causal one.

Within the QC loop, filters run cheapest-first — `variant.format()` allocates a numpy array per call, so `AD` is not fetched until `GQ`, `DP` and `GT` have rejected what they can.

### Fail-loud contract
The Step 3A and 3B scripts **abort** when the API exhausts retries or `AtlasClient` will not initialise. They used to write an unannotated VCF that looked complete. A half-annotated file is worse than an error when the next step ranks variants on those exact fields.

### Current state
Step 3 has not successfully completed. `snpeff_annotated.vcf` and both `track1_step3*.vcf.gz` in the repo root are empty placeholders from failed runs, and the snpEff GRCh38 database has not been downloaded.

`docs/01-procedure.md` records 591,339 surviving rare variants from Step 2 while `docs/02-findings.md` records 380,292. Re-run Step 2 and reconcile before relying on either figure.

## Working conventions

From `README.md` §"Project Workflow & Guidelines", enforced by team convention rather than tooling:

- **Experiment first, consolidate later.** Exploration lives in `notebooks/` as versioned experiments. Once an approach settles, the code moves into `src/` as reusable modules with matching tests in `tests/`. Both are currently empty; the standalone scripts in `notebooks/notebook/` are the consolidation staging ground.
- **Notebook naming:** `track{N}_step{S}_exp{NN}-{slug}.ipynb`. Git branches mirror the notebook being worked on.
- **Docs are living and must be updated as you go:**
  - `docs/01-procedure.md` — the plan; every step has an **Experiments** table to fill in as work completes.
  - `docs/02-findings.md` — record insights immediately, attributed to the notebook that produced them.
  - `docs/00-references.md` — every paper or external source consulted.
- Docs and README are in English; the pipeline scripts carry Spanish comments. Either is fine — don't rewrite existing ones.

## Submission constraints that shape the code

- **Scoring is F-max plus rank points.** Padding is punished hard: many variants at similarly high `epcr` collapses precision even when the true variant is present. The scoring function must produce *separation*, not just correct ordering — one confident row above a clearly lower-confidence tail beats ten near-tied rows.
- **Row #1 must be the answer.** A compound-heterozygous pair goes in the `_1` and `_2` column groups of a *single* row; homozygous and dominant/de novo cases fill only `_1`.
- Coordinates are **GRCh38 with numeric contigs** (`1`, not `chr1`) in the patient VCF and dbSNP. The submission template's example rows use `chr1`, and the Step 3 scripts add the prefix for API ids — check which convention applies before formatting output.
- Submission filename: `<username>_<model-name>.csv`. 6 Track 1 submissions allowed, 3 for Track 2 (only the latest Track 2 entry is reviewed).
- Track 2 requires an AI-usage disclosure (provider, tier, data-handling terms) in the report.
