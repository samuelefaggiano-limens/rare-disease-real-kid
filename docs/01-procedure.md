# High-Level Procedure: MVA Hackathon 2026

This document outlines the exhaustive, end-to-end bioinformatics and computational biology pipeline required to tackle both Track 1 (Variant Prediction) and Track 2 (Drug Repurposing) of the Rare Disease, Real Kid Hackathon.

## Track 1: Variant Prediction Pipeline
**Goal:** Isolate 1-2 causal genetic variants from ~5 million background variants present in the patient's Whole Genome Sequencing (WGS) VCF file.

### Step 0: Exploratory Data Analysis (EDA)
*   **Input:** All provided data files in `data/` (Clinical documents and VCF).
*   **Action:** Perform initial exploration to understand data structures, quality, and basic statistics before setting up the main pipeline.
*   **Tools:** Python, pandas, Jupyter Notebooks.
*   **VCF Column Definitions (Understanding the Dataset):**
    *   `CHROM`: Chromosome where the variant is located (e.g., '1', 'X').
    *   `POS`: The 1-based base position of the variant on the chromosome.
    *   `ID`: Variant identifier, usually a dbSNP rsID (e.g., 'rs123456'). Represented as '.' if unknown.
    *   `REF`: Reference allele (the DNA sequence found in the reference genome).
    *   `ALT`: Alternate allele (the DNA sequence found in the patient's genome).
    *   `QUAL`: Phred-scaled quality score for the assertion that an alternate allele exists.
    *   `FILTER`: Indicates if the variant passed quality filters ('PASS') or failed specific checks (e.g., 'LowQual', 'MQ40').
    *   `INFO_DP`: Combined read depth of coverage across all samples (from the INFO field).
    *   `INFO_MQ`: Root Mean Square (RMS) mapping quality of all reads spanning the variant.
    *   `INFO_FS`: Phred-scaled p-value using Fisher's exact test to detect strand bias (whether reads inappropriately favor one DNA strand over the other).
    *   `FORMAT_GT`: Genotype for the patient (e.g., '0/1' for heterozygous, '1/1' for homozygous alternate).
    *   `FORMAT_DP`: Read depth specifically for this patient's sample.
    *   `FORMAT_GQ`: Genotype quality (confidence in the specific GT call).
    *   `FORMAT_AD_REF`: Allelic depth (number of reads) specifically supporting the reference allele.
    *   `FORMAT_AD_ALT`: Allelic depth (number of reads) specifically supporting the alternate allele.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step0_exp01-data_exploration.ipynb` | Initial exploration of VCF and Clinical data formats, missingness, and general stats. | Done |

### Step 1: Phenotype Extraction (HPO Mapping)
*   **Input:** `data/Challenge_Clinical_Phenotype_1.docx`
*   **Action:** Read the clinical document and extract standardized medical symptoms. We will translate free-text symptoms (e.g., "small head") into standard **Human Phenotype Ontology (HPO)** terms (e.g., HP:0000252 - Microcephaly).
*   **Tools:** Manual review, or NLP tools like Doc2HPO / ClinPhen.
*   **Output:** A list of HPO terms representing the patient's exact clinical presentation.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step1_exp01-hpo_extraction.ipynb` | Extracts the HPO terms provided in the `Challenge_Clinical_Phenotype_1.docx` file into a list of IDs. | Done |

### Step 2: VCF Pre-Filtering
*   **Input:** `data/WGS_EX2312012_HGWCNDSX7.vcf.gz` (Unannotated VCF, 5,012,204 variants)
*   **Action:** Reduce the massive 5M row VCF file to a manageable size using intrinsic quality metrics and a fast extrinsic common-variant subtraction, prior to API-based annotation. This single, rapid step will safely strip away about 90-95% of the VCF leaving only rare/novel variants across the entire genome (coding and non-coding).
*   **Filter 1: Intrinsic QC Filters (Zero external downloads)**
    *   **What it does:** Relies on the sequencing machine's own confidence scores present in the VCF file.
    *   **Implementation details (based on `track1_step0_exp01-data_exploration.ipynb`):**
        *   **The PASS Filter:** The `FILTER` column contains flags for unreliable variants (e.g., `LowQual`, `MQ40`). We keep only rows where `FILTER` is missing (represented as `NaN` or `.` in pandas/cyvcf2) or `PASS`.
        *   **Genotype Quality (GQ) & Sample-Level Depth (DP):** Instead of using aggregate metrics like `INFO_DP` or `QUAL` which can be unreliable in difficult mapping regions, we prioritize sample-level confidence. We use `FORMAT_GQ` to measure the confidence that the specific 0/1 or 1/1 call for this child is correct, and `FORMAT_DP` to ensure the sequencer read this specific coordinate enough times *in this patient's sample*. We keep only rows where `FORMAT_GQ >= 20` (99% confidence) and `FORMAT_DP >= 10`.
        *   **Output Format:** The output of this stage must be saved as a standard VCF (e.g., `track1_step2_intrinsic_qc.vcf.gz`), injecting any custom metrics (like Allelic Balance) as VCF `INFO` tags. Outputting a CSV here breaks compatibility with the bioinformatics tools required for Filter 2 and Step 3.
*   **Filter 2: The "Common Variant" Subtraction (Lightweight Extrinsic)**
    *   **What it does:** Since MVA is an ultra-rare pediatric disease, the causal mutation cannot be common in the healthy public. We use a small "blacklist" database to perform a rapid subtraction.
    *   **Implementation details:**
        *   **Dataset (The Lightweight Blacklist):** We use **dbSNP Common (GRCh38)** (`common_all_20180418.vcf.gz`), which is a curated ~1.5 GB single file containing known common human variants (MAF >= 1% in at least one major population). This avoids the extreme computational overhead of streaming hundreds of gigabytes of gnomAD VCFs.
        *   **Chromosome Harmonization:** dbSNP Common natively uses numeric contigs (e.g., `1`, `2`), perfectly matching our patient's GRCh38 VCF.
        *   **Threshold:** dbSNP Common implicitly acts as a ~1% MAF filter. While slightly higher than a strict 0.5% cutoff, literature shows (e.g., Cho et al. 2025 on MOCS2) that recessive founder effects can push pathogenic alleles >3% in subpopulations, meaning a 1% threshold is actually a very safe boundary to prevent accidentally deleting a causal variant.
        *   **Tools:** We implement an *allele-aware* subtraction using `cyvcf2` directly in Python. By iterating over our surviving rare variants and performing a fast exact-match index lookup into the dbSNP VCF, we safely drop any variant that exactly matches both coordinates and alleles of a common variant, leaving us with our final "Rare" VCF.
*   **Experiments:**

| Notebook | Description | Results                                                                                                                                                                                 |
| :--- | :--- |:----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `track1_step2_exp01-vcf_prefiltering.ipynb` | Applies Intrinsic QC & AB filters, followed by a Common Variant subtraction (MAF >= 1%) using dbSNP GRCh38, one tabix lookup per variant. | Processed 5,012,204 variants -> 4,428,273 survived intrinsic QC. Subtraction dropped 3,836,934 common variants. Final result: 591,339 ultra-rare variants saved to `track1_step2_rare_variants.vcf.gz` (405 s). |
| `notebooks/notebook/run_step2_prefilter.py` | Same thresholds as exp01, re-implemented as three sequential passes (QC -> single scan of dbSNP -> rewrite) with bounded memory. Re-run on macOS (Apple Silicon) on 2026-09-17 to reconcile the 591,339 vs 380,292 (`02-findings.md`) discrepancy. | Reproduces exp01: 4,428,273 QC survivors, 3,836,941 dropped as common, **591,339 rare variants** (11.8%) in 1.4 min. 380,292 was wrong and has been corrected in `02-findings.md`. |

# Revised Procedure: Steps 3–5 (Track 1)

Step 3 is implemented as **one single notebook** (`track1_step3_exp01-annotation_and_filtering.ipynb`). It reads `track1_step2_rare_variants.vcf.gz` once, runs a cascade of annotation and filtering stages against it in memory, and writes **one output file**: `track1_step3_annotation_and_filtering.vcf.gz`. Nothing is saved to disk between the sub-stages below — they are steps inside the same notebook, not separate pipeline runs.

---

### Step 3: VCF Annotation & Filtering

*   **Input:** Pre-filtered VCF (`track1_step2_rare_variants.vcf.gz`, ~591,339 variants).
*   **Action:** Run a sequential cascade — a known-pathogenic lookup, local functional annotation, and population/clinical metadata enrichment — where each stage narrows or enriches the variant set for the next.
*   **Output of the whole step:** `track1_step3_annotation_and_filtering.vcf.gz`, ready for Step 4.

#### Stage 3.1: Known-Pathogenic Fast Pass

*   **What this does:** Before doing any of our own analysis, we check whether any of our 591,339 candidate variants have already been studied and classified by other clinical laboratories, using the public ClinVar database. Think of it as checking "has anyone already solved this specific line of the genome before?" before we spend hours building our own scoring system.
*   **Input:** `track1_step2_rare_variants.vcf.gz`.
*   **Implementation:**
    *   Download the ClinVar GRCh38 VCF (~100 MB) and its index.
    *   Normalise both our variant file and the ClinVar file with `bcftools norm -f GRCh38.fa -m -any -c w` so that the same mutation is always written the same way in both files (see the note on normalisation below — this is easy to skip and easy to regret skipping).
    *   Intersect the two files (`bcftools isec -n=2 -w1`) and keep only matches flagged "Pathogenic" or "Likely pathogenic."
    *   Record the ClinVar confidence level (the "star rating") for each hit.
*   **Why it matters:** this costs a few minutes of compute and no API calls. Since a team of medical geneticists already diagnosed this child, there is a real chance the causal variant is already sitting in a public clinical database. If we get a hit in a known MVA gene, we may have our answer immediately, and everything after becomes confirmation rather than discovery. If we get nothing, that is useful information too — it tells us the causal variant is probably novel and never recorded anywhere, so we must not treat "absent from ClinVar" as evidence that a variant is harmless later in the pipeline.

**A note on "normalisation" (relevant to every stage in this step):** the same genetic mutation can be written in more than one equivalent way in a VCF file — for example, a small deletion can be described starting from slightly different positions depending on the software that produced the file. If we don't force every file we use (our data, ClinVar, SnpEff's output) to write mutations the same standard way, two files can describe the identical mutation with different text, and our matching logic will conclude they are different variants — a silent, invisible error. Running `bcftools norm` once, up front, prevents this for every downstream comparison in Step 3.

#### Stage 3.2: Local Functional Annotation (SnpEff)

*   **What this does:** For every surviving variant, we determine which gene it falls inside, what it does to that gene's protein (does it change one amino acid? does it break the protein entirely?), and how severe that change is predicted to be. This is like translating a raw DNA typo into "this breaks gene X in a serious way" versus "this makes a harmless spelling change."
*   **Input:** `track1_step2_rare_variants.vcf.gz` (normalised).
*   **New annotations added:**
    *   `GENE_SYMBOL` — the gene the variant falls in (e.g. *BUB1B*).
    *   `CONSEQUENCE` — the type of protein-level change (e.g. missense, frameshift, stop-gain).
    *   `PUTATIVE_IMPACT` — a severity bucket: `HIGH`, `MODERATE`, `LOW`, or `MODIFIER`.
    *   `TRANSCRIPT_ID`, `HGVS_C`, `HGVS_P` — the exact transcript used for the call, and standard coding/protein notation (needed later to cross-check against ClinVar and to write a readable `notes` field in the final submission).
    *   `EXON_RANK` / `EXON_TOTAL` — which exon (out of how many) the variant falls in. This matters because a "stop" mutation near the beginning of a gene usually destroys the protein completely, while the same kind of mutation in the very last exon often does much less damage. This distinction becomes important for the drug-mechanism reasoning in Track 2.
*   **Tool:** SnpEff (pin the exact version and GRCh38 database build used, and record both in the methods report).
*   **Transcript choice — using MANE Select, not the longest transcript:** genes can have multiple "versions" (transcripts) of similar but not identical sequence. We must pick one to base our gene/consequence calls on. We use the **MANE Select** transcript — the transcript that NCBI and EBI have jointly designated as the clinical reference for each gene — rather than simply the longest one. This matters because clinical databases like ClinVar report variants against the MANE transcript; if we use a different one, our variant descriptions won't match theirs and cross-referencing will silently fail.
*   **Filtering — what we keep and why:**
    *   Keep only variants where `PUTATIVE_IMPACT` is `HIGH` or `MODERATE`. In plain terms: we discard variants predicted to have little or no effect on the protein (`LOW`/`MODIFIER`), and keep only those predicted to meaningfully change or break the protein. This is expected to reduce ~591,000 variants down to roughly 10,000–20,000.
    *   **Exception — always keep `splice_region_variant` regardless of its severity label.** SnpEff labels these as `LOW` impact by default, which would normally exclude them from the filter above. We override that, because these variants sit in the small stretch of DNA right next to where a gene's coding instructions are "spliced" together, and mutations there are a well-documented cause of serious inherited disease even though the tool underrates their severity. Skipping this exception is one of the easiest ways to accidentally throw away the correct answer.
*   **Output of this stage:** the surviving ~10,000–20,000 variants, now carrying the fields above, held in memory for the next stage.

#### Stage 3.3: Population Frequency & Clinical Metadata (MyVariant.info)

*   **What this does:** For each of the ~10,000–20,000 surviving variants, we ask a public annotation database three questions: "How common is this variant in healthy people?", "Has a computer model predicted how damaging it is?", and "Has any clinical lab already classified it?" A variant that is common in healthy adults essentially cannot be the cause of an ultra-rare disease, so this stage is largely about elimination by population frequency.
*   **Input:** the surviving variants from Stage 3.2.
*   **New annotations added:**
    *   `GNOMAD_AF` — how frequently this exact variant appears in a reference database of roughly 800,000 people (gnomAD). The lower, the rarer, and the more plausible as a cause of an ultra-rare disease.
    *   `GNOMAD_AF_POPMAX` — the highest frequency seen in any single population group. A variant that looks rare overall but is common in one specific ancestry group is more likely a normal population-specific variant than a disease cause.
    *   `GNOMAD_NHOMALT` — the number of healthy people in that same reference database who carry **two copies** of this variant (one from each parent). This is the single most decisive number in this stage for a recessive disease like MVA: if healthy adults exist who have two copies of this variant and are not sick, it cannot be the two-copy cause of this child's disease, no matter how "damaging" it looks by other measures.
    *   `REVEL` and `ALPHAMISSENSE_SCORE` / `ALPHAMISSENSE_PRED` — two independent computer predictions of how damaging a single-amino-acid change (missense variant) is likely to be. These only apply to missense variants; they will be empty for variants that delete or frameshift the gene, because "how damaging is this amino-acid substitution" isn't a meaningful question for a broken gene.
    *   `CLINVAR_SIG` and `CLINVAR_REVIEW_STATUS` — whether any clinical lab has already classified this exact variant (e.g. "Pathogenic," "Benign," "Uncertain"), and how much confidence to place in that classification (a rating from "one anonymous submission" up to "reviewed by an expert panel"). A single variant can have several conflicting classifications from different submitters, so we keep the highest-confidence one and note if there is disagreement.
    *   `GENE_LOEUF` / `GENE_PLI` — a gene-level (not variant-level) measure of how strongly that gene resists loss-of-function mutations across the human population. A gene where such mutations are almost never seen in healthy people is a gene where a loss-of-function mutation is more likely to cause disease.
*   **Practical details:**
    *   Query the MyVariant.info API in batches of 1,000 variants, with automatic retry if a request fails.
    *   Explicitly request the GRCh38 (`hg38`) version of the database — the service defaults to an older genome build (`hg19`) unless told otherwise, which would silently return the wrong information.
    *   Some fields (ClinVar classifications, REVEL/AlphaMissense scores) can come back as a list with more than one entry, because they are reported per-transcript or per-submission. We reduce each list to a single representative value (the most severe score, or the highest-confidence classification).
*   **Expected gaps in the data (plan for these now, don't be surprised later):**
    *   REVEL and AlphaMissense will be empty for every frameshift, deletion, insertion and splice-region variant — this is expected, not a data quality problem.
    *   Most variants will have no ClinVar entry at all. Absence of a ClinVar record means "nobody has studied this variant yet," not "this variant is harmless." These two must be scored differently later, in Step 5.
*   **Output of this stage — and of Step 3 as a whole:** `track1_step3_annotation_and_filtering.vcf.gz`, the final annotated candidate set, written to disk for use in Step 4.

#### Implementation status (branch `track1_opt01-macos-portable`, run 2026-09-17)

What was actually run diverges from the revised plan above in these points — to be reconciled by the team:

*   Implemented as **scripts + one results notebook** (`notebooks/notebook/track1_step3a_annotate.py`, `step3b_alphagenome.py`, `track1_step3_exp01-vcf_annotating.ipynb`), writing intermediate VCFs to disk (`snpeff_annotated.vcf`, `track1_step3a_snpeff_filter.vcf.gz`, `track1_step3b_regulatory_annotated.vcf.gz`), not a single in-memory notebook.
*   **Stage 3.1 (ClinVar fast pass) and `bcftools norm` were not run**; ClinVar is obtained per variant through MyVariant.info instead.
*   snpEff uses `-canon` (canonical transcript), **not MANE Select**; `TRANSCRIPT_ID`/`HGVS_C`/`HGVS_P`/`EXON_RANK` are available inside the raw `ANN=` field but not written as separate tags. `splice_region_variant` (LOW) is **not** kept.
*   MyVariant fields obtained: `GNOMAD_AF`, `GNOMAD_EXOME_AF`, `CADD_PHRED` (via dbNSFP), `REVEL`, `CLINVAR_SIG`, `CLINVAR_REVIEW_STATUS`, `OMIM_PHENOTYPE`. **Not yet obtained:** `GNOMAD_AF_POPMAX`, `GNOMAD_NHOMALT`, `ALPHAMISSENSE_*`, `GENE_LOEUF`/`GENE_PLI`. `assembly=hg38` confirmed as mandatory (see `02-findings.md`).
*   An additional **Step 3B (AlphaGenome Atlas AVI)** was run: `AVI_SCORE`, `AVI_QUANTILE`, `AVI_PHRED`, `REGULATORY_FEATURE`, `AVI_STATUS` (SNVs only).
*   Real outcome: 591,339 -> **1,526** HIGH/MODERATE (not 10–20k); 1,212 of them have gnomAD AF > 1 % because dbSNP Common misses most common coding variants.

*   **Implementation notes (`notebooks/notebook/step3b_alphagenome.py`, alphagenome 0.9.0, verified against the live API on 2026-09-17):**
    *   The client must be created with `atlas.create(api_key)`; a bare `AtlasClient()` has no gRPC stub. The key is read from `ALPHAGENOME_API_KEY` (env or `~/.env`).
    *   Scorers used: `AVI_SCORE` (raw score, plus a `quantiles` layer) and `AVI_SCORE_FEATURE_IMPORTANCE` (18 modalities). Tags written: `AVI_SCORE`, `AVI_QUANTILE`, `AVI_PHRED` (= -10 log10(1 - quantile); >= 20 is the top 1 % of genome-wide SNVs, >= 30 the top 0.1 %), `REGULATORY_FEATURE` (modality with the largest attribution, e.g. `MERGED_SPLICING`, `ALPHAMISSENSE`, `MAX_ABS_RNA_SEQ`), `AVI_STATUS`.
    *   **`TARGET_GENE` is not returned by the AVI scorers** and is therefore not written; the gene remains snpEff's `GENE_SYMBOL` from Step 3A.
    *   The server validates the REF allele against GRCh38 and one bad variant fails the whole batch; the script falls back to single-variant queries and records `AVI_STATUS=ref_mismatch` for the offenders instead of aborting. Infrastructure errors still abort after 5 retries (fail-loud contract).
    *   Downstream use: Step 5 feeds `AVI_PHRED` into the pathogenicity component on the same Phred scale as CADD.

**Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| snpEff `GRCh38.99` database (local build) | The versioned `snpEff.zip` is SnpEff 4.5covid19 (2020-04-15). Its `download GRCh38.99` fails: SourceForge has no `v4_5` databases (only `v4_3` GRCh38.86/92) and the newer Azure repository host no longer resolves (NXDOMAIN). The database is therefore **built locally** from Ensembl release 99 (GTF + primary-assembly FASTA + CDS + protein FASTA) with `snpEff build -gtf22 GRCh38.99`, which is the same Ensembl 99 gene set the prebuilt database was derived from. | Built 2026-09-17 in 1.5 min (`snpEffectPredictor.bin`, 117 MB). The build's protein check reports "Not found" for every transcript because Ensembl FASTA ids carry versions (`ENST…​.2`) and the GTF ids do not; it is a validation-only step. Sanity check passed: HBB rs334 chr11:5227002 T>A -> `missense_variant p.Glu7Val`, and a BUB1B intronic control -> `intron_variant c.180-62G>A`. |
| snpEff `-canon` annotation | `java -Xmx8g -jar snpEff/snpEff.jar -canon -noStats GRCh38.99 data/track1_step2_rare_variants.vcf.gz > snpeff_annotated.vcf` | 591,339 / 591,339 variants annotated in 15 s. Impact of the canonical annotation: 588,497 MODIFIER, 1,324 MODERATE, 1,316 LOW, **202 HIGH**. |
| `notebooks/notebook/track1_step3a_annotate.py` | Reads the snpEff `-canon` output, keeps HIGH/MODERATE, adds MyVariant.info fields (batched, threaded, fail-loud). **Fixed 2026-09-17:** queries now send `assembly=hg38` (the default hg19 lookup returned 97 % "notfound" and 3 % wrong-variant hits), indels are sent as valid HGVS `del`/`ins`/`delins`, CADD comes from `dbnsfp.cadd.phred` (the `cadd` dataset is hg19-only), OMIM from ClinVar conditions; LoFtool is not in MyVariant and is no longer emitted. | 1,526 HIGH/MODERATE variants (202 HIGH, 1,324 MODERATE; 1,219 SNVs, 307 indels) in 0.3 min. Coverage: MyVariant hit 1,525, `GNOMAD_AF` 1,481, `CADD_PHRED` 1,110, `REVEL` 864, `CLINVAR_SIG` 312, `OMIM_PHENOTYPE` 112. **1,212 of the 1,526 have gnomAD AF > 1 %** (see `02-findings.md`). |
| `notebooks/notebook/step3b_alphagenome.py` | AlphaGenome Atlas AVI enrichment (`AVI_SCORE`, `AVI_QUANTILE`, `AVI_PHRED`, `REGULATORY_FEATURE`, `AVI_STATUS`). Indels are marked `indel_unsupported` without a call (Atlas returns `NotImplementedError` for them); per-variant fallback on server-side validation errors. | 1,526 written, **1,220 scored** (all SNVs), 306 `indel_unsupported`, 0.9 min. `AVI_PHRED` >= 20 for 257 variants, >= 30 for 13. Top modality: `CACTUS_241_WAY` 562, `ALPHAMISSENSE` 509, `PHASTCONS_470_WAY` 73, `MERGED_SPLICING` 45, `PROTEIN_TERMINATION` 18. Output: `track1_step3b_regulatory_annotated.vcf.gz`. |
| `notebooks/notebook/track1_step3_exp01-vcf_annotating.ipynb` | Executed walkthrough of Step 3 for the team: what each tool adds, how to reproduce, and summaries computed from the three output VCFs (impact distribution, MyVariant coverage, gnomAD AF of the "rare" set, ClinVar classes, AVI status/modalities/Phred, gene-level MVA panel check). Reads the gitignored outputs; makes no API calls and prints no patient coordinates. | Impact: 202 HIGH / 1,324 MODERATE / 1,316 LOW / 588,497 MODIFIER. 1,212 of 1,526 HIGH/MODERATE have gnomAD AF > 1 %. AVI: 1,220 scored, 257 with Phred >= 20. Panel: two heterozygous *BUB1B* variants (stop-gain, ClinVar P/LP, OMIM 257300; missense absent from gnomAD), none in *CEP57*/*TRIP13*/*MAD2L1BP*. |

---

### Step 4: Gene Prioritization & Phenotype Matching

*   **Input:**
    *   The patient's extracted clinical HPO terms (from Step 1).
    *   `track1_step3_annotation_and_filtering.vcf.gz`.
*   **Action:** Compare the genes disrupted by our candidate variants against the child's actual documented symptoms, and boost variants in genes that are already known to cause this specific disease.
*   **Implementation Details:**
    *   **Gene mapping:** each variant is linked to a gene using SnpEff's `GENE_SYMBOL` from Stage 3.2.
    *   **Matching symptoms to genes (semantic similarity):** matching the child's symptoms to a gene's known symptoms cannot rely on exact text matching, because doctors and databases describe the same underlying problem in different words (e.g. "seizures" versus the broader "abnormality of the nervous system"). We instead use the **Human Phenotype Ontology (HPO)**, a structured map of how medical terms relate to one another, and measure how *close* the child's symptoms are to each candidate gene's known symptom profile, using standard similarity algorithms (Resnik, Lin) via the `pyhpo` library. Each variant receives a `PHENOTYPE_SIMILARITY_SCORE` based on this closeness.
        *   We also record a `PHENOTYPE_SCORE_IS_NULL` flag for genes that have **no** documented symptom profile at all in the database. This distinction matters: a gene with no documented symptoms is simply unstudied, which is different from a gene whose documented symptoms clearly don't match this child. Both would otherwise score as "0," which would unfairly bury variants in unstudied genes — and this disease is rare enough that an unstudied gene is a real possibility.
    *   **Reference data:** we load the HPO project's `genes_to_phenotype.txt` file locally (rather than querying a live web API), and record exactly which release version of this file we used, so the analysis can be reproduced later.
    *   **The MVA gene panel (targeted boost):** MVA is primarily caused by defects in a cellular quality-control system called the Spindle Assembly Checkpoint, plus related centrosome biology. We compile a list of relevant genes and give variants in those genes a boost — but not all genes on the list are equally certain, so we use two confidence tiers rather than one:
        *   **Tier 1 — genes with confirmed causal roles in MVA specifically:** *BUB1B*, *CEP57*, *TRIP13*, *CENATAC*/*CCDC84*. These receive a strong boost.
        *   **Tier 2 — genes involved in the same general biological system, but not established causes of MVA specifically:** *MAD2L1, CDC20, CENPE, SGO1, AURKB, PLK1, MAD1L1, BUB1, BUB3, TTK, ESPL1, PTTG1, NDC80, CENPF*. These receive a smaller boost. Splitting the list this way matters because some of these genes are already known to cause *other*, different diseases when disrupted — a large boost applied indiscriminately to all of them risks pushing an unrelated finding above the real answer.
        *   This gene panel, with a short note on why each gene is included, is kept in its own file (`data/mva_panel.tsv`) rather than buried in code, so it is easy to review and update.
    *   **Identifying compound-heterozygous candidates:** MVA is usually recessive, meaning it typically requires *both* copies of a gene to be broken — either the same broken copy inherited from both parents, or two different broken copies, one from each parent (a "compound heterozygous" case). We group our surviving variants by gene and flag any gene carrying two or more different candidate variants, since this is exactly the pattern the final submission format is designed to capture. Because we do not have the parents' DNA, we cannot directly prove that two such variants came from different parents — we can only note that a gene carries two plausible broken copies, and we state this limitation clearly rather than overstating our confidence.
*   **Tools:** Python (`pyhpo`, `networkx`, `pandas`), the local `genes_to_phenotype.txt` dataset, `data/mva_panel.tsv`.
*   **Output:** a master candidate table (CSV), enriched with `PHENOTYPE_SIMILARITY_SCORE`, `PHENOTYPE_SCORE_IS_NULL`, `MVA_PANEL_TIER`, and `COMPOUND_HET_CANDIDATE`, ready for final ranking in Step 5.

#### Implementation status (branch `track1_opt01-macos-portable`)

Implemented in `notebooks/notebook/track1_step4_prioritize.py` + `src/hpo.py` (tests in `tests/test_hpo.py`), **not yet run on the real data** by team decision. Divergences from the plan above: Resnik similarity is implemented in `src/hpo.py` (no `pyhpo` dependency); the panel is a single tier in code plus a curated SAC/centrosome list (no `data/mva_panel.tsv`, no Tier 1/Tier 2 boost yet; *CENATAC/CCDC84* not included, *MAD2L1BP* included); the null-phenotype flag exists as `has_hpo_annotations`. Outputs: `output/track1/track1_step4_gene_ranking.csv` and `data/track1_step4_candidates.tsv`.

*   **Implementation (`notebooks/notebook/track1_step4_prioritize.py`, logic in `src/hpo.py`):**
    *   **Gene-phenotype source:** HPO's `genes_to_phenotype.txt` + `hp.obo` (release 2026-09-01), downloaded to `data/`. 5,277 genes carry annotations. MVA itself is annotated (OMIM:257300 MVA1 -> *BUB1B*, OMIM:614114 MVA2 -> *CEP57*, OMIM:617598 MVA3 -> *TRIP13*, ORPHA:1052), *MAD2L1BP* has no HPO annotations at all.
    *   **Why semantic similarity and not exact overlap:** of the patient's 8 terms, *BUB1B* and *CEP57* match only 3 exactly (HP:0001518, HP:0002859, HP:0004322), *TRIP13* 2, *MAD2L1BP* 0. Diseases are annotated with siblings/parents of the patient's terms (e.g. Growth delay vs. Failure to thrive), so exact overlap under-scores the true genes. We use Resnik similarity: gene annotations are propagated to ontology ancestors, each term gets an information content `IC = -ln(genes annotated / all genes)`, and a patient term scores `IC(most informative common ancestor) / IC(term)` in [0, 1]. Gene score = mean over the patient's terms.
    *   **Inheritance:** per gene, `biallelic_evidence` = homozygous-alt variant or >= 2 heterozygous variants (compound-het candidates). MVA is autosomal recessive in every known gene.
    *   **Targeted checks:** the 4-gene MVA panel is always printed with every surviving variant, plus a curated SAC/kinetochore/centrosome list.
    *   **Outputs:** `output/track1/track1_step4_gene_ranking.csv` (gene-level, no coordinates) and `data/track1_step4_candidates.tsv` (variant-level; contains patient coordinates, so it stays in the gitignored `data/`).

**Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `notebooks/notebook/track1_step4_prioritize.py` | HPO Resnik similarity per gene + biallelic evidence + MVA panel / SAC check. Unit tests in `tests/test_hpo.py`. | Smoke-tested on a synthetic Step 3A VCF with the real ontology: *BUB1B*/*CEP57* score 0.632, *TRIP13* 0.558, a random gene (*TTN*) 0.438, an unannotated gene 0. PENDING on the real Step 3A output. |

---

### Step 5: Ranking and Submission

*   **Action:** Assign an Estimated Probability of Causal Relationship (`epcr`), a number between 0 and 1, to each candidate variant, representing how confident our pipeline is that this specific variant is the true cause of the disease.
*   **How we calculate it:** a transparent, hand-weighted scoring formula (not a trained machine-learning model — we only have one true answer in this entire project, and we don't know it, so there is no way to properly train or validate a model here; a manually weighted formula is also far easier to explain and double-check).
    *   We score **two separate types of variant differently**, because the evidence available for each type is completely different:
        *   **Loss-of-function variants** (frameshift, stop-gain, and similar): scored mainly on `PUTATIVE_IMPACT = HIGH`, whether the affected exon position means the broken protein is likely eliminated by the cell's own quality control (nonsense-mediated decay), and how strongly the gene resists this kind of damage in the general population (`GENE_LOEUF`). REVEL and AlphaMissense are not used here, because those scores don't apply to this type of variant.
        *   **Missense variants** (single amino-acid changes): scored mainly on REVEL and AlphaMissense.
    *   Both branches also factor in: extreme rarity in the general population and zero healthy people carrying two copies (`GNOMAD_AF`, `GNOMAD_NHOMALT`); a genotype pattern consistent with the expected inheritance pattern; a strong symptom match (`PHENOTYPE_SIMILARITY_SCORE`); membership in the Tier 1 MVA gene list; and any existing clinical classification, weighted by how confident that classification is.
*   **Handling missing information correctly:** several fields will legitimately be blank for many variants, and the scoring formula must treat "blank" correctly rather than as "bad":
    *   A missing REVEL/AlphaMissense score on a frameshift variant is expected — that variant is simply scored using the loss-of-function branch instead.
    *   No ClinVar entry means "unstudied," not "harmless" — it should not lower a variant's score, since the true cause of this child's disease may well be a variant nobody has ever recorded before.
    *   A missing phenotype-similarity score (flagged by `PHENOTYPE_SCORE_IS_NULL`) should fall back on the gene-constraint information instead of being treated as "no symptom match."
*   **Designing the submission for how it will be scored:** the competition rewards a *confident, well-separated* ranking far more than a long list of guesses. As the challenge describes it, submitting one very high-confidence correct guess alongside a handful of clearly low-confidence wrong guesses scores near-perfectly, while submitting many guesses that are all rated similarly highly scores poorly — even if the correct answer is technically included. In practical terms: submit a short list (roughly 5–10 rows, not 50), with a clear score gap between our top pick(s) and everything else, rather than a smooth spread of similar-looking scores.
*   **Building each row of the submission file:**
    *   A variant found as two identical copies (homozygous) → goes in the `_1` columns alone, `_2` left blank.
    *   A compound-heterozygous pair identified in Step 4 → both variants go on the **same row**, one in the `_1` columns and one in the `_2` columns, sharing a single `epcr` score. A proposed pair must never be split across two separate rows — the scoring system treats each row as one complete hypothesis, and splitting a pair turns one potentially fully-correct guess into two only-half-credit guesses.
    *   A single, extremely severe variant proposed as a dominant/spontaneous cause → `_1` columns alone.
    *   The `notes` field for each row should briefly state the gene, the exact mutation, what it does to the protein, and the one-line reason we believe it's relevant.
*   **Using our six submission attempts wisely:**
    *   Submit the Step 3.1 ClinVar fast-pass result first, early on, as a low-cost baseline.
    *   Use the next few submissions for the full pipeline's output as it is refined, ideally changing one meaningful thing at a time so we can tell what changed the leaderboard score.
    *   Keep one or two submissions in reserve for later refinements once the full pipeline's early results are in.
*   **Output:** the ranked candidate list written to `output/track1/track1_submission_template.csv`, named `<username>_<model-name>.csv`.

#### Implementation status (branch `track1_opt01-macos-portable`)

Implemented in `notebooks/notebook/track1_step5_submission.py` + `src/epcr.py` (tests in `tests/test_epcr.py`), **not yet run on the real data**. Divergences from the plan above: a single weighted geometric mean (phenotype 0.35, pathogenicity 0.30, inheritance 0.20, frequency 0.15) rather than separate LoF/missense branches; no `GNOMAD_NHOMALT`/`LOEUF` inputs yet; `AVI_PHRED` feeds pathogenicity alongside CADD/REVEL; gnomAD AF > 1 % is discarded (`--max-af`); separation is enforced with a geometric decay (`--tail-ratio 0.6`, default 5 rows). The Step 3.1 ClinVar-baseline submission has not been made.

*   **Implementation (`notebooks/notebook/track1_step5_submission.py`, logic in `src/epcr.py`):**
    *   **Leaderboard contract (from the Space's public `evaluation.py` / `tabs/submit_track1.py`):** `proband_id` must be `PROBAND01`; chromosomes carry the `chr` prefix (`chr15`); max 10 rows; `epcr` in (0, 1]; rank points 100/50/25/10 for rank 1/<=3/<=5/<=10 (half credit if only one variant of a compound-het pair is recovered); **F-max is computed per variant**, sweeping every distinct `epcr` as a threshold.
    *   **Components (each in [0, 1]):** phenotype (Step 4 score), pathogenicity (snpEff impact lifted/tempered by CADD_PHRED, REVEL, ClinVar), inheritance (1.0 for hom-alt or compound-het candidate, 0.45 for a lone het), frequency (linear in -log10 gnomAD AF between 1e-2 and 1e-5; missing = 1.0). Combined as a **weighted geometric mean** (0.35 / 0.30 / 0.20 / 0.15) so a weak axis costs multiplicatively; x1.2 bonus for MVA panel / MVA-annotated genes.
    *   **Row construction:** genes with >= 2 hets emit one compound-het row (two best variants, `_1` + `_2`); hom-alt variants emit a single row; lone hets emit a single penalised row.
    *   **Separation, by design:** keep `--max-rows` (default 5) and force `epcr_k <= 0.6 * epcr_{k-1}`. F-max punishes near-ties (10 rows at 0.9 cap precision at 0.1 even if the answer is among them); a forced geometric decay turns "equally good" into an explicit ranking. The script prints a what-if table (leaderboard result if row k were the truth) so the cost of every extra row is visible before submitting.

**Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `notebooks/notebook/track1_step5_submission.py` | epcr = weighted geometric mean + MVA bonus; compound-het pairing; geometric separation; local re-implementation of the leaderboard metrics. Unit tests in `tests/test_epcr.py`. | Smoke test (synthetic candidates): BUB1B compound-het row epcr 0.946, CEP57 hom row 0.567, tail <= 0.34. What-if: row 1 true -> 100 pts / F-max 1.0; row 2 true -> 50 pts / F-max 0.5. PENDING on real data. |

---

## Track 2: Drug Repurposing Pipeline
**Goal:** Propose existing, FDA-approved medications that could bypass or correct the biological mechanism disrupted by the variant found in Track 1.

### Step 1: Mechanism Characterization
*   **Input:** The causal gene and variant identified in Track 1.
*   **Action:** Conduct a deep literature review to understand the normal function of the gene and how the mutation breaks it.
*   **Key Questions to Answer:**
    *   Is it a **Loss of Function (LoF)** (e.g., missing protein, unstable protein) or **Gain of Function (GoF)** (toxic overactivity)?
    *   What is the specific pathway? (e.g., Kinetochore attachment, chromosome segregation).
    *   What are the downstream consequences? (e.g., premature anaphase, cell death, aneuploidy).
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

### Step 2: Identifying Druggable Targets (Finding a Workaround)
*   **Action:** A "druggable target" is a part of the cell that a medicine can physically latch onto. If the mutation completely destroyed a protein (it's missing), a drug can't latch onto something that isn't there. We have to find a workaround somewhere else in the biological assembly line.
*   **Strategies:**
    *   *The "Flood Control" Strategy (Downstream inhibition):* If the broken gene normally acts as a dam holding back a process, its absence causes a runaway biological flood. We look for a drug that targets and stops the runaway process further downstream.
    *   *The "Seesaw" Strategy (Compensation):* If the broken gene tipped the biological scales out of balance, is there *another* functioning protein we can drug to tip the scales back to normal?
    *   *The "Scaffolding" Strategy (Chaperones):* If the mutation didn't destroy the protein completely, but just bent it out of shape, we look for a drug (a "pharmacological chaperone") that acts like a brace, holding the protein in the correct shape so it works again.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

### Step 3: Computational Drug Screening
*   **Action (Finding the "Key" for our "Lock"):** In Step 2, we identified the specific protein we need to target (the "lock"). Now, we search through massive databases of known chemical compounds and FDA-approved medicines to find the specific drug (the "key") that fits into our lock and produces the desired effect (e.g., turning the protein off, or stabilizing it).
*   **Approaches:**
    *   **Knowledge Graphs:** Query databases like **DrugBank**, **ChEMBL**, or **DGIdb** (Drug Gene Interaction Database) for compounds interacting with our target pathway.
    *   **Transcriptomic Matching (LINCS L1000 / CMap):** If we know the gene expression signature of the disease, we can query the Connectivity Map to find a drug whose biological signature is the exact opposite of the disease signature.
    *   **Network Medicine:** Identify drugs that target protein-protein interaction (PPI) network neighbors of the causal gene.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

### Step 4: Rationalization and Selection
*   **Action:** Filter the list of potential drugs based on real-world constraints:
    *   Is it FDA approved?
    *   Is the toxicity profile acceptable for a pediatric patient?
    *   Can it cross the Blood-Brain Barrier (if neurological phenotypes are severe)?
*   **Selection:** Pick the 1 to 3 strongest candidates.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

### Step 5: Deliverables Generation
*   **Output:** 
    1. Write the rigorous scientific justification in `methods_description_form.xlsx` / Markdown report.
    2. Document all API calls/code used to query DrugBank/ChEMBL.
    3. Record a 3-minute pitch video explaining the biological logic.
