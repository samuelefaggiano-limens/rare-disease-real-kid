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

### Step 3: VCF Annotation & Filtering (Updated)

*   **Input:** Pre-filtered VCF (`data/track1_step2_rare_variants.vcf.gz`, 591,339 variants).
*   **Action:** Run a fully integrated pipeline that annotates functional consequences, clinical significance, and molecular impact, saving the result into a single VCF.
*   **Output:** `data/track1_step3_annotation_and_filtering.vcf.gz` (905 fully enriched HIGH/MODERATE variants), ready for Step 4.

#### Stage 3.1: Functional Annotation & Impact Filtering (SnpEff)
*   **What this does:** Predicts the consequence of each variant on the protein.
*   **Tool:** SnpEff using the `GRCh38.99` database (built locally from Ensembl release 99).
*   **Implementation:** Runs with the `-canon` flag to force a single canonical transcript per gene, avoiding duplicate outputs per isoform.
*   **Filtering:** We keep only variants with `PUTATIVE_IMPACT` of `HIGH` or `MODERATE`, filtering out `LOW` and `MODIFIER` variants (~99.7% of the dataset) since rare disease-causing mutations overwhelmingly occur in coding regions or splice sites.

#### Stage 3.2: Population Frequency & Clinical Metadata (MyVariant.info)
*   **What this does:** Fetches population frequencies, in-silico pathogenicity predictions, and known clinical classifications.
*   **Tool:** MyVariant.info API (`notebooks/utils/track1_step3a_annotate.py`).
*   **Key Fields Added:**
    *   `GNOMAD_AF`, `GNOMAD_EXOME_AF` (Population frequency; dbSNP Common leaves many common coding variants behind, so this is critical for filtering AF > 1%).
    *   `CADD_PHRED` and `REVEL` (Deleteriousness and missense pathogenicity).
    *   `CLINVAR_SIG`, `CLINVAR_REVIEW_STATUS` (Clinical classifications).
    *   `OMIM_PHENOTYPE`.
*   **Important Fix:** `assembly=hg38` is explicitly requested to prevent coordinates from being wrongly matched against `hg19`.

#### Stage 3.3: AlphaGenome Atlas Scoring (AVI)
*   **What this does:** Adds a unified molecular-impact score and identifies the most relevant regulatory feature.
*   **Tool:** AlphaGenome Atlas API (`notebooks/utils/step3b_alphagenome.py`).
*   **Key Fields Added:**
    *   `AVI_PHRED`: The Phred-scaled AlphaGenome Variant Impact score (≥20 is top 1%, ≥30 is top 0.1%).
    *   `REGULATORY_FEATURE`: The top contributing modality (e.g., `CACTUS_241_WAY`, `ALPHAMISSENSE`).
*   **Limitations:** Only supports SNVs (indels are flagged as `indel_unsupported`).

**Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step3_exp01-vcf_annotating.ipynb` | Runs the full SnpEff, MyVariant.info, and AlphaGenome pipeline sequentially, producing a final annotated VCF and analysing MVA panel genes (BUB1B has 2 rare variants). | Reduced 591,339 rare variants to 905 HIGH/MODERATE variants, fully annotated and saved in `track1_step3_annotation_and_filtering.vcf.gz`. |

### Step 4: Gene Prioritization & Phenotype Matching

*   **Input:**
    *   The patient's extracted clinical HPO terms (from Step 1) in `data/hpo_terms.txt`.
    *   `track1_step3_annotation_and_filtering.vcf.gz`.
*   **Action:** Compare the genes disrupted by our candidate variants against the child's actual documented symptoms, and boost variants in genes that are already known to cause this specific disease.
*   **Implementation Details:**
    *   **Gene mapping:** each variant is linked to a gene using SnpEff's `GENE_SYMBOL` from Stage 3.2.
    *   **Matching symptoms to genes (semantic similarity):** matching the child's symptoms to a gene's known symptoms cannot rely on exact text matching, because doctors and databases describe the same underlying problem in different words (e.g. "seizures" versus the broader "abnormality of the nervous system"). We instead use the **Human Phenotype Ontology (HPO)**, a structured map of how medical terms relate to one another, and measure how *close* the child's symptoms are to each candidate gene's known symptom profile, using standard similarity algorithms (Resnik, Lin) via the `pyhpo` library. Each variant receives a `PHENOTYPE_SIMILARITY_SCORE` based on this closeness.
        We also record a `PHENOTYPE_SCORE_IS_NULL` flag for genes that have **no** documented symptom profile at all in the database. This distinction matters: a gene with no documented symptoms is simply unstudied, which is different from a gene whose documented symptoms clearly don't match this child. Both would otherwise score as "0," which would unfairly bury variants in unstudied genes — and this disease is rare enough that an unstudied gene is a real possibility.
    *   **Reference data:** we load the HPO project's `genes_to_phenotype.txt` file locally (rather than querying a live web API), and record exactly which release version of this file we used, so the analysis can be reproduced later.
    *   **The MVA gene panel (targeted boost):** MVA is primarily caused by defects in a cellular quality-control system called the Spindle Assembly Checkpoint, plus related centrosome biology. We compile a list of relevant genes and give variants in those genes a boost — but not all genes on the list are equally certain, so we use two confidence tiers rather than one:
        *   **Tier 1 — genes with confirmed causal roles in MVA specifically:** *BUB1B*, *CEP57*, *TRIP13*, *CENATAC*/*CCDC84*, *MAD2L1BP*, *CEP192*, *SMC5*, *BUB1*. These receive a strong boost.
        *   **Tier 2 — genes involved in the same general biological system, but not established causes of MVA specifically:** *MAD2L1, CDC20, CENPE, SGO1, AURKB, PLK1, MAD1L1, BUB3, TTK, ESPL1, PTTG1, NDC80, CENPF*. These receive a smaller boost. Splitting the list this way matters because some of these genes are already known to cause *other*, different diseases when disrupted — a large boost applied indiscriminately to all of them risks pushing an unrelated finding above the real answer.
        *   This gene panel, with a short note on why each gene is included, is kept in its own file (`data/mva_panel.tsv`) rather than buried in code, so it is easy to review and update.
    *   **Identifying compound-heterozygous candidates:** MVA is usually recessive, meaning it typically requires *both* copies of a gene to be broken — either the same broken copy inherited from both parents, or two different broken copies, one from each parent (a "compound heterozygous" case). We group our surviving variants by gene and flag any gene carrying two or more different candidate variants, since this is exactly the pattern the final submission format is designed to capture. Because we do not have the parents' DNA, we cannot directly prove that two such variants came from different parents — we can only note that a gene carries two plausible broken copies, and we state this limitation clearly rather than overstating our confidence.
*   **Tools:** Python (`pyhpo`, `networkx`, `pandas`), the local `genes_to_phenotype.txt` dataset, `data/mva_panel.tsv`.
*   **Output:** a master candidate table (CSV) in `track1_step4_phenotype_matching.csv`, enriched with `PHENOTYPE_SIMILARITY_SCORE`, `PHENOTYPE_SCORE_IS_NULL`, `MVA_PANEL_TIER`, and `COMPOUND_HET_CANDIDATE`, ready for final ranking in Step 5.

**Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step4_exp01-phenotype_matching.ipynb` | Extracts phenotype similarities for all genes affected by variations and identifies compound heterozygous variants. Maps variants to MVA Tier 1 and 2 panels. | Prioritized variants with phenotype similarities and saved detailed findings to `docs/02-findings.md` |

### Step 5: Ranking and Submission

*   **Input:** `data/track1_step4_candidates.tsv` (the prioritized candidate variants from Step 4).
*   **Action:** Assign an Estimated Probability of Causal Relationship (`epcr`), a number between 0 and 1, to each candidate variant, representing how confident our pipeline is that this specific variant is the true cause of the disease.
*   **How we calculate it:** a transparent, hand-weighted scoring formula that uses a unified approach for all variant types (not a trained machine-learning model — we only have one true answer in this entire project, and we don't know it, so there is no way to properly train or validate a model here; a manually weighted formula is also far easier to explain and double-check).
    *   **Pathogenicity:** Sets a baseline using `PUTATIVE_IMPACT`, which is then elevated or lowered if in-silico scores (`CADD_PHRED`, `REVEL`, `AVI_PHRED`) are present. ClinVar classifications also heavily influence this component.
    *   **Frequency:** Variants are scored on their rarity in the general population using `GNOMAD_AF`.
    *   **Inheritance:** Scores how well the genotype pattern fits the expected MVA model (e.g., homozygous or compound heterozygous).
    *   **Phenotype:** Evaluates symptom match (`PHENOTYPE_SIMILARITY_SCORE`).
    *   **MVA Panel Boost:** A 1.2x bonus is applied to variants in genes matching the MVA panel.
*   **Handling missing information correctly:** several fields will legitimately be blank for many variants:
    *   Missing in-silico scores do not negatively impact a `HIGH` impact variant's baseline score.
    *   No ClinVar entry means "unstudied," not "harmless" — it does not lower a variant's score.
    *   A missing phenotype-similarity score defaults to 0.0.
*   **Designing the submission for how it will be scored:** the competition rewards a *confident, well-separated* ranking far more than a long list of guesses. As the challenge describes it, submitting one very high-confidence correct guess alongside a handful of clearly low-confidence wrong guesses scores near-perfectly, while submitting many guesses that are all rated similarly highly scores poorly — even if the correct answer is technically included. In practical terms: submit a short list (up to 10 rows), with a clear score gap between our top pick(s) and everything else, rather than a smooth spread of similar-looking scores.
*   **Building each row of the submission file:**
    *   A variant found as two identical copies (homozygous) → goes in the `_1` columns alone, `_2` left blank.
    *   A compound-heterozygous pair identified in Step 4 → both variants go on the **same row**, one in the `_1` columns and one in the `_2` columns, sharing a single `epcr` score. A proposed pair must never be split across two separate rows — the scoring system treats each row as one complete hypothesis, and splitting a pair turns one potentially fully-correct guess into two only-half-credit guesses.
    *   A single, extremely severe variant proposed as a dominant/spontaneous cause → `_1` columns alone.
    *   The `notes` field for each row should briefly state the gene, the exact mutation, what it does to the protein, and the one-line reason we believe it's relevant.
*   **Output:** the notebook `/notebooks/track1_step5_exp01-ranking.ipynb` displays the top 10 candidates with their EPCR scores. 

**Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step5_exp01-ranking.ipynb` | Applies the EPCR scoring logic to the prioritized variants, groups compound-heterozygous pairs, and formats the top 10 candidates. | Identified BUB1B compound-heterozygous pair as top candidate (EPCR=0.9379) and exported to `output/track1/templates/track1_submission_template.csv`. |

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
| `track1_step4_exp01-phenotype_matching.ipynb` | Extracts phenotype similarities for all genes affected by variations and identifies compound heterozygous variants. Maps variants to MVA Tier 1 and 2 panels. | Prioritized variants with phenotype similarities and saved detailed findings to `docs/02-findings.md` |

### Step 2: Identifying Druggable Targets (Finding a Workaround)
*   Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step4_exp01-phenotype_matching.ipynb` | Extracts phenotype similarities for all genes affected by variations and identifies compound heterozygous variants. Maps variants to MVA Tier 1 and 2 panels. | Prioritized variants with phenotype similarities and saved detailed findings to `docs/02-findings.md` |

### Step 3: Computational Drug Screening
*   **Action (Finding the "Key" for our "Lock"):** In Step 2, we identified the specific protein we need to target (the "lock"). Now, we search through massive databases of known chemical compounds and FDA-approved medicines to find the specific drug (the "key") that fits into our lock and produces the desired effect (e.g., turning the protein off, or stabilizing it).
*   **Approaches:**
    *   **Knowledge Graphs:** Query databases like **DrugBank**, **ChEMBL**, or **DGIdb** (Drug Gene Interaction Database) for compounds interacting with our target pathway.
    *   **Transcriptomic Matching (LINCS L1000 / CMap):** If we know the gene expression signature of the disease, we can query the Connectivity Map to find a drug whose biological signature is the exact opposite of the disease signature.
    *   **Network Medicine:** Identify drugs that target protein-protein interaction (PPI) network neighbors of the causal gene.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step4_exp01-phenotype_matching.ipynb` | Extracts phenotype similarities for all genes affected by variations and identifies compound heterozygous variants. Maps variants to MVA Tier 1 and 2 panels. | Prioritized variants with phenotype similarities and saved detailed findings to `docs/02-findings.md` |

### Step 4: Rationalization and Selection
*   **Action:** Filter the list of potential drugs based on real-world constraints:
    *   Is it FDA approved?
    *   Is the toxicity profile acceptable for a pediatric patient?
    *   Can it cross the Blood-Brain Barrier (if neurological phenotypes are severe)?
*   **Selection:** Pick the 1 to 3 strongest candidates.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step4_exp01-phenotype_matching.ipynb` | Extracts phenotype similarities for all genes affected by variations and identifies compound heterozygous variants. Maps variants to MVA Tier 1 and 2 panels. | Prioritized variants with phenotype similarities and saved detailed findings to `docs/02-findings.md` |

### Step 5: Deliverables Generation
*   **Output:** 
    1. Write the rigorous scientific justification in `methods_description_form.xlsx` / Markdown report.
    2. Document all API calls/code used to query DrugBank/ChEMBL.
    3. Record a 3-minute pitch video explaining the biological logic.
