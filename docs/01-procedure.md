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
| `track1_step0_exp01-data_exploration.ipynb` | Initial exploration of VCF and Clinical data formats, missingness, and general stats. | Done: 5,012,204 raw variants, GRCh38 numeric contigs, single sample. |

### Step 1: Phenotype Extraction (HPO Mapping)
*   **Input:** `data/Challenge_Clinical_Phenotype_1.docx`
*   **Action:** Read the clinical document and extract standardized medical symptoms. We will translate free-text symptoms (e.g., "small head") into standard **Human Phenotype Ontology (HPO)** terms (e.g., HP:0000252 - Microcephaly).
*   **Tools:** Manual review, or NLP tools like Doc2HPO / ClinPhen.
*   **Output:** A list of HPO terms representing the patient's exact clinical presentation.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step1_exp01-hpo_extraction.ipynb` | Extracts the HPO terms provided in the `Challenge_Clinical_Phenotype_1.docx` file into a list of IDs. | 8 HPO terms extracted (listed in `02-findings.md`). |

### Step 2: VCF Pre-Filtering
*   **Input:** `data/WGS_EX2312012_HGWCNDSX7.vcf.gz` (Unannotated VCF, 5,012,204 variants)
*   **Action:** Reduce the massive 5M row VCF file to a manageable size using intrinsic quality metrics and a fast extrinsic common-variant subtraction, prior to API-based annotation. This single, rapid step strips away 88.2% of the VCF (the original estimate was 90-95%), leaving only rare/novel variants across the entire genome (coding and non-coding).
*   **Filter 1: Intrinsic QC Filters (Zero external downloads)**
    *   **What it does:** Relies on the sequencing machine's own confidence scores present in the VCF file.
    *   **Implementation details (based on `track1_step0_exp01-data_exploration.ipynb`):**
        *   **The PASS Filter:** The `FILTER` column contains flags for unreliable variants (e.g., `LowQual`, `MQ40`). We keep only rows where `FILTER` is missing (represented as `NaN` or `.` in pandas/cyvcf2) or `PASS`.
        *   **Genotype Quality (GQ) & Sample-Level Depth (DP):** Instead of using aggregate metrics like `INFO_DP` or `QUAL` which can be unreliable in difficult mapping regions, we prioritize sample-level confidence. We use `FORMAT_GQ` to measure the confidence that the specific 0/1 or 1/1 call for this child is correct, and `FORMAT_DP` to ensure the sequencer read this specific coordinate enough times *in this patient's sample*. We keep only rows where `FORMAT_GQ >= 20` (99% confidence) and `FORMAT_DP >= 10`.
        *   **Allelic Balance (AB):** For heterozygous calls, the fraction of reads supporting the ALT allele (`AD_ALT / (AD_REF + AD_ALT)`) must lie between 0.30 and 0.70. A true germline het should sit near 0.5; values far from it usually indicate mapping artefacts or contamination. The computed value is written as the `AB_FILTER` INFO tag. This is the largest single QC filter (248,144 variants dropped).
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
| `track1_step2_exp01-vcf_prefiltering.ipynb` | Applies Intrinsic QC & AB filters, followed by a Common Variant subtraction (MAF >= 1%) using dbSNP GRCh38, one tabix lookup per variant. | Processed 5,012,204 variants -> 4,428,273 survived intrinsic QC. Subtraction dropped 3,836,934 common variants (3,836,941 when counted per allele on multi-allelic sites, as `run_step2_prefilter.py` does; see `02-findings.md`). Final result: 591,339 ultra-rare variants saved to `track1_step2_rare_variants.vcf.gz` (405 s). |

### Step 3: VCF Annotation & Filtering (Updated)

*   **Input:** Pre-filtered VCF (`data/track1_step2_rare_variants.vcf.gz`, 591,339 variants).
*   **Action:** Run a fully integrated pipeline that annotates functional consequences, clinical significance, and molecular impact, saving the result into a single VCF.
*   **Output:** `data/track1_step3_annotation_and_filtering.vcf.gz` (905 fully enriched HIGH/MODERATE variants), ready for Step 4.

#### Stage 3.1: Functional Annotation & Impact Filtering (SnpEff)
*   **What this does:** Predicts the consequence of each variant on the protein.
*   **Tool:** SnpEff using the `GRCh38.99` database (built locally from Ensembl release 99).
*   **Implementation:** Runs with the `-canon` flag to force a single canonical transcript per gene, avoiding duplicate outputs per isoform.
*   **Filtering:** We keep only variants with `PUTATIVE_IMPACT` of `HIGH` or `MODERATE`, filtering out `LOW` and `MODIFIER` variants (99.85% of the dataset in the notebook run) since rare disease-causing mutations overwhelmingly occur in coding regions or splice sites.

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
    *   **Gene mapping:** each variant is linked to a gene using SnpEff's `GENE_SYMBOL` from Stage 3.1.
    *   **Matching symptoms to genes (semantic similarity):** matching the child's symptoms to a gene's known symptoms cannot rely on exact text matching, because doctors and databases describe the same underlying problem in different words (e.g. "seizures" versus the broader "abnormality of the nervous system"). We instead use the **Human Phenotype Ontology (HPO)**, a structured map of how medical terms relate to one another, and measure how *close* the child's symptoms are to each candidate gene's known symptom profile, using standard similarity algorithms (Resnik, Lin) via the `pyhpo` library. Each variant receives a `PHENOTYPE_SIMILARITY_SCORE` based on this closeness.
        We also record a `PHENOTYPE_SCORE_IS_NULL` flag for genes that have **no** documented symptom profile at all in the database. This distinction matters: a gene with no documented symptoms is simply unstudied, which is different from a gene whose documented symptoms clearly don't match this child. Both would otherwise score as "0," which would unfairly bury variants in unstudied genes — and this disease is rare enough that an unstudied gene is a real possibility.
    *   **Reference data:** we load the HPO project's `genes_to_phenotype.txt` file locally (rather than querying a live web API), and record exactly which release version of this file we used, so the analysis can be reproduced later.
    *   **The MVA gene panel (targeted boost):** MVA is primarily caused by defects in a cellular quality-control system called the Spindle Assembly Checkpoint, plus related centrosome biology. We compile a list of relevant genes and give variants in those genes a boost — but not all genes on the list are equally certain, so we use two confidence tiers rather than one:
        *   **Tier 1 — genes with confirmed causal roles in MVA specifically:** *BUB1B*, *CEP57*, *TRIP13*, *CENATAC*/*CCDC84*, *MAD2L1BP*, *CEP192*, *SMC5*, *BUB1*. These receive a strong boost.
        *   **Tier 2 — genes involved in the same general biological system, but not established causes of MVA specifically:** *MAD2L1, CDC20, CENPE, SGO1, AURKB, PLK1, MAD1L1, BUB3, TTK, ESPL1, PTTG1, NDC80, CENPF*. These receive a smaller boost. Splitting the list this way matters because some of these genes are already known to cause *other*, different diseases when disrupted — a large boost applied indiscriminately to all of them risks pushing an unrelated finding above the real answer.
        *   This gene panel, with a short note on why each gene is included, is kept in its own file (`data/mva_panel.tsv`) rather than buried in code, so it is easy to review and update.
    *   **Identifying compound-heterozygous candidates:** MVA is usually recessive, meaning it typically requires *both* copies of a gene to be broken — either the same broken copy inherited from both parents, or two different broken copies, one from each parent (a "compound heterozygous" case). We group our surviving variants by gene and flag any gene carrying two or more different candidate variants, since this is exactly the pattern the final submission format is designed to capture. Because we do not have the parents' DNA, we cannot directly prove that two such variants came from different parents — we can only note that a gene carries two plausible broken copies, and we state this limitation clearly rather than overstating our confidence.
*   **Tools:** Python (`pyhpo`, `networkx`, `pandas`), the local `genes_to_phenotype.txt` dataset, `data/mva_panel.tsv`.
*   **Output:** a master candidate table, `data/track1_step4_candidates.tsv` (the notebook also writes a CSV copy, `track1_step4_phenotype_matching.csv`), enriched with `PHENOTYPE_SIMILARITY_SCORE`, `PHENOTYPE_SCORE_IS_NULL`, `MVA_PANEL_TIER`, and `COMPOUND_HET_CANDIDATE`, ready for final ranking in Step 5.

**Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step4_exp01-phenotype_matching.ipynb` | Extracts phenotype similarities for all genes affected by variations and identifies compound heterozygous variants. Maps variants to MVA Tier 1 and 2 panels. | *BUB1B* ranked #2 among candidate genes (Resnik 0.632) and was the only Tier 1 panel gene with rare HIGH/MODERATE variants (two heterozygous variants, flagged compound-het). Details in `docs/02-findings.md`. |

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
*   **Output:** the notebook `notebooks/track1_step5_exp01-ranking.ipynb` saves the top 10 candidates with their EPCR scores. 

**Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| `track1_step5_exp01-ranking.ipynb` | Applies the EPCR scoring logic to the prioritized variants, groups compound-heterozygous pairs, and formats the top 10 candidates. | Identified BUB1B compound-heterozygous pair as top candidate (EPCR=0.9379) and exported to `output/track1/lostanos_hpo_matcher_v1.csv`. |

---

# Track 2: Drug Repurposing Pipeline

### Step 1: Establish the variant-specific disease mechanism

**Input:** The two BUB1B variants identified in Track 1 (one stop-gain
and one missense variant).

Characterize each allele independently and then assess their combined
effect. Determine:

-   the exact transcript/protein consequence and genomic context;
-   the location of each variant relative to known BUBR1 functional
    domains;
-   whether the stop-gain is expected to cause nonsense-mediated decay,
    truncation, or another loss-of-function mechanism;
-   whether the missense variant affects protein stability,
    localization, interactions, or function, using clinical/literature
    evidence, functional assays, domain information, and computational
    predictions as **separate evidence types**;
-   whether the two alleles plausibly result in reduced BUBR1 function
    and how this relates to the known biology of MVA1.

Build the mechanistic chain:

`BUB1B variants → BUBR1 molecular defect → spindle assembly checkpoint dysfunction → chromosome-segregation/aneuploidy phenotype`

Review evidence on BUB1B/BubR1, the spindle assembly checkpoint (SAC),
kinetochore--microtubule attachment, APC/C--CDC20 regulation, chromosome
segregation, and MVA1.

**Important:** Do not infer the functional effect of the missense
variant from the variant class alone. The strength and direction of the
evidence must be explicitly documented.

------------------------------------------------------------------------

### Step 2: Map the disease mechanism to therapeutic targets

**Input:** The mechanistic hypothesis from Step 1.

Construct an evidence-weighted network linking BUB1B dysfunction to
potentially modifiable biological processes. Examine:

-   **Direct/near-direct SAC components:** BUB3, MAD2L1, BUB1, CDC20,
    MAD1L1, TTK/MPS1;
-   **APC/C and cell-cycle regulation:** APC/C components, CDC20, and
    relevant downstream regulators;
-   **Kinetochore/chromosome segregation:** CENPE, NDC80, SGO1, AURKB,
    PLK1, ESPL1, and other experimentally supported components;
-   **aneuploidy/chromosome-instability vulnerabilities:** dependencies
    or synthetic-lethal/selective-vulnerability mechanisms reported in
    relevant models.

Distinguish clearly between:

1.  direct BUB1B interactors;
2.  members of the same biological pathway;
3.  genes associated with MVA1;
4.  indirect relationships through broader cell-cycle, aneuploidy, or
    chromosome-instability biology.

For each potential target, define the **desired direction of
modulation** and the rationale. Consider three therapeutic concepts:

-   **functional compensation** --- preserve or enhance a process that
    could compensate for reduced BUBR1 activity;
-   **downstream phenotype correction** --- reduce chromosome
    mis-segregation or restore cell-cycle/checkpoint control;
-   **phenotype/vulnerability modulation** --- exploit a
    disease-associated dependency.

Do not assume that inhibiting a pathway component is beneficial simply
because BUB1B regulates it. In particular, APC/C--CDC20 inhibition
should be treated as a testable hypothesis, not an established
therapeutic strategy.

**Output:** A short list of prioritized therapeutic nodes, each linked
to a specific mechanistic rationale and evidence level.

------------------------------------------------------------------------

### Step 3: Map targets to approved drugs and integrate computational evidence

**Input:** Prioritized therapeutic targets from Step 2.

Identify existing approved medications capable of modulating the
selected targets using resources such as **ChEMBL, DrugCentral, DGIdb,
Open Targets, PubChem, DrugBank (where accessible), ClinicalTrials.gov,
and the primary literature**.

For every drug--target hypothesis, record:

-   drug and target;
-   mechanism and direction of modulation;
-   approved indication and regulatory status;
-   evidence for the drug--target interaction;
-   evidence connecting the target to the BUB1B/MVA1 mechanism;
-   evidence that the proposed direction of modulation could produce a
    beneficial cellular effect;
-   relevant human safety/pharmacology information.

Apply the key filtering rule:

`BUB1B dysfunction → biological consequence → therapeutic target → drug mechanism → expected beneficial effect`

A drug should **not** be prioritized merely because it interacts with a
protein somewhere in the BUB1B pathway.

Use transcriptomic and network-based analyses as **supporting
evidence**, not as the primary basis for candidate selection. Where
suitable, use CMap/LINCS to test whether drug-induced expression
signatures are compatible with reversing disease-relevant molecular
states, alongside pathway enrichment, gene-set analysis, PPI networks,
and drug--target network analysis.

Computational similarity alone is insufficient to establish therapeutic
efficacy; CMap/LINCS results should be treated as hypothesis-supporting
evidence because transcriptomic signatures can be context-, dose-, and
cell-type-dependent.

> **Status:** the CMap/LINCS check was **not performed**. Within the hackathon timeline we prioritised the drug-target mapping and the MVA-wide scan, and the candidate compounds (NAD+ precursors in particular) have little or no L1000 profiling. No hypothesis in the final ranking relies on transcriptomic evidence. The plan above is kept as the intended next step (see the report, §9.1).

**Output:** A candidate-drug table with a transparent evidence trail
from BUB1B dysfunction to the proposed drug effect.

------------------------------------------------------------------------

### Step 4: Integrate evidence, assess feasibility, and define validation

Integrate the candidates using a transparent, non-clinical
prioritization framework covering:

-   **mechanistic relevance** --- proximity of the target to the BUB1B
    defect;
-   **directionality** --- evidence that the drug acts in the desired
    direction;
-   **experimental evidence** --- BUB1B-deficient,
    SAC/chromosome-segregation, aneuploidy, or CIN models;
-   **human evidence** --- clinical use, pharmacology, and safety;
-   **drug feasibility** --- established mechanism, pharmacokinetics,
    formulation, route, interactions, contraindications, and pediatric
    evidence where available.

Evidence from direct BUB1B-deficient models should be distinguished from
evidence in related SAC/aneuploidy/cancer models and from purely
computational predictions.

Select approximately **1--3 final hypotheses**. For each, report:

1.  BUB1B defect and affected process;
2.  therapeutic target and drug mechanism;
3.  rationale for the proposed direction of modulation;
4.  supporting evidence and its limitations;
5.  human/pediatric safety and translational considerations;
6.  a proposed experiment.

Where feasible, compare **BUB1B-mutant and control cellular systems**
under vehicle and candidate-drug treatment. Potential readouts include
mitotic progression, chromosome-segregation errors, micronuclei,
aneuploidy, cell viability, SAC activation, and relevant molecular
biomarkers. Explicitly distinguish **phenotypic rescue** from
**selective toxicity**.

------------------------------------------------------------------------

### Track 2 Experiments

Track 2 is implemented as scripts in `notebooks/utils/` rather than as
`.ipynb` experiments. The pipeline is a deterministic sequence of public
database queries with no exploratory or patient-data component, so the
notebook format adds no reproducibility benefit here; the scripts are
runnable end-to-end and their outputs are committed. Steps 1 and 2 were
produced as a literature analysis with no code.

| Step | Implementation | Description | Results |
| :--- | :--- | :--- | :--- |
| 1 · Mechanism | *(literature analysis)* | Characterises each allele and the combined effect: L751* stop-gain expected to undergo NMD (null allele), N1016K missense in the C-terminal pseudokinase domain required for PP2A-B56 recruitment. | Mechanistic chain established and recorded in `docs/02-findings.md`. |
| 2 · Targets | *(literature analysis)* | Evidence-weighted network from BUB1B dysfunction to modifiable processes, across SAC components, APC/C regulation, kinetochore biology and aneuploidy vulnerabilities. | Five prioritised therapeutic nodes. **Four were subsequently revised or discarded in Step 3** (see below). The two hypotheses finally prioritised in Step 4 (NAD+/SIRT2, senescence) were not among these five: they come from the BubR1-dosage and senescence literature in `00-references.md` §4-5, added during Steps 3-4. |
| 3 · Drugs | `notebooks/utils/track2_step3_drug_mapping.py` + `track2_target_panel.tsv` | Maps a 35-protein panel to approved drugs via ChEMBL (primary), DGIdb (cross-check) and STRING (network distance). Records direction concordance against the desired modulation. | **All 18 proteins at network distance 1 from BUB1B have zero approved concordant drugs**; all 26 approved concordant drugs found act on proteins at distance 2-3 or with no STRING path to BUB1B. Outputs in `output/track2/`. |
| 3 · CMap/LINCS | — | Transcriptomic signature check proposed above as supporting evidence. | **Not performed.** No candidate in the current ranking relies on it; the Step 4 scores use only literature, ChEMBL/DGIdb and STRING evidence. Not run within the hackathon timeline; kept as a proposed follow-up (report §9.1). |
| 3 · Figure | `docs/figures/build_track2_mechanisms_figure.py` | Documents all 12 routes with the reason each is retained or discarded, in five categories. EN/ES, same build-script convention as the Track 1 figure. | `docs/figures/track2_mechanisms_explained.html` + PNG. |
| 4 · Prioritisation | `notebooks/utils/track2_step4_prioritize.py` + `track2_step4_candidates.tsv` | Transparent weighted scoring on the five axes above, with mechanistic relevance and directionality dominant so drug availability alone cannot carry a candidate. | Two prioritised hypotheses: NAD+/SIRT2 BubR1 stabilisation (0.817) and senescent-cell clearance (0.667). Output in `output/track2/`. |
| 5 · Deliverables | `output/track2/lostanos_track2_report.md`, `lostanos_track2_methods_description.xlsx`, `lostanos_track2_pitch_script.md` | Report, methods form, AI-usage disclosure, pitch video. | Report and methods form complete (AI-usage disclosure filled in; ~20 person-hours); pitch video recorded and edited (3:22), upload to YouTube/Vimeo pending. |
| 6 · MVA-wide scan | `notebooks/utils/track2_step6_mva_druggability_scan.py` | Scalability check: the same drug lookup run automatically (no curated panel, direction-agnostic) on up to 15 high-confidence STRING interactors of all eight genes with a reported MVA phenotype. | Only 2 of 108 first-shell interactors have any approved drug, neither relevant to chromosome segregation; no MVA gene is itself druggable. ~4 min. Outputs `output/track2/track2_step6_mva_scan_*.csv`. |

**Reproducing Track 2:**

```bash
python notebooks/utils/track2_step3_drug_mapping.py   # ~5 min, queries ChEMBL/DGIdb/STRING
python notebooks/utils/track2_step4_prioritize.py     # instant, offline
python notebooks/utils/track2_step6_mva_druggability_scan.py  # ~4 min, all 8 MVA genes
python docs/figures/build_track2_mechanisms_figure.py # regenerates the figure
```

Step 3 depends on external APIs. UniProt accessions are cached in
`notebooks/utils/track2_uniprot_accessions.json` because the service
throttles bursts and returned HTTP 500 mid-run; ChEMBL and DGIdb are
queried live with exponential backoff. No patient data is involved at any
point in Track 2.
