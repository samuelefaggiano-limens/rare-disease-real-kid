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
| `track1_step2_exp01-vcf_prefiltering.ipynb` | Applies Intrinsic QC & AB filters, followed by a Common Variant subtraction (MAF >= 1%) using dbSNP GRCh38. | Processed 5M variants -> 2.64M survived intrinsic QC. Subtraction dropped 4.4M common variants. Final result: 591,339 ultra-rare variants saved to `track1_step2_rare_variants.vcf.gz`. |

### Step 3: VCF Annotation & Scoring (Parallel Tracks)
*   **Input:** Pre-filtered VCF (`track1_step2_rare_variants.vcf.gz` containing ~591,339 variants).
*   **Action:** We will run two independent annotation pipelines in parallel to evaluate the strengths of traditional functional annotation versus deep-learning regulatory scoring, and then compare their top candidates.

#### Step 3A: Traditional Functional Annotation (SnpEff & MyVariant.info)
*   **Input:** Pre-filtered VCF.
*   **Action:** We must annotate each surviving rare variant to understand its biological context.
*   **Data Added During Annotation:**
    To successfully evaluate each variant, we integrate existing VCF metrics with new biological context annotations:
    *   **Original VCF Columns Maintained:** 
        *   `CHROM`, `POS`, `REF`, `ALT`: The fundamental coordinates and alleles.
        *   `FORMAT_GQ`, `FORMAT_DP`, `FILTER`: Our quality bouncers for filtering.
        *   `FORMAT_GT`, `FORMAT_AD_REF`, `FORMAT_AD_ALT`: Genotype and allele depths, crucial for proving recessive vs. dominant inheritance and confirming variant legitimacy.
    *   **New Annotations Added:**
        *   `GENE_SYMBOL`: The affected gene (e.g., *BUB1B*). Found by mapping the variant coordinates to a reference genome.
        *   `CONSEQUENCE`: The physical alteration to the protein (e.g., *Missense*, *Frameshift*). Calculated dynamically based on how the `REF` changes to `ALT` inside the transcript.
        *   `PUTATIVE_IMPACT`: Categorizes consequences into `HIGH`, `MODERATE`, `LOW`, or `MODIFIER`. Sourced from functional annotator calculations.
        *   `GNOMAD_AF`: The overall allele frequency in the human population. 
        *   `CADD_PHRED` & `REVEL`: Numeric pathogenicity scores evaluating how destructive the mutation is. 
        *   `CLINVAR_SIG` & `CLINVAR_REVIEW_STATUS`: Known clinical significance and the evidence star-rating. 
        *   `OMIM_PHENOTYPE`: Inherited syndromes linked to the gene. 
        *   `LOFTOOL_SCORE`: Gene's intolerance to loss-of-function. 
*   **Tools:**
    *   **Functional Annotation:** We will use a tool like **SnpEff** locally. These "live" annotators use a reference transcriptome to dynamically calculate the `GENE_SYMBOL`, `CONSEQUENCE`, and `PUTATIVE_IMPACT` on the fly from the `REF`/`ALT` alleles, ensuring 100% coverage even for entirely novel mutations. SnpEff databases for GRCh38 are roughly ~5GB.
    *   **Filter:** Keep only variants where `PUTATIVE_IMPACT` is `HIGH` or `MODERATE` (drastically reducing the ~591k variants to roughly ~10k-20k).
    *   **Local Transcriptome Annotation (SnpEff):** Run SnpEff locally using the GRCh38 database with canonical transcript enforcement `(-canon)` to avoid transcript duplication and retain the single most severe biological consequence per variant. Dynamically derive: `GENE_SYMBOL`, `CONSEQUENCE`, and `PUTATIVE_IMPACT`.
    *    **High-Throughput Metadata Cross-Referencing:** Query the MyVariant.info REST API strictly via its batch endpoint `(POST /v1/variant`, batch size = 1,000) using a Python pipeline with retry and backoff logic to prevent rate-limiting. To append the remaining metadata (`GNOMAD_AF`, `CADD_PHRED`, `REVEL`, `CLINVAR_SIG`, `OMIM_PHENOTYPE`, etc.), we will use query an external REST API (like **MyVariant.info**) to fetch and append these fields programmatically via Python. Using APIs is preferred. SnpEff will provide `GENE_SYMBOL`, `CONSEQUENCE` and `PUTATIVE_IMPACT`. `MyVariant.info` will provide the remaining fields. 
*   **Output:**
    The final output will be a VCF file containing exactly one row per variant (prioritizing the canonical transcript's most severe consequence to avoid duplicates). The final columns will be:
    *   **From Input VCF:** `CHROM`, `POS`, `REF`, `ALT`, `FILTER`, `FORMAT_DP`, `FORMAT_GQ`, `FORMAT_GT`, `FORMAT_AD_REF`, `FORMAT_AD_ALT`
    *   **Added at this Step:** `GENE_SYMBOL`, `CONSEQUENCE`, `PUTATIVE_IMPACT`, `GNOMAD_AF`, `CADD_PHRED`, `REVEL`, `CLINVAR_SIG`, `CLINVAR_REVIEW_STATUS`, `OMIM_PHENOTYPE`, `LOFTOOL_SCORE`
    *   **Format:**  Multi-sample annotated VCF (metadata written into standard VCF `INFO` tags: `CSQ`, `CADD`, `CLINVAR_SIG`, etc.) alongside a flattened, unified tabular file `(track1_step3a_coding_annotated.tsv)` for downstream ranking.

#### Step 3B: Deep Learning Regulatory Scoring (AlphaGenome Atlas)
*   **Input:** Pre-filtered VCF (`track1_step2_rare_variants.vcf.gz`).
*   **Action:** Evaluate the rare variants across the *entire* genome (both coding and non-coding regions) for structural and regulatory impact using Google DeepMind's AlphaGenome Atlas API.
*   **Optimization & Pre-filtering:** To bypass API latency across ~591k raw variants, apply a regulatory region pre-filter using Bedtools against candidate cis-regulatory elements (ENCODE cCREs, UTRs, promoters, and +/-250 bp canonical/cryptic splice junctions), reducing input variants to biologically actionable regions (~40k–60k variants).
    *   *Why this complements 3A:* Traditional tools like SnpEff are excellent at identifying variants that change an amino acid (protein-coding changes). However, AlphaGenome excels at predicting regulatory effects (e.g., altered transcription factor binding, chromatin accessibility, or RNA splicing). Note that AlphaGenome is not restricted to the non-coding genome it also evaluates variants *inside* coding regions to see if they cause hidden regulatory disruptions, like breaking an exonic splice site or altering expression.
*   **Data Added During Annotation:**
    *   `AVI_SCORE` (AlphaGenome Variant Impact): A unified numeric score predicting the pathogenicity and regulatory disruption of the variant.
    *   `TARGET_GENE`: The specific gene(s) whose expression or splicing is predicted to be altered by the variant (crucial for Step 4 Phenotype Matching).
    *   `REGULATORY_FEATURE`: The specific mechanism or non-coding element disrupted (e.g., Splicing, Enhancer, Promoter) to inform Track 2 mechanism characterization.
*   **Execution:** Query the AlphaGenome Atlas client `(alphagenome.atlas.atlas.AtlasClient)` in chunks of 5,000 variants:
    *   client.query_variants(variants_chunk, requested_scorers=['avi'])
    *   Extract scores from the resulting `AnnData` objects.
*   **Tools:**
    *   **AlphaGenome Atlas API** (`alphagenome.atlas.atlas.AtlasClient`).
    *   We will use `AtlasClient.query_variants(variants, requested_scorers=['avi'])` to retrieve the scores from the returned `AnnData` matrix and merge them back into our pipeline.
    *   **Filter:** Retain only variants that meet a high `AVI_SCORE` threshold (e.g., top 1%) to isolate the most destructive candidates with significant regulatory impact.
*   **Output:** 
    An annotated dataset (VCF file) containing our rare variants appended with their respective AlphaGenome scores, target genes, and regulatory metadata.
    Tabular dataset `(track1_step3b_regulatory_annotated.tsv)` filtered to top-tier candidates (e.g., `AVI_SCORE` ≥99thpercentile of scored variants).

#### Step 3C: Comparison & Candidate Selection
*   **Action:** We will merge the top-ranked candidates from both 3A (ranked by CADD/REVEL + Phenotype) and 3B (ranked by AVI score). This will allow us to see if AlphaGenome identifies a highly impactful non-coding regulatory variant that the traditional SnpEff pipeline completely ignored. 
* **Classification Logic:**  
    * **Concordant High-Confidence Candidates (Coding + Regulatory):** Variants flagged with `PUTATIVE_IMPACT` (HIGH/MODERATE) that also display high `AVI_SCORE`, highlighting coding alterations that simultaneously disrupt local secondary structures, exonic splice enhancers (ESEs), or transcription factor binding.
    *  **Non-Coding Regulatory Rescues:** Variants discarded by Step 3A (annotated as `LOW` or `MODIFIER`) but scoring in the top tier of `AVI_SCORE` with a clear, predicted `TARGET_GENE`.
    *  **Clinical & Conservation Concordance:** Pathogenic/Likely Pathogenic ClinVar records or variants exceeding thresholds (`CADD_PHRED`≥ 20,`REVEL`≥0.5) advance directly to priority status.

  Final Deliverable:A unified candidate table `(track1_step3_prioritized_candidates.parquet / .tsv)` containing sample genotypes, variant coordinates, combined pathogenicity/regulatory scores, and mapped gene symbols ready for clinical phenotype intersection.


* **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

### Step 4: Gene Prioritization & Phenotype Matching
*   **Action:** We now have a small list of rare, protein-altering mutations. We cross-reference the genes containing these mutations with the HPO terms generated in Step 1.
*   **Matching:** We compare the patient's extracted HPO terms against established "Gene-to-Phenotype" databases (like OMIM or the HPO database itself). These databases contain lists of symptoms (HPO terms) that are scientifically proven to occur when specific genes are mutated. By calculating the overlap (or semantic similarity) between the child's symptoms and a gene's known disease profile, we can rank the candidate genes.
    *   *Example of Cross-referencing:* 
        1.  **Patient Data:** The clinical text (Step 1) gives us `HP:0000252` (Microcephaly) and `HP:0001250` (Seizures).
        2.  **Variant Data:** Step 4 filtering leaves us with rare mutations in three candidate genes: *Gene A*, *Gene B*, and *BUB1B*.
        3.  **Database Lookup:** We query the HPO database for these three genes.
            *   *Gene A* is only known to cause `HP:0010442` (Polydactyly - extra fingers). *Match: 0% overlap*.
            *   *Gene B* is known to cause `HP:0000365` (Deafness). *Match: 0% overlap*.
            *   *BUB1B* is known to cause `HP:0000252` (Microcephaly), `HP:0001250` (Seizures), and `HP:0000518` (Cataracts). *Match: High overlap!*
        4.  **Conclusion:** Because the known symptoms of a broken *BUB1B* gene perfectly align with our patient's symptoms, the mutation in *BUB1B* is prioritized as the most likely cause, effectively discarding the "noisy" mutations in *Gene A* and *Gene B*.
*   **Targeted Search:** For this specific hackathon, we also explicitly check if the surviving variants land in known MVA genes (e.g., *BUB1B*, *CEP57*, *TRIP13*) or genes related to the **Spindle Assembly Checkpoint (SAC)** or centrosomes.
*   **Tools:** Exomiser, Phenolyzer, or manual literature querying via PubMed.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

### Step 5: Ranking and Submission
*   **Action:** Assign an Estimated Probability of Causal Relationship (`epcr`) between 0 and 1 to each candidate variant. 
    *   *What is EPCR?* This is not a standard biological metric; it is a custom confidence score we must engineer for this hackathon. It represents our pipeline's belief that a specific variant is the true cause of the disease.
    *   *How we calculate it:* We will build a scoring function (either a weighted heuristic formula or an ML model) that combines the variant's features into a single score. For example, a variant gets a higher `epcr` if it has:
        1.  Near-zero population frequency (gnomAD).
        2.  A high computational pathogenicity score (e.g., CADD, REVEL).
        3.  A high semantic overlap with the patient's HPO terms (Phenotype matching).
        4.  A zygosity that perfectly fits an autosomal recessive or dominant inheritance model.
*   **Output:** Format the top candidates, ranked from highest `epcr` to lowest, into `output/track1/track1_submission_template.csv`.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

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
