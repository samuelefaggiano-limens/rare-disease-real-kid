# High-Level Procedure: MVA Hackathon 2026

This document outlines the exhaustive, end-to-end bioinformatics and computational biology pipeline required to tackle both Track 1 (Variant Prediction) and Track 2 (Drug Repurposing) of the Rare Disease, Real Kid Hackathon.

## Track 1: Variant Prediction Pipeline
**Goal:** Isolate 1-2 causal genetic variants from ~5 million background variants present in the patient's Whole Genome Sequencing (WGS) VCF file.

### Step 0: Exploratory Data Analysis (EDA)
*   **Input:** All provided data files in `data/` (Clinical documents and VCF).
*   **Action:** Perform initial exploration to understand data structures, quality, and basic statistics before setting up the main pipeline.
*   **Tools:** Python, pandas, Jupyter Notebooks.
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
*   **Action:** Reduce the massive 5M row VCF file to a manageable size using intrinsic quality metrics and a fast extrinsic common-variant subtraction, prior to API-based annotation.
*   **Filter 1: Intrinsic QC Filters (Zero external downloads)**
    *   **What it does:** Relies on the sequencing machine's own confidence scores present in the VCF file.
    *   **Implementation details (based on `track1_step0_exp01-data_exploration.ipynb`):**
        *   **The PASS Filter:** The `FILTER` column contains flags for unreliable variants (e.g., `LowQual`, `MQ40`). We keep only rows where `FILTER` is missing (represented as `NaN` or `.` in pandas/cyvcf2) or `PASS`.
        *   **Genotype Quality (GQ) & Sample-Level Depth (DP):** Instead of using aggregate metrics like `INFO_DP` or `QUAL` which can be unreliable in difficult mapping regions, we prioritize sample-level confidence. We use `FORMAT_GQ` to measure the confidence that the specific 0/1 or 1/1 call for this child is correct, and `FORMAT_DP` to ensure the sequencer read this specific coordinate enough times *in this patient's sample*. We keep only rows where `FORMAT_GQ >= 20` (99% confidence) and `FORMAT_DP >= 10`.
    *   **Estimated Impact:** Reduces the 5,012,204 variants down to about 3.5 to 4 million. This is ultra-safe and clears the garbage, but the file is still too big for an API.
*   **Filter 2: The "Common Variant" Subtraction (Lightweight Extrinsic)**
    *   **What it does:** Since MVA is an ultra-rare pediatric disease, the causal mutation cannot be common in the healthy public. We use a small "blacklist" database to perform a rapid subtraction.
    *   **Implementation details:**
        *   **Dataset:** dbSNP Common or gnomAD Common VCFs (usually just 1-2 GB).
        *   **Threshold:** We filter out variants with a frequency >0.5% (or 0.1%). *Crucially, we do not use a strict 1% cutoff.* As seen in literature (e.g., Cho et al. 2025 on MOCS2, see 00-references.md file), causal recessive variants can reach >3% in specific subpopulations due to founder effects. A more conservative threshold prevents accidentally dropping localized causal mutations.
        *   **Tools:** We must perform *allele-aware* subtraction using tools like `bcftools isec -c all` or `slivar`. We must explicitly avoid coordinate-only tools like `bedtools intersect`, which would erroneously delete a novel pathogenic mutation if it simply shares a coordinate with a common harmless variant.
    *   **Estimated Impact:** Humans share most of their variation. This single, rapid step will safely strip away about 90-95% of the VCF, reducing the ~4 million high-quality variants down to roughly 100,000 to 200,000 rare/novel variants across the entire genome (coding and non-coding).
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

### Step 3: VCF Annotation
*   **Input:** Pre-filtered VCF (approx. 100,000 to 200,000 rare variants).
*   **Action:** We must annotate each surviving rare variant to understand its biological context.
*   **VCF Column Utilities (For Filtering & Verification):**
    *   *Essential:* `CHROM`, `POS`, `REF`, `ALT` (The basic GPS coordinates).
    *   *Useful for Pre-Filtering:* `FORMAT_GQ`, `FORMAT_DP`, `FILTER` (The quality bouncers).
    *   *Useful for Inheritance/Verification:* `FORMAT_GT` (Genotype, proves recessive vs dominant), `FORMAT_AD` (Allele Depths, proves variant is not an artifact).

*   **Data Added During Annotation:**
    *   **Gene Mapping:** Identifying exactly which gene (or regulatory region) the mutation is located in.
    *   **Consequence (Protein Effect):** Determining how the mutation physically alters the resulting protein (e.g., *Missense*, *Nonsense*, *Frameshift*, *Splice-site*).
    *   **Population Frequency:** Checking how common this mutation is in the general, healthy public (using massive databases like *gnomAD*).
    *   **Pathogenicity Scores:** Using computational AI and algorithms (like CADD, REVEL, SIFT, or PolyPhen) that analyze the mutation and predict how damaging it is.
*   **Tools (Annotation Engines):** We will use a "live" functional annotator like **SnpEff** or **Ensembl VEP**. 
    *   *Why not an API?* REST APIs (like `MyVariant.info`) act as lookup tables. If the child has a *de novo* completely novel mutation that has never been seen in any human before, it won't be in the API's database. The API would return empty data, causing the pipeline to blindly drop the variant later because it lacks a "Consequence". `SnpEff` and `VEP` use a reference transcriptome to *calculate* the biological consequence on-the-fly from the `REF`/`ALT` letters, guaranteeing 100% coverage even for entirely novel mutations. 
    *   APIs like `MyVariant.info` can still be used as a *secondary* step to pull down specific pre-calculated pathogenicity scores (like REVEL or CADD) if local databases are too large to host locally.
*   **Output:** 
    *   **Resulting DataFrame Structure:** The output of this annotation step will be a flattened Pandas DataFrame containing **every row (variant) from the pre-filtered VCF**, expanded with the newly fetched biological context. The expected columns include:
        *   `CHROM`, `POS`, `REF`, `ALT`, `FORMAT_GQ`, `FILTER`, `FORMAT_GT`: The original variant coordinates and alleles. We use these to track exactly where we are in the genome and to perfectly match our findings back to the submission template format required for Track 1.
        *   `GENE_SYMBOL`: The name of the affected gene (e.g., *BUB1B*, *CEP57*).
        *   `CONSEQUENCE`: The most severe predicted functional impact (e.g., `missense_variant`, `frameshift_variant`). We will use this in Step 4 to filter out "harmless" changes and keep only those predicted to severely break the protein.
        *   `GNOMAD_AF`: The overall allele frequency (AF) in the gnomAD database. We will use this in Step 4 to aggressively filter out common mutations and keep only those that are extremely rare (< 0.5%).
        *   `CADD_PHRED` & `REVEL`: Numeric pathogenicity prediction scores generated by advanced computational algorithms. We will use these in Step 4 and Step 6 to prioritize mutations that the algorithms flag as highly destructive.
        *   `CLINVAR_SIG`: Known clinical significance, if previously reported in ClinVar.
        *   `dbSNP_RSID`: The standard identifier for the variant (like a social security number, often formatted as `rs...`). Useful for quickly searching PubMed or Google Scholar later on for Track 2.
        *   `CONSERVATION_SCORE`: Measures how unchanged this specific piece of DNA has been across millions of years of evolution.
        *   `PUTATIVE_IMPACT`: Groups all mutation consequences into categories like `HIGH`, `MODERATE`, `LOW`, and `MODIFIER` (e.g., via `snpeff`). Simplifies filtering in Step 4.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

### Step 4: Variant Filtering (The Funnel)
*   **Action:** Apply a series of strict filters to reduce the millions of variants down to a few dozen candidates.
*   **Filter 1: Quality Control (Removing "Typos").** 
    *   *How/Why:* DNA sequencing machines aren't perfect; they make reading errors. This filter ensures we only keep mutations that the machine is highly confident are real, throwing out low-quality reads.
*   **Filter 2: Population Rarity (Removing Common Traits).** 
    *   *How/Why:* MVA is an ultra-rare disease. If a genetic mutation is found in 5% of the healthy public, it cannot be the cause. We use databases like gnomAD to filter out these "common" harmless variants, keeping only those that are extremely rare (e.g., < 0.1% frequency).
*   **Filter 3: Functional Impact (Removing "Harmless" Changes).** 
    *   *How/Why:* Many mutations act like harmless spelling differences (e.g., "color" vs "colour"). We filter these out, keeping only the mutations that land in crucial protein-coding regions and are predicted to severely damage or scramble the resulting protein.
*   **Filter 4: Inheritance Model (Matching the CSV Scenarios).** 
    *   *How/Why:* Humans have two copies of every gene. Depending on the exact biology of this child's MVA, we are looking for variants that fit one of the three scenarios required for submission:
        *   *Scenario 1 (Homozygous - Recessive):* The patient inherited the exact same rare mutation from both parents, breaking both copies of a gene in the same way.
        *   *Scenario 2 (Compound Heterozygous - Recessive):* The patient inherited two *different* rare mutations that landed in the exact same gene (both copies are broken, but differently).
        *   *Scenario 3 (Dominant):* A single, extraordinarily severe mutation that causes the disease all by itself, even if the other copy of the gene is healthy.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

### Step 5: Gene Prioritization & Phenotype Matching
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

### Step 6: Ranking and Submission
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
