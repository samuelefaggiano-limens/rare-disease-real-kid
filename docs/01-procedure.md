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

### Step 2: VCF Pre-Filtering & Annotation
*   **Input:** `data/WGS_EX2312012_HGWCNDSX7.vcf.gz` (Unannotated VCF, ~5M variants)
*   **Action:** The current VCF only has coordinates. We must annotate each variant to understand its biological context. 
*   **VCF Column Utilities (For Filtering & Verification):**
    *   *Essential:* `CHROM`, `POS`, `REF`, `ALT` (The basic GPS coordinates).
    *   *Useful for Pre-Filtering:* `QUAL`, `INFO_DP`, `FILTER` (The quality bouncers).
    *   *Useful for Inheritance/Verification:* `FORMAT_GT` (Genotype, proves recessive vs dominant), `FORMAT_AD` (Allele Depths, proves variant is not an artifact).
    *   *Useless/Redundant:* `INFO_MQ`, `INFO_FS`, `FORMAT_DP`, `FORMAT_GQ`, `ID` (handled implicitly by our QC bouncers or fetched later).

*   **Data Added During Annotation:**
    *   **Gene Mapping:** Identifying exactly which gene (or regulatory region) the mutation is located in. Most of our DNA doesn't code for proteins, so we need to know if the mutation hit a critical "instruction manual" (a gene) or just empty space.
    *   **Consequence (Protein Effect):** Determining how the mutation physically alters the resulting protein. For example:
        *   *Missense:* Changes a single "word" (amino acid) in the instructions, which might mildly or severely affect the protein's shape.
        *   *Nonsense:* Introduces an abrupt "stop sign," resulting in a truncated, usually useless protein.
        *   *Frameshift:* Inserts or deletes letters in a way that scrambles the entire rest of the sentence (highly destructive).
        *   *Splice-site:* Disrupts the delicate process of how the gene's instructions are cut and pasted together before being read.
    *   **Population Frequency:** Checking how common this mutation is in the general, healthy public (using massive databases like *gnomAD*). Since MVA is an ultra-rare disease, the mutation causing it should be extremely rare or completely novel. If 5% of healthy people have a mutation, it is harmless background noise, not the cause of MVA.
    *   **Pathogenicity Scores:** Using computational AI and algorithms (like CADD, REVEL, SIFT, or PolyPhen) that analyze the mutation and predict how damaging it is. These tools look at evolutionary conservation and protein structure to give a score indicating whether the mutation is likely "benign" or highly "pathogenic" (damaging).
*   **Tools (Annotation Engines):** **MyVariant.info API (Python)**, A REST API that aggregates data from ClinVar, gnomAD, dbNSFP, and Ensembl into a single JSON response. Because this hackathon focuses on a single patient VCF and emphasizes reproducible Python pipelines without massive local database dependencies, we will programmatically fetch the annotations using the **MyVariant.info API (via Python)**. This allows us to query coordinates directly and retrieve all necessary biological context in a lightweight, easily parseable manner.
*   **Output:** 
    *   **Resulting DataFrame Structure:** The output of this annotation step will be a flattened Pandas DataFrame containing **every row (variant) from the original VCF**, but expanded with the newly fetched biological context as additional columns. The expected columns include:
        *   `CHROM`, `POS`, `REF`, `ALT`, `QUAL`, `FILTER`, `GT`: The original variant coordinates and alleles. These act as the fundamental "GPS coordinates" for the mutation. **Why we need it:** We use these to track exactly where we are in the genome and to perfectly match our findings back to the submission template format required for Track 1.
        *   `GENE_SYMBOL`: The name of the affected gene (e.g., *BUB1B*, *CEP57*). **Why we need it:** Most of human DNA does not contain instructions for making proteins. We need to know which specific "instruction manual" (gene) the mutation landed in so we can cross-reference it with the child's symptoms in Step 4.
        *   `CONSEQUENCE`: The most severe predicted functional impact (e.g., `missense_variant`, `frameshift_variant`). **Why we need it:** Not all mutations are equally destructive. This column tells us if the mutation just changes a single "letter" (missense) or completely scrambles the rest of the sentence (frameshift). We will use this in Step 3 to filter out "harmless" changes and keep only those predicted to severely break the protein.
        *   `GNOMAD_AF`: The overall allele frequency (AF) in the gnomAD database, which represents how common the mutation is in the healthy general public. **Why we need it:** MVA is an ultra-rare disease. If a mutation is found in 5% of healthy people, it is just harmless background noise. We will use this in Step 3 to aggressively filter out common mutations and keep only those that are extremely rare (< 0.1%).
        *   `CADD_PHRED` & `REVEL`: Numeric pathogenicity prediction scores generated by advanced computational algorithms. **Why we need it:** These AI-driven scores estimate how "damaging" a mutation is to the resulting protein's structure or function. We will use these in Step 3 and Step 5 to prioritize mutations that the algorithms flag as highly destructive.
        *   `CLINVAR_SIG`: Known clinical significance, if previously reported in ClinVar (a public archive of human genetic variants and their relationship to disease). **Why we need it:** If this exact mutation has already been seen in another patient and classified as "Pathogenic" by geneticists, it strongly supports our case. Conversely, if it is classified as "Benign," we can likely ignore it.
        *   `dbSNP_RSID` (`dbsnp.rsid`): The standard identifier for the variant (like a social security number, often formatted as `rs...`). **Why we need it:** Having a standard ID makes it much easier to quickly search PubMed, Google Scholar, or clinical databases for existing literature, which is highly useful later on for Track 2 (Drug Repurposing).
        *   `CLINVAR_PHENOTYPE` (`clinvar.rcv.conditions.name`): The name of the specific disease associated with this variant in the ClinVar database. **Why we need it:** If the variant's condition is literally listed as "Mosaic variegated aneuploidy", it's an immediate smoking gun for our pipeline. 
        *   `CONSERVATION_SCORE` (e.g., `dbnsfp.gerp++_rs`): Measures how unchanged this specific piece of DNA has been across millions of years of evolution (across species like mice, dogs, fish). **Why we need it:** If a genetic coordinate hasn't changed across 100 vertebrate species, it implies it is biologically critical. A mutation there is highly likely to be pathogenic (bad).
        *   `SPLICE_SCORE` (e.g., `dbnsfp.ada_score`): Predicts if the mutation will disrupt the cell's ability to "splice" (cut and paste) the gene's instructions together. **Why we need it:** These help catch sneaky mutations that don't directly change an amino acid word but mess up the overall sentence structure of the gene, which is a frequent cause of rare genetic diseases.
        *   `CONSENSUS_PATHOGENICITY` (`dbnsfp.sift.pred`, `dbnsfp.polyphen.pred`, `dbnsfp.mutationtaster.pred`): Pathogenicity predictions from other independent models (SIFT, PolyPhen, MutationTaster). **Why we need it:** While CADD and REVEL are excellent, relying on a single algorithm can lead to false positives. By checking multiple tools, if they *all* agree a mutation is "Damaging," we can trust the prediction much more.
        *   `PUTATIVE_IMPACT` (`snpeff.ann.putative_impact`): Groups all mutation consequences into four categories: `HIGH`, `MODERATE`, `LOW`, and `MODIFIER`. **Why we need it:** Simplifies filtering in Step 3. Instead of writing complex rules for 30 different consequence strings, we can drop all variants where impact is 'LOW' or 'MODIFIER'.
        *   `GNOMAD_HOMOZYGOUS_COUNT` (`gnomad_genome.hom.hom` / `gnomad_exome.hom.hom`): The exact number of healthy individuals in gnomAD who possess two copies (homozygous) of this mutation. **Why we need it:** MVA is an ultra-rare, severe pediatric disease. If healthy adults in gnomAD are homozygous for a mutation, it cannot be the cause of our patient's disease under Scenario 1 (Homozygous Recessive).
        *   `PROTEIN_DOMAIN` (`dbnsfp.interpro_domain`): Indicates if the mutation lands inside a known, critical working part of the protein (like a "kinase domain"). **Why we need it:** Crucial for Track 2 (Drug Repurposing). Knowing exactly which part of the protein broke helps us identify what kind of biological activity is missing and what kind of drug could compensate.
*   **Experiments:**

| Notebook | Description | Results |
| :--- | :--- | :--- |
| | | |

### Step 3: Variant Filtering (The Funnel)
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

### Step 4: Gene Prioritization & Phenotype Matching
*   **Action:** We now have a small list of rare, protein-altering mutations. We cross-reference the genes containing these mutations with the HPO terms generated in Step 1.
*   **Matching:** We compare the patient's extracted HPO terms against established "Gene-to-Phenotype" databases (like OMIM or the HPO database itself). These databases contain lists of symptoms (HPO terms) that are scientifically proven to occur when specific genes are mutated. By calculating the overlap (or semantic similarity) between the child's symptoms and a gene's known disease profile, we can rank the candidate genes.
    *   *Example of Cross-referencing:* 
        1.  **Patient Data:** The clinical text (Step 1) gives us `HP:0000252` (Microcephaly) and `HP:0001250` (Seizures).
        2.  **Variant Data:** Step 3 filtering leaves us with rare mutations in three candidate genes: *Gene A*, *Gene B*, and *BUB1B*.
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
