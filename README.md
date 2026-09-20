# Rare Disease, Real Kid: MVA Hackathon 2026

Repository for participating in the [SageBio/rare-disease-real-kid-mva-hackathon-2026](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026) Hugging Face Space.

## Overview
The genome and clinical story here belong to a real child living with Mosaic Variegated Aneuploidy (MVA), an ultra-rare genetic condition affecting fewer than 50 people worldwide.

---

## Track 1 - Variant Prediction
**Goal:** Predict the specific genetic variant(s) driving the child's Mosaic Variegated Aneuploidy (MVA) using their clinical symptoms (phenotype) and genomic data (genotype). 

This is a foundational track designed to validate that our pipeline can recover a clinically confirmed variant from real patient data among millions of possibilities. 
**The Ground Truth:** The organizers already know the exact genetic variant causing this child's disease because a team of medical geneticists previously analyzed the child's DNA and confirmed the diagnosis. Our challenge is to build an algorithm to find that exact same "needle in the haystack" using only the provided data.

### The Input Data
The necessary data has been downloaded and is available in the local `data/` directory.

#### 1. Clinical Data
*   **File:** `data/Challenge_Clinical_Phenotype_1.docx` (17 KB)
*   **What it is:** A detailed clinical description of the child's observable symptoms and medical history (likely containing standardized HPO terms). This is crucial for guiding our search toward specific biological pathways.

#### 2. Processed Genomic Data (VCF)
*   **Files:** `data/WGS_EX2312012_HGWCNDSX7.vcf.gz` (301 MB) and its index `.tbi` (2.3 MB)
*   **What it is:** A highly refined list of genetic differences. The variants in this file were processed and filtered using **GATK VariantFiltration (v4.2.4.0)** with standard hard filters (e.g., `QD < 2.0`, `MQ < 40.0`). The coordinates match the **GRCh38** reference genome assembly (using numeric contig notation like `1`, `2` rather than `chr1`).
*   **How we use it:** **This is our primary dataset.** We will filter these millions of variants down to find the ultra-rare, pathogenic mutation causing the disease.

#### 3. Raw Genomic Data (FASTQ)
*   **Files:** `data/WGS_EX2312012_HGWCNDSX7_S16_L00[1-4]_R[1-2]_001.fastq.gz` (8 paired-end read files spanning 4 lanes, totaling ~82 GB)
*   **What it is:** The raw, unaligned output directly from an Illumina DNA sequencing machine. 
*   **How we use it:** We will likely ignore this for Track 1, as the VCF has already done the heavy lifting of identifying the variants.

### The Output (Deliverables)
We are permitted up to **6 submissions**. Each submission requires three components:
#### 1. Predictions File Format (`track1_submission_template.csv`)
The CSV must contain one row per proposed causal variant (or compound-heterozygous pair). 
*   **Ranking:** Rows must be ranked by our Estimated Probability of Causal Relationship (`epcr`), where 1 is the most likely and 0 is the least.
*   **File Naming:** `<username>_<model-name>.csv`.
*   **Key Fields:** 
    *   `proband_id`: Provided in the dataset.
    *   `chrom_1`, `pos_1`, `ref_1`, `alt_1`: Coordinates (GRCh38) and alleles for the primary variant.
    *   `chrom_2`, `pos_2`, `ref_2`, `alt_2`: Used *only* for compound-heterozygous pairs.
    *   `epcr`: Our confidence score (0 to 1].
    *   `finding_type`: "primary" (MVA candidate) or "secondary" (incidental findings).
    *   `notes`: Brief rationale.

##### Mapping Biology to the CSV Rows
MVA is typically an **autosomal recessive** condition, meaning it usually requires *both* copies of a gene to be broken. Here is how the genetics translate to our submission:
*   **Scenario 1: Homozygous (Same mutation from both parents):** Both parents passed down the exact same mutation. We put this single variant into `chrom_1`, `pos_1`, `ref_1`, `alt_1` and leave the `_2` columns completely blank.
*   **Scenario 2: Compound Heterozygous (Two different mutations in the same gene):** Mom passed down a mutation at Position A, and Dad passed down a different mutation at Position B inside the *same gene*. We put Mom's variant in the `_1` columns and Dad's variant in the `_2` columns on the **exact same row**. This tells the judges: *"These two variants are working together to break the gene."*
*   **Scenario 3: Dominant / De Novo (One extremely bad mutation):** A single heterozygous mutation is so severe it causes the disease all by itself. We put this single variant in the `_1` columns and leave the `_2` columns blank.

The ultimate goal is for **Row #1** of our CSV to contain either the single causal variant (Scenario 1 or 3) OR the pair of compound-heterozygous variants (Scenario 2).

#### 2. Methods Report (PDF or Markdown)
 A write-up explaining our approach, any AI tools used, and datasets referenced (format template: `output/track1/methods_description_form.xlsx`).
#### 3. Code Repository
 This GitHub repository containing the reproducible pipeline.

### Evaluation & Scoring
The condition might be caused by a single variant (dominant) or two variants (compound heterozygous recessive). Our CSV is automatically scored live on the leaderboard against the clinically confirmed true variant(s).
1.  **Rank Points:** If our top-ranked row (our #1 guess) exactly matches the true variant(s), we get 100 points. If the correct answer is lower in our list (e.g., #3 or #5), we get partial credit. If the true answer involves two variants (heterozygous pair) and our row only correctly guesses one of them, we get half credit.
2.  **F-Max:** The maximum F1-score across all possible confidence thresholds (`epcr`). It evaluates the best balance of **precision** (how many of our submitted guesses were actually right) and **recall** (did we find the true variant). The scoring system checks our ranked list and finds the "sweet spot" that maximizes this balance. This rewards pipelines that confidently isolate the right answer and heavily penalizes submissions that dump lots of noisy, low-confidence guesses. *Note: Secondary/incidental findings do not negatively impact the automated score.*

    *   **How F-max is calculated (Example):**
        *   The F1-score combines Precision and Recall. F-max calculates the F1-score at *every possible cutoff* of your `epcr` confidence values and keeps the highest score.
        *   **Team 1 (Good F-max):** Submits the true variant with `epcr: 0.99`, and 9 wrong variants with `epcr: 0.05`. At a threshold of 0.90, they have 1 guess, which is 100% accurate. Precision = 1.0, Recall = 1.0. **F-max = 1.0**.
        *   **Team 2 (Bad F-max):** Submits the true variant with `epcr: 0.99`, and 9 wrong variants with `epcr: 0.98`. At a threshold of 0.90, they have 10 guesses. Precision is only 10% (0.1). Because they couldn't mathematically separate the signal from the noise, their F-max drops significantly.

---

## Track 2 - Drug Repurposing
**Goal:** Identify existing, approved medications that are plausible candidates for further investigation based on the biological mechanisms disrupted in this case. Submissions are hypotheses for follow-up, not evidence that a medicine works.

Once we know (or have a strong guess at) what biological mechanisms are involved from Track 1, the next step is determining if an existing medication acts on this biological problem. We must characterize the mechanism behind the variant(s) and write up our reasoning and proposed candidate medication(s).

### The Input Data & Mechanism Characterization
The starting point for Track 2 is the output (the answer) from Track 1. Once we have confidently identified the specific causal genetic variant(s) and the affected gene, we use that knowledge as our input to find a potential treatment.

**What is Mechanism Characterization?**
Before we can pick a drug, we must understand *how* the genetic mutation broke the cell's biology. This step is the "Mechanism Characterization". It involves researching the gene we found in Track 1 to answer:
*   **Loss-of-Function (LOF):** Did the mutation destroy the gene, leaving the protein missing or broken? If so, we might need a drug that replaces its function, stabilizes the remaining broken protein, or bypasses the broken step in the pathway.
*   **Gain-of-Function (GOF):** Did the mutation make the protein hyperactive, or cause it to do something toxic? If so, we might need a drug (like an inhibitor) to block it or turn it off.
*   **Pathways & Downstream Consequences:** What biological cascade does this gene normally control? (For MVA, this is often the "Spindle Assembly Checkpoint" or cell division). If Gene A is broken, does it cause Protein B to accumulate dangerously? We map this out to find "druggable targets"—often proteins *downstream* of the broken gene that a medication can actually interact with.

Once this mechanism is clearly mapped out using medical literature (e.g., PubMed), we search pharmacological databases (like DrugBank or ChEMBL) for existing, FDA-approved medications that interact with this specific pathway to correct the downstream consequences.

### The Output (Deliverables)
We are permitted up to **3 submissions**, but **only the latest entry is reviewed** by the panel. *(Note: Only one designated team member should submit on behalf of the team; duplicate submissions will not be reviewed.)*

Each submission requires three components:

#### 1. Detailed Report (PDF or Markdown)
A comprehensive written report proposing repositioned drug candidates. 
*   **Mechanism Characterization:** The report *must* characterize the variant's mechanism (e.g., loss-of-function / gain-of-function, pathway disrupted, downstream biological consequence) to form the basis of the repurposing rationale.
*   **AI Usage Disclosure:** If an LLM or AI assistant is used, we must record the provider, tier, and data-handling settings (e.g., "Anthropic API, Claude Sonnet, commercial terms, no training on customer content").
*   **Format:** We should use the provided methods description template (`output/track2/methods_description_form.xlsx`) to organize the analysis.
*   **File Naming:** `<username>_track2_report.md` (or `.pdf`)

#### 2. Code Repository
This GitHub repository containing any computational methods or scripts used to find/validate the drug candidates.
*   Must contain documented, reproducible code and configuration files.
*   May remain private during the hackathon, but must be made public once it ends for judges to review.

#### 3. Pitch Video
A 3-minute pitch video walking through our reasoning and proposed candidate(s), uploaded to YouTube or Vimeo.

### Evaluation & Scoring
Unlike Track 1, Track 2 uses **qualitative evaluation** rather than automated scoring. An independent panel of expert judges will review all entries over a ~2-3 month window. Submissions are evaluated on scientific merit and real-world potential across four weighted criteria:

1.  **Scientific Rigor (35%):** Is the variant mechanism characterization sound? Is the proposed drug repurposing candidate well-supported by that mechanism?
2.  **Potential Impact (25%):** If independently validated, how much could this contribute to understanding or diagnosing MVA, for this child or for others living with it?
3.  **Innovation (25%):** Did you bring a genuinely creative angle, method, or tool to the problem?
4.  **Scalability (15%):** Could your approach realistically be applied to other individuals, other diseases, or larger datasets?

---

## Glossary & Notes
* **Phenotype:** The observable, physical expression of a disease or trait (the actual symptoms the doctor sees, like microcephaly or growth restriction).
* **Genotype:** The underlying genetic code or DNA sequence of the individual.
* **Homozygous:** Inheriting the exact same genetic variant from both parents (both copies of a gene have the identical mutation).
* **Heterozygous:** Inheriting a mutated variant from one parent, and a normal, healthy variant from the other.
* **Compound Heterozygous:** Inheriting two *different* mutated variants in the *same* gene (one broken copy from the mother, and a differently broken copy from the father). Both copies of the gene are broken, but in different ways. The challenge asks us to pair these up in a single row in our CSV if we suspect this is the cause.
* **FASTQ:** The raw, uninterpreted output straight from the DNA sequencing machine. It consists of billions of tiny puzzle pieces (short strings of DNA letters like A, C, T, G) with no context of where they belong in the human genome.
* **VCF (Variant Call Format):** How we find the differences. Bioinformaticians use massive computing power to take the billions of FASTQ "puzzle pieces" and match (align) them to a standard human reference blueprint. Once aligned, the computer compares the child's DNA to the reference and throws away everything that is identical, recording *only the differences*. A VCF file lists these differences (e.g., "At chromosome 1, position 100,000, a normal reference has an 'A', but this child has a 'T'").

---

## Project Workflow & Guidelines

To maintain organization and reproducibility, please follow this step-by-step workflow when contributing to the project:

### 1. Planning & Research
* **Procedure (`docs/01-procedure.md`)**: Start here. This file contains the high-level overview of our action plan. As you work, continuously update it with implementation details, further info, and results.
* **References (`docs/00-references.md`)**: Whenever you encounter interesting papers, literature, or external information, jot them down here to maintain a centralized knowledge base.

### 2. Experimentation & Discovery
* **Experiments (`notebooks/`)**: Conduct all exploratory and versioned experiments in this directory. Ensure these experiments follow the structured steps outlined in `docs/01-procedure.md`.
* **Findings (`docs/02-findings.md`)**: Whenever an experiment yields an interesting result, breakthrough, or insight, record it here immediately so the team can track discoveries over time.

### 3. Consolidation & Testing
* **Finalized Code (`src/`)**: Once a notebook experiment is complete and the optimal approach is identified, extract the best version of the code and move it here as reusable, production-ready modules.
* **Unit Tests (`tests/`)**: Every time a new piece of logic is added to the `src/` directory, you must write corresponding tests here to ensure code reliability.

### 4. Data & Configuration
* **Datasets (`data/`)**: All raw, processed, and downloaded data required for the project are stored here.
* **Agent Rules (`AGENTS.md`)**: Contains the core instructions, rules, and context for the LLM running in opencode. Refer to this for AI-specific operational guidelines.
