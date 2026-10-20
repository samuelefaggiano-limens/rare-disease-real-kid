# Methods Description — MVA Hackathon 2026, Track 1 (Variant Prediction)

*This document contains the required answers for the Track 1 submission, integrating the findings and methodology from our computational pipeline.*

### 1. Team name
Los Tanos

### 2. Model number
Model 1

### 3. Please describe your model/approach in detail.
Our computational approach consists of a 5-step automated pipeline designed for high-precision variant prioritization:

*   **Step 1 (Phenotype Extraction):** We extracted 8 core HPO terms from the clinical document, including HP:0002859 (Rhabdomyosarcoma) and HP:0200067 (Recurrent spontaneous abortion). Rhabdomyosarcoma is the strongest anchor: it is an established cancer predisposition in MVA. Recurrent pregnancy loss is suggestive of chromosomal instability but has no direct HPO annotation to any MVA gene.
*   **Step 2 (VCF Pre-Filtering):** To eliminate noise, we applied strict intrinsic QC filters (GQ>=20, DP>=10, 'PASS') and an Allelic Balance filter (0.30 - 0.70 for heterozygous variants). We then performed an allele-aware subtraction against dbSNP Common (MAF >= 1%). This step efficiently stripped 88.2% of the raw VCF, reducing the initial 5,012,204 variants to 591,339 ultra-rare and high-quality variants.
*   **Step 3 (Annotation):** We annotated the surviving variants using a locally built snpEff GRCh38.99 database, MyVariant.info (assembly=hg38), and AlphaGenome Atlas. This isolated 905 protein-altering (HIGH/MODERATE) variants, of which 193 have gnomAD AF ≤ 1% or no gnomAD record. Because dbSNP Common (2018) misses most common coding variants, gnomAD AF is the effective frequency filter; it is applied in Step 5.
*   **Step 4 (Phenotype Matching):** We applied an ontology-aware Resnik semantic similarity model to rank genes against the patient's HPO profile. *BUB1B* ranked #2 among candidate genes (Resnik similarity 0.632) and was the only Tier 1 MVA panel gene carrying rare HIGH/MODERATE variants.
*   **Step 5 (Ranking & EPCR Scoring):** We unified phenotype similarity, pathogenicity scores (CADD, REVEL, AVI_PHRED, ClinVar), allele frequencies, and inheritance modeling into an Estimated Probability of Causal Relationship (EPCR). The pipeline successfully isolated a compound-heterozygous pair in *BUB1B* (a stop-gain and a missense variant) as the top candidate. It achieved an EPCR of 0.9379, a clear separation from the second-ranked candidate (LZTR1, EPCR=0.5627).

### 4. (required) If commercially-available Generative AI was used, record the provider, the plan or tier, and the relevant data-handling setting.
- Google Cloud Vertex AI (Gemini API), training on our data disabled by Google's data-use policy, no 👍/👎 feedback provided.
- Anthropic, Claude Team plan (Claude Code desktop app, Claude Opus), commercial terms, no training on customer content, no 👍/👎 feedback provided.

### 5. Please state whether the submission file is the automated output of your computational approach (preferred), or if it has undergone downstream manual review and curation.
The submission file is the unedited, automated output of our computational approach. The EPCR scoring and final ranking were generated programmatically. The manual checks described in Q6 were done after generation and did not change any score or rank.

### 6. Please describe any downstream manual review in detail.
The submission file was not edited. After generation, we reviewed the top finding and the filtering steps for plausibility:
* **Clinical Literature Curation:** We manually reviewed the primary *BUB1B* finding against established MVA syndrome literature. Rhabdomyosarcoma is a documented cancer in MVA1 patients, which is consistent with the computational phenotype score. Recurrent pregnancy loss is consistent with chromosomal instability but is not, by itself, specific to MVA.
* **Structural & Functional Interpretation:** We performed a manual review of the specific molecular impact of the top *BUB1B* compound heterozygous pair. We checked the predicted consequences (a TTA→TGA stop-gain in exon 17 creating a premature stop at residue 751, and an AAT→AAG missense change in exon 23 replacing asparagine 1016 with lysine; numbering follows the Ensembl release-99 canonical transcript used by snpEff, and corresponds to Leu737 and Asn1002 in UniProt O60566) and drew a figure of the two phase scenarios (`docs/figures/bub1b_variants_explained.html`).
* **Quality Control Auditing:** To ensure the integrity of our filtering, we manually spot-checked 300 common SNVs against a local dbSNP dataset via `tabix`. This manual audit revealed critical gaps in the 2018 dbSNP file, which is why gnomAD AF ≤ 1% is enforced downstream.


### 7. Please describe any public data used in detail (e.g., reference panels, population databases, literature).
Our approach utilized several public databases and reference panels to filter, annotate, and prioritize variants:

*   **Ensembl (Release 99):** GRCh38 primary assembly DNA, GTF, and FASTA files were used to locally build a snpEff database for variant effect prediction.

*   **dbSNP (b151, 2018):** The common_all_20180418.vcf.gz dataset was used in our pre-filtering step for allele-aware subtraction of common variants.

*   **gnomAD:** Used for population allele frequency (AF) filtering, allowing us to confidently exclude variants with AF > 1% across genomes.

*   **Human Phenotype Ontology (HPO):** Release 2026-09-01 (the ontology and `genes_to_phenotype.txt`) was used to compute the ontology-aware Resnik semantic similarity scores for gene prioritization. The 8 HPO terms themselves came pre-tabulated in the clinical document.

*   **MyVariant.info:** Queried via API (using assembly=hg38) to aggregate variant annotations.

*   **Pathogenicity Predictors:** AlphaGenome Atlas (for AVI_PHRED, AlphaMissense, and CACTUS conservation scores), CADD, and REVEL were used to assess the functional impact of coding variants.

*   **Clinical & Disease Databases:** ClinVar and OMIM were utilized for clinical significance and disease-gene associations (e.g., mapping BUB1B to OMIM 257300 / MVA syndrome 1).

*   **Literature:** Existing literature on Mosaic Variegated Aneuploidy (MVA) syndrome was used to contextualize the phenotypic findings, specifically the established susceptibility to childhood cancers like rhabdomyosarcoma.

### 8. Please describe any proprietary data used in detail.
None. Our approach strictly relied on publicly available databases and open-source reference data (such as gnomAD, ClinVar, dbSNP, and the Human Phenotype Ontology) for variant annotation and phenotype matching.

### 9. Is your approach able to output proposed pairs of compound heterozygous candidate variants, or only single candidate variants?
Both. Our approach outputs proposed pairs of compound heterozygous candidate variants and single candidate variants. Our pipeline structures the output to explicitly distinguish between the two based on biological impact.
A single candidate variant (such as our *LZTR1* finding) leaves the secondary variant columns blank. When two rare heterozygous variants affect the same gene, the pipeline pairs them as a candidate compound-heterozygous configuration and populates the secondary variant columns (`chrom_2`, `pos_2`, `ref_2`, `alt_2`).
This capability is clearly demonstrated by our top-ranked primary finding in the *BUB1B* gene. The system paired a stop-gain (position 40209701) and a missense variant (position 40220612). Phase cannot be determined from a single-sample VCF; if the two are *in trans*, the expected result is a strong reduction of BubR1 function (the stop-gain allele is likely removed by nonsense-mediated decay and the missense allele is probably hypomorphic), which is the established mechanism of MVA1.

### 10. How did you handle secondary or incidental findings, if any? (See the finding_type column in your submission.)
Our pipeline was specifically optimized to identify the genetic driver of the patient's primary clinical presentation. Because our EPCR calculation heavily relies on an ontology-aware Resnik semantic similarity model to match variants against the patient's specific HPO profile, we intentionally did not curate or flag secondary findings. Consequently, all variants in our top 10 submission are designated as `primary`.
Coming from a foundational background in biological sciences rather than specialized genomic bioinformatics, we approached this dataset with a profound sense of ethical responsibility. We firmly believe that computational power must be guided by strict clinical criteria. It is often tempting to leverage AI to mine the genome and report secondary or incidental findings merely to showcase technical capability. However, doing so without a direct clinical mandate often generates erroneous assumptions and harmful "noise" rather than actionable medical insight.
This commitment to responsible analysis also guided our technical decisions. For instance, we faced a methodological dilemma regarding non-coding regions. We understand that these are vital regulatory zones, yet our extensive literature review into standard diagnostic workflows revealed that they are often excluded in this specific type of investigation to manage false positives. By acknowledging our boundaries and rigorously researching established protocols to build our clinical criteria, we chose to align with standard practices rather than blindly processing data. Therefore, variants that might be pathogenic but have weak phenotype similarity (such as those in Tier 2 genes *SGO1* and *MAD1L1*) were naturally penalized by our scoring system, ensuring our output remains strictly focused on the patient's true diagnostic need.


### 11. Please provide an estimate of run time and cost for your approach (if possible).
With all reference data downloaded and the snpEff database built: Step 2 pre-filtering takes about 1.5 min on a 10-core laptop (about 7 min in the notebook on 2 cores); Step 3 annotation is bounded by MyVariant.info and AlphaGenome API response times (minutes); Steps 4–5 take under a minute. There is no cost beyond the hardware (at least 2 cores, 8 GB RAM, 8 GB free disk; snpEff runs with `-Xmx8g`) and free API access.


### 12. Provide a method abstract (up to 500 words). 
Our computational approach for Track 1 consists of a fully automated, five-step pipeline designed to prioritize causal genetic variants from whole-genome sequencing data by integrating strict quality control, multi-source annotation, and quantitative phenotype matching.

Methodologically, the pipeline begins by extracting eight core Human Phenotype Ontology (HPO) terms from the patient’s clinical documentation. We then perform rigorous VCF pre-filtering, applying intrinsic quality control thresholds (GQ>=20, DP>=10, and a 0.30–0.70 allelic balance for heterozygous calls) followed by an allele-aware subtraction of common variants using dbSNP (MAF >= 1%). This step efficiently strips 88.2% of the initial 5 million raw variants. Surviving ultra-rare variants are annotated using a locally built snpEff GRCh38.99 database, MyVariant.info (hg38), and AlphaGenome Atlas to isolate protein-altering (HIGH/MODERATE) changes. To evaluate clinical relevance, we apply an ontology-aware Resnik semantic similarity model to rank genes against the patient's specific HPO profile. Finally, we calculate an Estimated Probability of Causal Relationship (EPCR) by unifying phenotype similarity, pathogenicity scores (CADD, REVEL, AVI_PHRED, ClinVar), allele frequencies (gnomAD), and inheritance modeling.

Using this approach, we successfully isolated a compound heterozygous pair of variants in the BUB1B gene (a stop-gain and a missense mutation) as the primary finding. The model achieved a highly confident separation for this candidate, yielding an EPCR of 0.9379.

A major strength of our method is its end-to-end automation and strict reliance on publicly available databases, minimizing subjective curation bias. Furthermore, the dynamic EPCR scoring mechanism excels at recognizing compound heterozygous configurations, which suits autosomal recessive conditions.

The primary limitation of our approach is its dependency on single-sample VCFs, meaning the phase of compound heterozygous variants cannot be definitively confirmed computationally and requires familial segregation analysis. Additionally, the pipeline's reliance on HPO semantic similarity is constrained by the inherent incompleteness of existing gene-to-phenotype databases, which occasionally lack full coverage for recently established panel genes.

### 13. Please include strengths and limitations of your method as well as methodological detail
**Methodological Strengths:**
* **End-to-End Automation:** The pipeline minimizes subjective human curation bias by programmatically executing all steps, from VCF pre-filtering to the final EPCR ranking.
* **Compound Heterozygous Detection:** The method pairs rare heterozygous variants within a gene into candidate compound-heterozygous configurations (as in our *BUB1B* finding).
* **Multidimensional Scoring:** The EPCR metric effectively integrates multiple distinct dimensions of data, including pathogenicity scores (CADD, REVEL, AVI_PHRED), population frequencies (gnomAD), and ontology-aware phenotype similarity (Resnik).
* **Public Data Reliance:** The exclusive use of open-source databases ensures high reproducibility and accessibility.

**Methodological Limitations:**
* **Phase Ambiguity:** Because the pipeline operates on single-sample VCFs rather than familial trios, the phase of compound heterozygous variants cannot be definitively confirmed computationally; true *in trans* configuration requires downstream familial segregation analysis.
* **Ontology Gaps:** The phenotype matching step relies heavily on the completeness of existing HPO databases. As we noted during our analysis, database coverage for some MVA panel genes is partial, meaning novel or uncurated gene-to-phenotype associations might be artificially penalized.
