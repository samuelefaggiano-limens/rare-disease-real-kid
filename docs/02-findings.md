# Findings

This document tracks all new interesting findings and insights (Track 1, Track 2, or secondary), ensuring traceability to the specific notebooks where they were generated.

## Track 1: Phenotype Findings
**Source:** `track1_step1_exp01-hpo_extraction.ipynb`
*   The clinical document directly provided the following HPO terms that outline the patient's severe presentation:
    *   **HP:0002859**: Rhabdomyosarcoma (Cancer)
    *   **HP:0000121**: Nephrocalcinosis (Renal)
    *   **HP:0004322**: Short stature (Growth)
    *   **HP:0001508**: Failure to thrive (Growth)
    *   **HP:0003202**: Skeletal muscle atrophy (Muscle)
    *   **HP:0001622**: Premature birth
    *   **HP:0001518**: Small for gestational age
    *   **HP:0200067**: Recurrent spontaneous abortion (Family History)
*   **Insight:** The presence of *Recurrent spontaneous abortion* in the parents strongly hints at a heritable chromosomal instability issue or aneuploidy running in the family. Furthermore, the *Rhabdomyosarcoma* perfectly matches the literature (found in `00-references.md`) stating that MVA patients are highly susceptible to childhood cancers like Wilms' tumor and rhabdomyosarcoma.

## Track 1: VCF Pre-Filtering Findings
**Source:** `track1_step2_exp01-vcf_prefiltering.ipynb`
*   Applying intrinsic QC filters (GQ>=20, DP>=10, and 'PASS') coupled with a strict Allelic Balance filter (0.30 - 0.70 for heterozygous variants) successfully removes nearly half of the raw variants.
*   **Insight:** The raw VCF contained 5,012,204 variants. After basic intrinsic QC, 2,645,065 variants survived (a ~47.2% reduction in noise). These filtered variants are saved in `data/WGS_EX2312012_HGWCNDSX7_intrinsic_qc.vcf.gz`. 