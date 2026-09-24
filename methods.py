import os

markdown_content = """# Methods Description — MVA Hackathon 2026, Track 1 (Variant Prediction)

*This document contains the required answers for the Track 1 submission, integrating the findings and methodology from our computational pipeline.*

### 1. Team name
Lara Pugnaloni

### 2. Model number 
Model 1 

### 3. Please describe your model/approach in detail.
Our computational approach consists of a 5-step automated pipeline designed for high-precision variant prioritization:

*   **Step 1 (Phenotype Extraction):** We extracted 8 core HPO terms from the clinical document, including HP:0002859 (Rhabdomyosarcoma) and HP:0200067 (Recurrent spontaneous abortion), which provided strong phenotypic signals matching MVA syndrome susceptibility literature.
*   **Step 2 (VCF Pre-Filtering):** To eliminate noise, we applied strict intrinsic QC filters (GQ>=20, DP>=10, 'PASS') and an Allelic Balance filter (0.30 - 0.70 for heterozygous variants). We then performed an allele-aware subtraction against dbSNP Common (MAF >= 1%). This step efficiently stripped 88.2% of the raw VCF, reducing the initial 5,012,204 variants to 591,339 ultra-rare and high-quality variants.
*   **Step 3 (Annotation):** We annotated the surviving variants using a locally built snpEff GRCh38.99 database, MyVariant.info (assembly=hg38), and AlphaGenome Atlas. This isolated 1,526 protein-altering (HIGH/MODERATE) variants. gnomAD AF was utilized as the effective frequency filter (max AF 1%).
*   **Step 4 (Phenotype Matching):** We applied an ontology-aware Resnik semantic similarity model to rank genes against the patient's HPO profile. *BUB1B* ranked #2 genome-wide (similarity score 0.632) and was uniquely identified as the only Tier 1 MVA panel gene carrying rare HIGH/MODERATE variants.
*   **Step 5 (Ranking & EPCR Scoring):** We unified phenotype similarity, pathogenicity scores (CADD, REVEL, AVI_PHRED, ClinVar), allele frequencies, and inheritance modeling into an Estimated Probability of Causal Relationship (EPCR). The pipeline successfully isolated a compound-heterozygous pair in *BUB1B* (a stop-gain and a missense variant) as the top candidate. It achieved an EPCR of 0.9379, showing a massive and confident separation from the second-ranked candidate (LZTR1, EPCR=0.5627).

### 4. (required) If commercially-available Generative AI was used, record the provider...
[Reemplazar con tu uso de IA, por ejemplo: "OpenAI API, GPT-4o, commercial terms, no training on customer content" o "No commercially-available Generative AI was used for the computational pipeline."]

### 5. Please state whether the submission file is the automated output of your computational approach (preferred), or if it has undergone downstream manual review and curation.
The submission file is the completely automated output of our computational approach. The EPCR scoring and final ranking (including the confident separation of the BUB1B compound-heterozygous variants) were generated programmatically without downstream manual review or expert curation.
"""

# Nombre del archivo a generar
file_name = "methods_description_submission.md"

# Crear y escribir el archivo
with open(file_name, "w", encoding="utf-8") as file:
    file.write(markdown_content)

print(f"✅ Archivo '{file_name}' generado exitosamente desde methods.py.")