import re

with open('/home/samuele/Projects/Learning/rare-disease-real-kid/docs/01-procedure.md', 'r', encoding='utf-8') as f:
    content = f.read()

start_marker = "### Step 3: VCF Annotation & Filtering"
end_marker = "### Step 4: Gene Prioritization & Phenotype Matching"

start_idx = content.find(start_marker)
end_idx = content.find(end_marker)

if start_idx == -1 or end_idx == -1:
    print("Could not find markers.")
    exit(1)

old_step3 = content[start_idx:end_idx]

# Function to safely strike through lines while keeping markdown structure
striked_lines = []
for line in old_step3.split('\n'):
    stripped = line.strip()
    if not stripped:
        striked_lines.append(line)
    elif stripped.startswith('#'):
        m = re.match(r'^(#+)\s+(.*)', line)
        if m:
            striked_lines.append(f"{m.group(1)} ~~{m.group(2)}~~")
        else:
            striked_lines.append(f"~~{line}~~")
    elif stripped.startswith('*'):
        m = re.match(r'^(\s*\*\s+)(.*)', line)
        if m:
            striked_lines.append(f"{m.group(1)}~~{m.group(2)}~~")
        else:
            striked_lines.append(f"~~{line}~~")
    elif stripped.startswith('|'):
        if '---' in line:
            striked_lines.append(line)
        else:
            # strike through each cell
            parts = line.split('|')
            new_parts = []
            for p in parts:
                if p.strip():
                    new_parts.append(f" ~~{p.strip()}~~ ")
                else:
                    new_parts.append(p)
            striked_lines.append('|'.join(new_parts))
    else:
        striked_lines.append(f"~~{line}~~")

striked_step3 = '\n'.join(striked_lines)

new_step3 = """### Step 3: VCF Annotation & Filtering (Updated)

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

---
"""

# Now write back
new_content = content[:start_idx] + new_step3 + "\n\n" + striked_step3 + "\n" + content[end_idx:]

with open('/home/samuele/Projects/Learning/rare-disease-real-kid/docs/01-procedure.md', 'w', encoding='utf-8') as f:
    f.write(new_content)

print("Procedure updated.")
