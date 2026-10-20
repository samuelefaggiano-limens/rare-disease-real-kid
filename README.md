# Rare Disease, Real Kid: MVA Hackathon 2026 — team Los Tanos

Our submission to the [SageBio "Rare Disease, Real Kid" MVA Hackathon 2026](https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026). The genome and clinical story belong to a real child living with Mosaic Variegated Aneuploidy (MVA), an ultra-rare condition affecting fewer than 50 people worldwide, with no established treatment.

**Team:** Samuele Faggiano, Lara Pugnaloni.

> Everything here is a research hypothesis for follow-up, not medical advice. No drug proposed in this repository has been tested for MVA.

---

## Results at a glance

### Track 1 — Variant prediction
From 5,012,204 variants in the child's whole-genome VCF, the pipeline ranks a **compound-heterozygous pair in *BUB1B*** (the MVA1 gene) first:

| | Variant | Consequence |
|---|---|---|
| 1 | p.Leu751\* (stop-gain, exon 17) | premature stop; transcript predicted to undergo nonsense-mediated decay |
| 2 | p.Asn1016Lys (missense, exon 23) | change in the C-terminal pseudokinase domain |

`epcr` = **0.9379**, against 0.5627 for the next candidate (*LZTR1*). Phase is not determinable from a single-sample VCF; the pair is consistent with a compound-heterozygous configuration. Protein numbering follows the Ensembl release-99 canonical transcript used by snpEff (Leu737 and Asn1002 in UniProt O60566).

### Track 2 — Drug repurposing
- **Mechanism:** partial loss of function. Reduced BubR1 weakens the spindle assembly checkpoint, which causes mosaic aneuploidy, cancer predisposition and an accumulation of senescent cells.
- **Central finding:** no approved drug acts on any protein at STRING distance 1 from BubR1. Running the same lookup automatically for all eight genes with an MVA or MVA-like phenotype shows the gap holds across MVA: 2 of 108 high-confidence interactors have any approved drug, and neither acts on chromosome segregation.
- **Prioritised hypotheses:**
  - **H1 (0.817):** NAD+ precursors → SIRT2 → BubR1 stabilisation, framed as remission maintenance.
  - **H2 (0.667):** targeting BubR1-driven senescent cells with low-dose rapalogs as senomorphics, with senolytics as a comparator.
- Both hypotheses fail or stand on the same untested link, and one patient-cell experiment tests it (report §7).

---

## Deliverables

| Track | File |
|---|---|
| 1 | Submission CSV: [`output/track1/lostanos_hpo_matcher_v1.csv`](output/track1/lostanos_hpo_matcher_v1.csv) |
| 1 | Methods description: [`output/track1/methods_description_submission.md`](output/track1/methods_description_submission.md) · report: [`output/track1/lostanos_track1_report.pdf`](output/track1/lostanos_track1_report.pdf) |
| 2 | Report: [`output/track2/lostanos_track2_report.md`](output/track2/lostanos_track2_report.md) |
| 2 | Methods description (organisers' template): [`output/track2/lostanos_track2_methods_description.xlsx`](output/track2/lostanos_track2_methods_description.xlsx) |
| 2 | Pitch video: *link added on upload* · script: [`output/track2/lostanos_track2_pitch_script.md`](output/track2/lostanos_track2_pitch_script.md) |
| both | Figures for non-specialists (EN/ES): [`docs/figures/`](docs/figures/) |

## Repository layout

```
docs/
  00-references.md     every paper, database and tool consulted
  01-procedure.md      the plan, step by step, with an Experiments table per step
  02-findings.md       insights as they were found, attributed to the notebook/script that produced them
  figures/             build_*.py -> bilingual HTML + PNG
notebooks/
  track{N}_step{S}_exp{NN}-{slug}.ipynb   versioned experiments (Track 1, steps 0-5)
  utils/               consolidated, runnable pipeline scripts (Track 1 steps 2-5, Track 2 steps 3, 4, 6)
output/
  track1/              submission CSV, methods description, report; templates/ holds the organisers' blanks
  track2/              report, methods description, pitch script, Step 3/4/6 result tables
requirements.txt
```

## How to reproduce

### Setup
```bash
pip install -r requirements.txt
```
Python 3.12. `bcftools` and `tabix` are not needed: VCF I/O is pure Python (`cyvcf2`, `gzip`).

### Track 2 (public data only, no patient data needed)
Run from the repository root:
```bash
python notebooks/utils/track2_step3_drug_mapping.py            # ~5 min; live ChEMBL, DGIdb, STRING, UniProt queries
python notebooks/utils/track2_step4_prioritize.py              # instant, offline; scores the hypotheses
python notebooks/utils/track2_step6_mva_druggability_scan.py   # ~4 min; same lookup for all 8 MVA genes
python docs/figures/build_track2_mechanisms_figure.py          # regenerates the mechanisms figure
```
Outputs go to `output/track2/`. Database contents change over time, so a re-run on a later date can differ slightly from the committed tables (queried 2026-10-14).

### Track 1 (requires the gated patient data)
1. **Get the data.** The VCF and clinical document are gated: request access on the hackathon's Hugging Face Space and accept the data-use terms. Put the files in `data/` at the repository root (`data/` is gitignored). Also place there:
   - dbSNP Common `common_all_20180418.vcf.gz` (NCBI, GRCh38);
   - HPO `hp.obo` and `genes_to_phenotype.txt` (release 2026-09-01).
2. **snpEff.** Install snpEff and Java 17. The `GRCh38.99` database cannot be downloaded with the current snpEff release, so build it from Ensembl release 99 with `snpEff build -gtf22 GRCh38.99` (details in [`docs/02-findings.md`](docs/02-findings.md), "Environment & Reference-Data Findings").
3. **Run the steps** from the repository root. Each step reads the previous step's output from `data/` by default:
   ```bash
   python notebooks/utils/run_step2_prefilter.py           # QC, allelic balance, dbSNP Common subtraction -> 591,339 variants
   java -Xmx8g -jar snpEff/snpEff.jar -canon GRCh38.99 data/track1_step2_rare_variants.vcf.gz > data/snpeff_annotated.vcf
   python notebooks/utils/track1_step3a_annotate.py        # HIGH/MODERATE filter + MyVariant.info (gnomAD, CADD, REVEL, ClinVar)
   python notebooks/utils/step3b_alphagenome.py            # AlphaGenome Atlas scores; needs ALPHAGENOME_API_KEY in the environment or .env
   python notebooks/utils/track1_step4_prioritize.py       # HPO phenotype matching (Resnik) + MVA panel tiers
   python notebooks/utils/track1_step5_submission.py       # epcr ranking -> output/track1/lostanos_hpo_matcher_v1.csv
   ```
   The same steps are available as notebooks in `notebooks/`, run with `notebooks/` as the working directory. The submitted CSV came from the notebook run of Step 3 (905 HIGH/MODERATE variants); a separate script run of Step 3 kept 1,526, and the cause of the difference is recorded as open in [`docs/02-findings.md`](docs/02-findings.md).

Filtering thresholds (GQ ≥ 20, DP ≥ 10, allelic balance 0.30–0.70 for heterozygous calls, gnomAD AF ≤ 1 %) and the rationale for each are in [`docs/01-procedure.md`](docs/01-procedure.md).

## Data privacy

The child's genome and clinical record are not in this repository. The VCF, FASTQ, clinical document and every intermediate file carrying genotypes stay in the gitignored `data/` folder and are deleted at the end of the hackathon, as the data-use agreement requires. The only patient-derived content committed is what the deliverables themselves need: the ten ranked variants of the Track 1 submission and the named *BUB1B* variants discussed in the reports.

## AI usage

AI assistants were used for coding, literature search and writing. The provider, plan and data-handling setting for each tool are recorded in the methods descriptions of both tracks, as the organisers require.

## License

The organisers release participant submissions, code and reports under [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). This repository is released under the same license; see [`LICENSE`](LICENSE).
