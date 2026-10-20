# Drug Repurposing Hypotheses for BUB1B-Related Mosaic Variegated Aneuploidy

**MVA Hackathon 2026 — Track 2 (Drug Repurposing)**
**Team:** Los Tanos · **Report status:** final draft (2026-10-20)
**Code:** https://github.com/samuelefaggiano-limens/rare-disease-real-kid

> These are hypotheses for experimental follow-up, not evidence that any medicine works, and not clinical recommendations. Any decision about this child's care belongs to the treating team.

---

## Summary

The child carries two rare heterozygous *BUB1B* variants: a stop-gain (p.Leu751\*) and a missense change (p.Asn1016Lys). Together they fit the compound-heterozygous pattern of MVA type 1. We think the stop-gain allele is a null allele (removed by nonsense-mediated decay). The missense allele is the child's only source of BubR1 protein, and the change sits in the C-terminal pseudokinase domain, which is needed to recruit the PP2A-B56 phosphatase to kinetochores. The net expected effect is **partial loss of BubR1 function**. That weakens the spindle assembly checkpoint (SAC) and leads to chromosome mis-segregation and constitutional mosaic aneuploidy.

We mapped this mechanism to a panel of 35 candidate proteins and queried ChEMBL, DGIdb and STRING for approved drugs.

**Central finding:** no approved drug acts in the desired direction on any of the 18 proteins that interact directly with BubR1 (STRING distance 1). All 26 approved, direction-concordant drugs act on proteins further away (distance 2–3) or with no network path to BUB1B. No approved drug corrects the defect itself. The useful question becomes which *consequence* of the defect is worth addressing, and how honest the evidence for addressing it is.

We scored the surviving routes with a transparent, weighted framework in which mechanistic relevance and direction of effect dominate. Two hypotheses result:

| | Hypothesis | Concept | Score |
|---|---|---|---|
| **H1** | **NAD+ precursors (niacin, nicotinamide riboside) → SIRT2 → BubR1 stabilisation** | functional compensation | 0.817 |
| **H2** | **Targeting BubR1-driven senescent cells (low-dose rapalogs as senomorphics; senolytics as comparator)** | phenotype modulation | 0.667 |

H1 is the only route we found that acts on the amount of BubR1 rather than on its consequences. It carries an open contradiction that we state up front: rhabdomyosarcoma cells depend on BUB1B, and this child has had a rhabdomyosarcoma. H2 has the most direct match between model and phenotype in the analysis: senescent-cell clearance delayed muscle wasting in BubR1-deficient mice, and this child has muscle atrophy. For each hypothesis we propose an experiment that distinguishes phenotypic rescue from selective toxicity or proliferative side effects.

Both hypotheses carry a tumour-related drawback, and we examined how to mitigate each one (§5):
- **H1:** the mitigation is *timing*, not a second drug. Co-administering an NAD+- or BubR1-lowering agent would cancel the rescue itself. H1 is therefore framed as a remission-maintenance hypothesis, supported by nicotinamide riboside safety data from two other genome-instability syndromes.
- **H2:** the mitigation is a *different drug class within the same node*. Low-dose rapalogs act as **senomorphics**: they suppress the inflammatory secretome of senescent cells without killing them, keep the anti-tumour arrest in place, and have paediatric rhabdomyosarcoma exposure data.

**Beyond this child.** Running the same lookup automatically for all eight genes with a reported MVA phenotype shows the druggability gap is not specific to *BUB1B*: only 2 of 108 high-confidence first-shell interactors have any approved drug, and neither acts on the mitotic defect (§9). The method itself is gene-agnostic and ran in about four minutes. §10 sets out what the proposed experiments would contribute to understanding and diagnosing MVA, including functional classification of the child's unclassified missense variant.

---

## 1. Starting point

**From Track 1.** Our variant-prioritisation pipeline ranked a *BUB1B* compound-heterozygous pair first, with an epcr of 0.94; the next candidate scored 0.56. *BUB1B* is a Tier 1 MVA gene (MVA1, OMIM 257300). No other established MVA gene (*CEP57, TRIP13, MAD2L1BP, CENATAC, CEP192, SMC5, BUB1*) carried rare protein-altering variants.

**The child's phenotype** (from the clinical document, mapped to HPO):
- Rhabdomyosarcoma.
- Nephrocalcinosis.
- Short stature, failure to thrive, small for gestational age, premature birth.
- Skeletal muscle atrophy.
- Recurrent spontaneous abortion in the parents.

**No microcephaly and no congenital heart defect are reported.** This matters below: two well-evidenced MVA mechanisms target exactly those two features, and we set them aside for this child.

---

## 2. Mechanism characterisation

### 2.1 Allele 1: p.Leu751\* (stop-gain, exon 17 of 23)

| Evidence type | Observation |
|---|---|
| Variant class / transcript | A TTA→TGA change creates a premature stop at residue 751 of 1,064. The stop lies well upstream of the last exon–exon junction, so the transcript is predicted to undergo **nonsense-mediated decay (NMD)**. |
| Population | gnomAD genomes AF ≈ 3×10⁻⁵ |
| Clinical | ClinVar: Pathogenic/Likely pathogenic (single submitter) |
| Computational | CADD 36; AlphaGenome AVI Phred 33.8 (top 0.1 %), dominant modality PROTEIN_TERMINATION |
| Literature | In PCS/MVA families, transcripts from truncating and splice-site *BUB1B* alleles were significantly reduced, consistent with NMD (Matsuura 2006, PMID 16411201) |

**Interpretation: a null allele.** If NMD escape did occur, the truncated product (residues 1–750) would keep the N-terminal TPR, KEN boxes and KARD, but would lack the pseudokinase domain. We treat NMD as the expected outcome, not a demonstrated one.

### 2.2 Allele 2: p.Asn1016Lys (missense, exon 23)

The procedure requires that we do not infer function from variant class, and that each type of evidence is kept separate:

| Evidence type | Observation | Strength |
|---|---|---|
| Population | Absent from gnomAD | supports rarity only |
| Clinical | No ClinVar record | none |
| Domain | Lies in the C-terminal **pseudokinase domain**. This domain promotes kinetochore recruitment of PP2A-B56, SAC silencing and chromosome alignment (Gama Braga 2020, PMID 33207204). | moderate (position, not tested for this residue) |
| Computational | CADD 24.5; REVEL 0.47 (borderline); AlphaGenome AVI Phred 25.6, dominant modality ALPHAMISSENSE | weak–moderate, discordant across predictors |
| Functional assay | **None exists for N1016K** | — |

**Interpretation: probably hypomorphic, not proven.** The functional effect of N1016K is the least certain link in the chain, and every hypothesis below inherits that uncertainty. Asparagine to lysine adds a positive charge inside a folded domain. That is plausible as a destabilising change, but it has not been demonstrated.

### 2.3 Combined effect and phase

The two variants are heterozygous with allelic balance 0.54 and 0.46. **Phase cannot be determined** from a single-sample VCF without parental DNA. If the two variants are in trans, the child's BubR1 comes entirely from the N1016K allele, at roughly half dose or less. That matches the MVA1 model: one null allele plus one hypomorphic allele, consistent with the literature on BubR1 dosage (§2.4).

**Mechanistic chain:**

```
BUB1B L751* (NMD, null)  +  N1016K (pseudokinase domain, hypomorphic?)
        → reduced BubR1 quantity, and impaired PP2A-B56 recruitment by what remains
        → weakened spindle assembly checkpoint + unstable kinetochore–microtubule attachment
        → premature anaphase, chromosome mis-segregation
        → constitutional mosaic aneuploidy  → growth restriction, cancer predisposition (rhabdomyosarcoma)
        ↘ (separate arm) accumulation of p16-positive senescent cells → tissue degeneration (muscle)
```

### 2.4 Why BubR1 dosage is the key variable

- BubR1-hypomorphic mice develop progressive aneuploidy, dwarfism, muscle wasting and short lifespan (Baker 2004, PMID 15208629).
- Mice that over-express BubR1 are protected from aneuploidy and cancer, and keep better muscle function (Baker 2013, PMID 23242215).
- In patients, a decrease of more than 50 % in BubR1 expression or activity was proposed as sufficient for the PCS/MVA phenotype. Restoring normal BubR1 levels by chromosome-15 transfer corrected the checkpoint defect in patient cells (Matsuura 2006, PMID 16411201).

**BubR1 behaves as a dosage-sensitive protein across its whole range, not as an on/off switch.** This is what makes "raise the BubR1 that remains" a coherent therapeutic goal (H1).

---

## 3. From mechanism to targets

### 3.1 Approach

We built a **35-protein target panel** (`notebooks/utils/track2_target_panel.tsv`). Each protein was annotated with:
- its relationship to BUB1B: direct interactor, same pathway, MVA gene, or indirect chromosome-instability (CIN) biology;
- the **desired direction of modulation** (increase / decrease / context-dependent);
- the therapeutic concept it serves: functional compensation, downstream correction, or vulnerability/phenotype modulation;
- the part of *this child's* phenotype it addresses.

Keeping the panel in a TSV rather than in code lets anyone review or dispute an inclusion without touching the pipeline.

Each gene was resolved to its reviewed UniProt accession. Approved drugs were then mapped through **ChEMBL** (primary source; curated mechanisms of action), cross-checked against **DGIdb**, and network distance to BUB1B was taken from **STRING** as supporting evidence only. Two implementation details change the result and are documented in the code:
1. **ChEMBL has to be queried through every target that contains the protein**, not only single-protein targets. Bortezomib, for example, is annotated against "26S proteasome", not against PSMB5.
2. **Salts and hydrates are collapsed to the parent molecule.** Without this, approval phase is under-reported.

ChEMBL annotates rapalogs against FKBP1A rather than MTOR, so we added FKBP1A to the panel explicitly.

### 3.2 Central finding: proximity and druggability are inversely related

| Relationship to BUB1B | Genes | Genes with an approved, direction-concordant drug |
|---|---|---|
| Direct interactor | 10 | **0** |
| Same pathway | 7 | **0** |
| MVA gene (BUB1B, TRIP13) | 2 | **0** |
| Indirect CIN biology | 16 | 8 |

All 18 proteins at STRING distance 1 from BUB1B have **no** approved drug acting in the desired direction: CDC20, TRIP13, TTK, AURKB, PLK1, BUB1, BUB3, MAD2L1, CENPE, NDC80, SGO1, ESPL1, UBE2C, ANAPC2, MAD1L1, KIF18A, PPP2CA, FOXM1.

The 26 approved concordant drugs act on:
- **distance 2:** PSMB5, BCL2, SRC, MTOR;
- **distance 3:** ABL1, FKBP1A;
- **no STRING path to BUB1B:** DHODH, POLE.

Mitotic machinery is essentially undrugged in the approved pharmacopoeia, because inhibiting it is only ever desirable in oncology.

**Direction matters as much as availability.** PLK1 and CENPE inhibitors appear prominently in the network. But PLK1 phosphorylates BubR1's KARD to recruit PP2A-B56, which is exactly the function the N1016K allele is predicted to weaken. Inhibiting PLK1 would remove the recruitment that remains. **A drug that exists but pushes in the wrong direction is not a candidate.**

### 3.3 Routes examined and why each was retained or set aside

| Reason | Routes |
|---|---|
| **No approved drug exists** | APC/C–CDC20 inhibition; TRIP13 inhibition |
| **Wrong direction** | PLK1 / CENPE inhibition; spindle poisons as a therapy |
| **Phenotype not present in this child** | Wnt inhibition (congenital heart defect, PMID 40055864); neural stem cell protection (microcephaly, PMID 41820377) |
| **Published negative result** | P53 blockade (PMID 37900274); growth hormone (CEP57-MVA, PMID 35434947) |
| **Toxicity, not rescue** | MPS1 + KIF18A synthetic lethality (requires a competent SAC, PMID 42647139); targeting the aneuploid state |
| **Retained** | NAD+/SIRT2 → BubR1 (H1); senescent-cell clearance (H2) |

The two hypotheses we kept come from arms of BubR1 biology that do not run through mitosis: protein **dosage** and **senescence**. Our initial target network was built around the mitotic machinery and missed both. They entered the analysis during drug mapping, when it became clear that nothing proximal to the mitotic defect was druggable.

---

## 4. Prioritisation framework

Each route was scored from 0 to 3 on the five axes the procedure specifies, then combined as a weighted sum (`notebooks/utils/track2_step4_prioritize.py`):

| Axis | Weight | Question |
|---|---|---|
| Mechanistic relevance | 0.35 | How close is the target to the BUB1B defect? |
| Directionality | 0.20 | Is there evidence the drug acts in the desired direction? |
| Experimental evidence | 0.20 | BUB1B-deficient models > SAC/CIN models > computational |
| Human exposure | 0.15 | Clinical use, pharmacology, safety, paediatric data |
| Tractability | 0.10 | Approved, route, formulation, availability |

Mechanistic relevance and directionality together carry 55 % of the weight. **Drug availability on its own cannot carry a hypothesis.** An axis that does not apply is dropped from the denominator rather than scored zero. No model was trained, because there are no labelled outcomes for MVA repurposing. The weights are published so a reader can disagree and recompute.

| Rank | Hypothesis | Score | Mech. | Dir. | Exp. | Human | Tract. |
|---|---|---|---|---|---|---|---|
| 1 | **H1 NAD+ → SIRT2 → BubR1** | **0.817** | 3 | 2 | 2 | 2 | 3 |
| 2 | **H2 Senescent-cell clearance** | **0.667** | 2 | 2 | 2 | 2 | 2 |
| 3 | H5 Targeting the aneuploid state | 0.583 | 1 | 2 | 2 | 2 | 3 |
| 4 | H4 APC/C–CDC20 braking | 0.483 | 3 | 1 | 1 | 0 | 0 |
| — | H3 Spindle-poison interaction | *(0.854)* | 3 | n/a | 2 | 2 | 3 |

H3 has the highest number but is **excluded by design**: it is a research question, not a repurposing proposal (§6.1).

---

## 5. Prioritised hypotheses

### 5.1 H1 — NAD+ precursors to stabilise the remaining BubR1

1. **BUB1B defect and affected process.** BubR1 quantity is reduced: one allele is null, and the other probably produces a hypomorphic protein. The result is SAC weakness and mis-segregation.
2. **Target and drug mechanism.**
   - SIRT2 keeps BubR1 deacetylated at lysine 668 and protects it from degradation. The acetyltransferase CBP opposes this.
   - SIRT2 activity depends on NAD+.
   - NAD+ precursors raise cellular NAD+. Candidates: **niacin** (an approved drug), **nicotinamide riboside** and **nicotinamide mononucleotide** (sold as supplements).
3. **Rationale for the direction.** Raising BubR1 is protective in mice (PMID 23242215). NMN treatment raised BubR1 *in vivo*, and SIRT2 overexpression extended lifespan in BubR1-hypomorphic mice (North 2014, PMID 24825348). K668 lies outside both patient variants. Because the stop-gain allele is expected to undergo NMD, the N1016K-derived protein is the plausible substrate. **This is an inference from domain position, not an experimental result for N1016K.**
4. **Supporting evidence and limitations.**
   - All BubR1 evidence is murine. There is no MVA-specific test.
   - Stabilising a hypomorphic protein helps only if that protein keeps partial function, which is unknown for N1016K.
   - **Open contradiction:** rhabdomyosarcoma cells depend on the FoxM1/BUB1B axis for growth (Wan 2012, PMID 23002205), and this child has had a rhabdomyosarcoma. NAD+ availability also supports proliferative metabolism. Raising BubR1 systemically in a cancer-predisposition background could stabilise mitosis in healthy tissue *or* support tumour-cell proliferation. We record this as an unresolved contradiction that could invalidate the hypothesis, not as a caveat.
   - **Why a co-administered "safety" drug does not resolve it.** Tumour cells depend on the NAD+ pool, and NAD+ depletion (e.g. NAMPT inhibition) is itself studied as an antitumour strategy. Supplying NAD+ precursors restored tumourigenic properties to NAMPT-inhibited cancer cells (Navas 2023, PMID 36864434). Adding an NAD+- or BubR1/FoxM1-lowering drug would cancel the rescue in exactly the healthy tissue H1 is meant to protect, and would add systemic toxicity in a child with nephrocalcinosis and failure to thrive.
   - **Mitigation: timing.** We therefore frame H1 as a **remission-maintenance hypothesis**. It is to be considered only after documented complete oncological remission, never during active tumour burden, and only after the paired go/no-go cell screen below. This does not remove the contradiction; it turns it into an explicit precondition plus a laboratory decision point.
5. **Human and paediatric safety.**
   - Nicotinamide riboside was given to 24 people with ataxia-telangiectasia (17 under 18), another genome-instability syndrome with cancer predisposition. It was tolerated, and ataxia scores and immunoglobulin levels improved (Veenhuis 2021, PMID 34515380).
   - A 2-year open-label extension in ataxia-telangiectasia patients of all ages found NR safe and well tolerated, with stable NAD+ elevation in blood (Presterud 2024, PMID 37899683).
   - In Werner syndrome, another constitutional genome-instability disorder, a randomised placebo-controlled crossover trial of NR (1,000 mg/day, 26 + 26 weeks) reported no serious adverse events and a significant decrease in blood creatinine (Shoji 2025, PMID 40459998). A renal-function signal is directly relevant to a child with nephrocalcinosis, but it was a secondary endpoint in a different disease.
   - Oral NMN was tolerated in healthy adults (Irie 2020, PMID 31685720).
   - Niacin corrected NAD+ deficiency and improved muscle performance in adult mitochondrial myopathy (Pirinen 2020, PMID 32386566). It has well-known dose-dependent adverse effects (flushing; hepatotoxicity with high-dose or sustained-release forms).
   - None of these studies is in MVA, none tested the BubR1 mechanism, and none was powered for tumour incidence. They establish tolerability and a monitoring precedent, not efficacy.
6. **Proposed experiment.** Use patient-derived cells (fibroblasts or iPSC) and isogenic *BUB1B* L751\*/N1016K lines against controls, with vehicle versus nicotinamide riboside.
   - **Readouts:** BubR1 protein (immunoblot), K668 acetylation (IP–MS), micronuclei, mis-segregation by live imaging, aneuploidy by single-cell karyotyping.
   - **Run the same exposure in matched rhabdomyosarcoma lines** to test whether it increases proliferation. This separates phenotypic rescue from a pro-tumour effect.

### 5.2 H2 — Targeting BubR1-driven senescent cells

1. **BUB1B defect and affected process.** Low BubR1 drives accumulation of p16-positive senescent cells. This arm is independent of mis-segregation and contributes to tissue degeneration, including sarcopenia.
2. **Target and drug mechanism.** Senolytics kill senescent cells by disabling the pro-survival network those cells depend on (SRC/ABL family kinases, BCL-xL, PI3K nodes; Zhu 2015, PMID 25754370). The best-studied combination is **dasatinib** (an approved drug) **+ quercetin** (a supplement); fisetin is an alternative.
3. **Rationale for the direction.** Genetic clearance of p16-positive cells in BubR1-hypomorphic mice delayed sarcopenia, fat loss and cataracts. Late clearance still slowed progression (Baker 2011, PMID 22048312). **This experiment was done in the BubR1 model itself, and the phenotype it rescued (muscle wasting) is documented repeatedly in this child.**
4. **Supporting evidence and limitations.**
   - The rescue was *genetic*, not pharmacological, and not in MVA.
   - Senescence is a tumour-suppressive programme, so removing senescent cells in a cancer-predisposition background cuts both ways.
   - This route addresses a consequence; it does not correct the BubR1 defect.
   - **Mitigation: senomorphics instead of senolytics.** Senomorphics suppress senescence markers, above all the senescence-associated secretory phenotype (SASP), without killing the cell (Zhang 2023, PMID 35015337). They keep the proliferative arrest, which is itself a barrier to malignant progression (Di Micco 2021, PMID 33328614; Colucci 2025, PMID 40513577), while silencing the IL-1α/IL-6/NF-κB secretome that drives both tissue wasting and paracrine tumour promotion.
     - **Rapalogs (sirolimus, everolimus, temsirolimus) are the strongest option.** mTOR inhibition blunts the SASP by suppressing IL1A translation (Laberge 2015, PMID 26147250) and via MAPKAPK2 translation control (Herranz 2015, PMID 26280535). In mice, rapamycin also **suppressed the ability of senescent cells to stimulate tumour growth** (Laberge 2015), which is the opposite of the senolytic risk.
     - The rapalogs entered our drug mapping through FKBP1A, which ChEMBL uses as their annotated target. Once FKBP1A was resolved, the Step 3 re-run lists all three as approved, direction-concordant drugs.
     - **Two other approved senomorphics, both deprioritised for this child.** Metformin blocks SASP expression via IKK/NF-κB (Moiseeva 2013, PMID 23521863), but carries a black-box warning for lactic acidosis, which is a concern given the renal background (nephrocalcinosis). Ruxolitinib (JAK1/2) alleviates the SASP and frailty in aged mice (Xu 2015, PMID 26578790), but its black-box warning (infection, malignancy, cytopenias) weighs against tumour surveillance.
     - **No senomorphic has been tested in a BubR1-deficient system.** This refines the drug choice within the node; it does not upgrade the node's evidence class or its score.
5. **Human and paediatric safety.**
   - Intermittent senolytic dosing has been given to adults in early-phase trials (Hickson 2019, PMID 31542391). A trial in adult survivors of childhood cancer is ongoing (NCT04733534).
   - Dasatinib is approved in paediatric haematological indications. Its known toxicities include myelosuppression and pleural effusion.
   - There are no paediatric senolytic data.
   - **Rapalogs:** temsirolimus was added to chemotherapy in a Children's Oncology Group phase 3 trial in intermediate-risk rhabdomyosarcoma (ARST1431, 297 evaluable patients). It did **not** improve 3-year event-free survival, and it increased grade 3–4 anaemia (Gupta 2024, PMID 38936378). That is a negative efficacy result in RMS, not a safety failure, and it means a paediatric RMS safety dataset exists for this class. Everolimus has oral paediatric approval (SEGA in tuberous sclerosis). At high doses, rapalogs are immunosuppressive.
6. **Proposed experiment.**
   - **In cells:** quantify senescence burden (p16/p21, SA-β-gal, SASP cytokines) in *BUB1B*-mutant versus control primary cells. Then compare two arms:
     - intermittent dasatinib + quercetin (senolytic): does it selectively deplete senescent cells without affecting proliferating controls?
     - low-dose sirolimus (senomorphic): does it lower IL-6/IL-1α and the other SASP markers **while the proliferative arrest is retained rather than released**?
   - Run a conditioned-medium assay on rhabdomyosarcoma lines to test whether each arm reduces or increases the paracrine growth support from senescent cells.
   - **In BubR1-hypomorphic mice:** compare pharmacological senolysis against the published genetic-clearance readouts, with **tumour incidence as a co-primary endpoint**.

### 5.3 Evidence chain at a glance

Each hypothesis is only as strong as its weakest link. The table grades every link separately, by the system the evidence comes from, so the gaps are visible rather than averaged away.

| Hypothesis | Link in the chain | Evidence (system) | Strength |
|---|---|---|---|
| both | L751\* removes one allele's protein | Truncating *BUB1B* alleles reduced by NMD in PCS/MVA patient cells (human cells; PMID 16411201) | strong |
| both | N1016K impairs the remaining protein | Pseudokinase-domain position (PMID 33207204); discordant in-silico predictors; no assay | **weak** |
| both | Less BubR1 causes aneuploidy and the progeroid / muscle phenotype | Hypomorphic mice (PMID 15208629); >50 % BubR1 decrease sufficient in patients (PMID 16411201) | strong |
| H1 | More BubR1 is protective | Overexpression mice (PMID 23242215) | strong (mouse) |
| H1 | NAD+ precursors raise BubR1 via SIRT2 | NMN in mice; SIRT2 overexpression in BubR1-hypomorphic mice (PMID 24825348) | moderate (one study, mouse) |
| H1 | The route works on the N1016K protein | K668 lies outside both variants (inference) | **absent — §7 experiment** |
| H1 | Tolerated in humans, including children | NR in ataxia-telangiectasia (PMIDs 34515380, 37899683); RCT in Werner syndrome (PMID 40459998) | moderate (other diseases) |
| H1 | Does not promote rhabdomyosarcoma | RMS cells depend on FoxM1/BUB1B (PMID 23002205) | **open risk — §7, §9.1** |
| H2 | Low BubR1 drives senescent-cell accumulation | BubR1-hypomorphic mice (PMID 22048312) | strong (mouse) |
| H2 | Removing that burden delays muscle wasting | Genetic p16+ clearance in the same mice (PMID 22048312) | moderate (genetic, not drug) |
| H2 | Rapalogs suppress the SASP and its tumour-promoting effect | Human cells and mouse xenografts (PMIDs 26147250, 26280535) | moderate |
| H2 | Paediatric exposure in RMS | ARST1431 phase 3 (PMID 38936378) | strong for safety, none for efficacy |
| H2 | Effective in a BubR1-deficient system | Not tested | **absent — §7 experiment** |

Both hypotheses fail or stand on the links marked **absent**, and both are tested by the same experiment (§7). That is deliberate: the cheapest decisive experiment targets the weakest link.

---

## 6. Other routes, reported for completeness

### 6.1 H3 — Spindle-poison interaction (a research question, not a candidate)

Vinca alkaloids and taxanes kill cells by triggering a mitotic arrest that **depends on a working SAC**. Cells with a weakened checkpoint are reported to slip through mitosis instead. Both drug classes are standard in rhabdomyosarcoma protocols.

A published MVA case with rhabdomyosarcoma was managed with reduced-intensity chemotherapy (Nishitani-Isa 2019, PMID 31184400). Genomic-instability disorders as a group show altered treatment sensitivity (Nakano 2024, PMID 39264246).

**How a constitutionally weakened SAC shifts spindle-poison dose–response is a directly testable question:** live-cell imaging of mitotic arrest duration, slippage and viability in *BUB1B*-mutant versus control cells. **We make no clinical claim here; treatment decisions belong to the oncology team.**

### 6.2 H4 — APC/C–CDC20 braking (a documented druggability gap)

This is the most direct mechanistic compensation for a failed mitotic checkpoint complex. Preclinical work with the APC/C antagonist proTAME after SAC disruption is limited to mouse oocytes (Yonekura 2025, PMID 40887771). But no approved drug acts on CDC20, ANAPC2 or UBE2C in any direction; only tool compounds exist (proTAME, apcin). We report the gap rather than drop the node silently.

### 6.3 H5 — Targeting the aneuploid state (selective toxicity, not rescue)

This route has the best drug availability in the analysis:
- proteasome inhibitors;
- DHODH inhibitors, including paediatric safety data for teriflunomide (TERIKIDS, PMID 34800398).

Aneuploid cells carry proteotoxic and metabolic burdens that make them selectively vulnerable (Tang 2011, PMID 21315436). But in a body where aneuploidy is constitutional, killing aneuploid cells corrects nothing. H5 is the worked example of why tractability cannot carry a hypothesis: it ranks third.

---

## 7. Shared validation design

Both prioritised hypotheses can be tested in one cellular system:

| Arm | Cells | Treatment | Primary readouts | What it distinguishes |
|---|---|---|---|---|
| Rescue | *BUB1B*-mutant (patient-derived / isogenic) vs control | vehicle / NR / D+Q / low-dose sirolimus | micronuclei, mis-segregation, aneuploidy, SAC arrest duration, BubR1 level | phenotypic rescue |
| Selectivity | same | same | viability and proliferation of non-senescent cells | rescue vs selective toxicity |
| Tumour safety | rhabdomyosarcoma lines (± BUB1B knock-down) | same | proliferation, colony formation | rescue vs pro-tumour effect |

**A positive result in the rescue arm alone is not sufficient.** The tumour-safety arm is required before any translational claim.

---

## 8. Limitations

- **The N1016K effect is unproven.** Every hypothesis assumes it produces a partially functional protein.
- **Phase is unknown.** Without parental DNA we cannot confirm the variants are in trans.
- **Most mechanistic evidence comes from BubR1-hypomorphic mice**, not from human MVA.
- **Network distance (STRING) is supporting evidence only**, and depends on the confidence threshold.
- **Drug–target mapping depends on the curation status of ChEMBL and DGIdb** at query time. Results were re-run on 2026-10-13; UniProt accessions are cached for reproducibility.
- **CMap/LINCS transcriptomic analysis**, proposed in our procedure as supporting evidence, **was not performed**. No hypothesis depends on it.
- **Scores are hand-weighted**, not learned. The weights are explicit so they can be challenged.

---

## 9. Innovation and scalability

- **Direction-aware drug mapping.** Every drug–target link records whether the drug pushes the target in the *desired* direction. This is how PLK1 and CENPE inhibitors, the most prominent network hits, were correctly rejected.
- **Discards are categorised, not just filtered.** The five reasons for setting a route aside (no drug, wrong direction, phenotype mismatch, published negative result, toxicity rather than rescue) are themselves an output that other teams and clinicians can reuse.
- **Phenotype-anchored filtering.** Mechanisms were checked against *this child's* documented phenotype. As a result, two well-evidenced MVA mechanisms (microcephaly, heart defects) were discarded for this child.
- **Scalability.** The pipeline is driven by a gene-level target panel (a TSV) and public APIs. Re-running it for another MVA gene (*CEP57*, *TRIP13*) or another monogenic disorder only requires a new panel. No patient data is used anywhere in Track 2.

**Scalability, demonstrated rather than claimed.** We re-ran the drug lookup automatically for all eight genes with a reported MVA or MVA-like phenotype (`notebooks/utils/track2_step6_mva_druggability_scan.py`). For each gene, it takes the top 15 high-confidence STRING interactors (score ≥ 0.7), and queries ChEMBL for approved drugs acting on any target that contains each protein, in any direction. The run takes about four minutes, needs no curation, and uses no patient data.

| MVA gene | Gene itself druggable | First-shell interactors | With any approved drug |
|---|---|---|---|
| *BUB1B* (MVA1) | no | 15 | 0 |
| *CEP57* (MVA2) | no | 5 | 0 |
| *TRIP13* (MVA3) | no | 15 | 0 |
| *MAD2L1BP* | no | 13 | 0 |
| *CEP192* | no | 15 | 1 (EGLN2: HIF prolyl-hydroxylase inhibitors) |
| *SMC5* | no | 15 | 1 (DDB1: lenalidomide, pomalidomide via CRBN–DDB1) |
| *BUB1* | no | 15 | 0 |
| *CENATAC* | no | 15 | 0 |

**2 of 108** interactors (90 unique proteins) have any approved drug, and neither acts on chromosome segregation: both hits come from complexes that the partner protein joins for unrelated functions. Positive controls behave as expected (PSMB5 → proteasome inhibitors, ABL1 → tyrosine-kinase inhibitors, TUBB → taxanes and vinca alkaloids). The finding of §3.2 therefore generalises: **across MVA, the defect sits in an undrugged part of the proteome**, and repurposing for any MVA child will have to act on consequences — dosage, senescence, treatment sensitivity — as we did here. This is a disease-level result, not only a case-level one.

Because the scan is direction-agnostic it is an upper bound on druggability; choosing the direction of modulation remains the expert step, and the curated-panel method of Step 3 is how it is done. Outputs: `output/track2/track2_step6_mva_scan_summary.csv` and `track2_step6_mva_scan_proteins.csv`.

### 9.1 A fast route to the H1 tumour-safety question

**Therapeutic-window pre-screen, then a single 48-hour co-culture readout.** Of all the open problems in this report, the H1 contradiction — whether raising NAD+ helps the healthy tissue *and* whether it feeds residual rhabdomyosarcoma — is the one where a proposed method could compress months of work into days without leaving the set of tools the field already trusts. Both steps rest on resources that are standard practice rather than novel: the NIH **LINCS L1000 / Connectivity Map** (the reference catalogue of transcriptional signatures for approved drugs and their doses, including niacin, nicotinamide riboside and sirolimus) and the Broad **DepMap**, which carries CRISPR and drug-sensitivity measurements for real human rhabdomyosarcoma lines. Fine-tuning a standard transfer-learning drug-response model on the DepMap rhabdomyosarcoma profiles answers the H1 question computationally before any reagent is ordered: at each candidate dose, does the signature turn on the SIRT2/NAD+ rescue arm *without* turning on the proliferation programme this tumour depends on, above all the **FoxM1/BUB1B axis** (Wan 2012, PMID 23002205)? Only the single dose-combination that survives that filter goes to the bench, where a live-cell imaging plate (Incucyte-type phase-contrast) analysed by **Cellpose**, the de-facto open-source segmentation model in cell biology, measures both arms simultaneously — micronuclei and mis-segregation in the *BUB1B*-mutant cells, and division rate in fluorescently tagged rhabdomyosarcoma cells in the same well. This collapses the shared validation design in §7 into one plate and turns "does H1 help or harm?" into a 48-hour camera run that any university lab can perform, and it produces the therapeutic window as an explicit deliverable: the dose range in which rescue is separable from tumour promotion.

> **Status:** this is a proposal, not a result. No LINCS, DepMap or imaging analysis was run for this submission (§8), and the method reduces the cost of the tumour-safety arm rather than replacing it. Its value to this case is that the H1 precondition in §5 — safe only in oncological remission, and only after a laboratory go/no-go decision — currently depends on an experiment nobody has scheduled; this route makes that decision a two-day quantity instead of a two-month one. It is also the cheapest way to have a laboratory answer ready for the *next* MVA child, whose phenotype may differ from this one.

---

## 10. Potential impact on understanding and diagnosing MVA

If the experiments proposed here are run, their value does not depend on either drug working.

- **A functional classification for the child's missense variant.** p.Asn1016Lys has no ClinVar record, and in-silico predictors disagree (§2.2). The rescue arm of §7 measures BubR1 level, checkpoint function and segregation errors in cells carrying it. That is functional evidence in the ACMG/AMP framework (criteria PS3/BS3; Richards 2015, PMID 25741868). Submitting it to ClinVar would turn an unclassified variant into a classified one for the next family that carries it. Confirming that the two variants are on different copies of the gene still needs parental samples or long-read sequencing.
- **A reusable cellular readout for new *BUB1B* variants.** The same assay panel (BubR1 protein, micronuclei, aneuploidy by single-cell karyotyping, SAC arrest duration, p16/SASP senescence burden) is a quantitative phenotype for any patient-derived or isogenic *BUB1B* line. It could serve as a functional test for future variants of uncertain significance and as an endpoint for any MVA intervention study.
- **Chemotherapy sensitivity in MVA children with cancer (H3).** Whether a constitutionally weakened checkpoint shifts the response to vinca alkaloids and taxanes is directly relevant to how children with MVA and rhabdomyosarcoma or Wilms tumour are dosed. Today the only guidance is case reports (PMID 31184400).
- **A disease-level map of where intervention is possible.** The MVA-wide scan (§9) shows that the proteins closest to every known MVA gene are essentially undrugged. Future work on any MVA subtype can start from the consequence layer (BubR1 dosage, senescence, treatment sensitivity) instead of rediscovering the gap.
- **A testable senescence component of MVA.** If senescence burden is elevated in patient cells and responds to senomorphics, part of the growth failure and muscle wasting seen in MVA would be attributable to a modifiable process, not only to aneuploidy itself.

---

## 11. Reproducibility

```bash
python notebooks/utils/track2_step3_drug_mapping.py   # ~5 min; live ChEMBL / DGIdb / STRING queries
python notebooks/utils/track2_step4_prioritize.py     # offline, instant
python notebooks/utils/track2_step6_mva_druggability_scan.py  # ~4 min; same lookup for all 8 MVA genes
python docs/figures/build_track2_mechanisms_figure.py # explanatory figure (EN/ES)
```

**Outputs** (`output/track2/`):
- `track2_step3_drug_target_table.csv`
- `track2_step3_target_summary.csv`
- `track2_step3_string_network.csv`
- `track2_step4_prioritised_hypotheses.csv`
- `track2_step6_mva_scan_summary.csv`, `track2_step6_mva_scan_proteins.csv`

**Full analysis trail:** `docs/01-procedure.md`, `docs/02-findings.md`, `docs/00-references.md`.

**Public data sources:**
- ChEMBL (EBI REST API)
- DGIdb (GraphQL)
- STRING (REST)
- UniProt (REST)
- Europe PMC (citation verification)
- PubMed
- ClinicalTrials.gov

Open Targets was evaluated and not used: its `knownDrugs` field was unavailable at query time. No proprietary data were used.

---

## 12. AI usage disclosure

> **TODO (team — required):** provider, plan/tier and data-handling setting for each AI tool used, e.g. "Anthropic, Claude Opus via Claude Code, <plan>, <training/retention setting>" and "Google, Gemini Pro, <plan>, <setting>".

AI assistants were used for literature triage, code writing and drafting text. Every PMID cited in this report was checked against Europe PMC. The scientific choices — the target panel, the directions of modulation, the scoring weights and the final hypotheses — were reviewed by the team.

---

## References

1. Baker DJ et al. 2004. BubR1 insufficiency causes early onset of aging-associated phenotypes and infertility in mice. *Nat Genet.* PMID 15208629
2. Baker DJ et al. 2011. Clearance of p16Ink4a-positive senescent cells delays ageing-associated disorders. *Nature.* PMID 22048312
3. Baker DJ et al. 2013. Increased expression of BubR1 protects against aneuploidy and cancer and extends healthy lifespan. *Nat Cell Biol.* PMID 23242215
4. Gama Braga L et al. 2020. BUBR1 pseudokinase domain promotes kinetochore PP2A-B56 recruitment, spindle checkpoint silencing, and chromosome alignment. *Cell Rep.* PMID 33207204
5. Hanks S et al. 2004. Constitutional aneuploidy and cancer predisposition caused by biallelic mutations in BUB1B. *Nat Genet.* PMID 15475955
6. Hickson LJ et al. 2019. Senolytics decrease senescent cells in humans: preliminary report from a clinical trial of dasatinib plus quercetin in individuals with diabetic kidney disease. *EBioMedicine.* PMID 31542391
7. Irie J et al. 2020. Effect of oral administration of nicotinamide mononucleotide on clinical parameters and nicotinamide metabolite levels in healthy Japanese men. *Endocr J.* PMID 31685720
8. Matsuura S et al. 2006. Monoallelic BUB1B mutations and defective mitotic-spindle checkpoint in seven families with premature chromatid separation (PCS) syndrome. *Am J Med Genet A.* PMID 16411201
9. Nakano Y et al. 2024. Update on recommendations for cancer screening and surveillance in children with genomic instability disorders. *Clin Cancer Res.* PMID 39264246
10. Nishitani-Isa M et al. 2019. Rhabdomyosarcoma with premature chromatid separation-mosaic variegated aneuploidy syndrome: reduced-intensity chemotherapy. *Pediatr Int.* PMID 31184400
11. North BJ et al. 2014. SIRT2 induces the checkpoint kinase BubR1 to increase lifespan. *EMBO J.* PMID 24825348
12. Pirinen E et al. 2020. Niacin cures systemic NAD+ deficiency and improves muscle performance in adult-onset mitochondrial myopathy. *Cell Metab.* PMID 32386566
13. Tang YC et al. 2011. Identification of aneuploidy-selective antiproliferation compounds. *Cell.* PMID 21315436
14. Veenhuis SJG et al. 2021. Nicotinamide riboside improves ataxia scores and immunoglobulin levels in ataxia telangiectasia. *Mov Disord.* PMID 34515380
15. Wan X et al. 2012. Identification of FoxM1/Bub1b signaling pathway as a required component for growth and survival of rhabdomyosarcoma. *Cancer Res.* PMID 23002205
16. Chitnis T et al. 2021. Safety and efficacy of teriflunomide in paediatric multiple sclerosis (TERIKIDS). *Lancet Neurol.* PMID 34800398
17. Sterling NA et al. 2023. P53 independent pathogenic mechanisms contribute to BubR1 microcephaly. *Front Cell Dev Biol.* PMID 37900274
18. González-Blanco A et al. 2026. Proteostasis failure and mitochondrial dysfunction contribute to chromosomal instability-induced microcephaly. *Nat Commun.* PMID 41820377
19. Pun R et al. 2025. BubR1 controls heart development by promoting expression of cardiogenesis regulators. *J Am Heart Assoc.* PMID 40055864
20. Tabe C et al. 2026. Spindle assembly checkpoint competency determines sensitivity to KIF18A inhibition in small cell lung cancer. *Cancer Res Commun.* PMID 42647139
21. Feng B et al. 2022. A novel CEP57 variant associated with mosaic variegated aneuploidy syndrome in a Chinese female presenting with short stature. *Mol Genet Genomic Med.* PMID 35434947
22. Yonekura S et al. 2025. Effects of reversine and proTAME treatment on chromosome segregation during mouse oocyte maturation. *Zygote.* PMID 40887771
23. SEN-SURVIVORS trial (senolytics in adult survivors of childhood cancer). ClinicalTrials.gov NCT04733534
24. Navas LE et al. 2023. NAD pool as an antitumor target against cancer stem cells in head and neck cancer. *J Exp Clin Cancer Res.* PMID 36864434
25. Presterud R et al. 2024. Long-term nicotinamide riboside use improves coordination and eye movements in ataxia telangiectasia. *Mov Disord.* PMID 37899683
26. Shoji M et al. 2025. Nicotinamide riboside supplementation benefits in patients with Werner syndrome: a double-blind randomized crossover placebo-controlled trial. *Aging Cell.* PMID 40459998
27. Di Micco R et al. 2021. Cellular senescence in ageing: from mechanisms to therapeutic opportunities. *Nat Rev Mol Cell Biol.* PMID 33328614
28. Colucci M et al. 2025. Senescence in cancer. *Cancer Cell.* PMID 40513577
29. Zhang L et al. 2023. Targeting cellular senescence with senotherapeutics: senolytics and senomorphics. *FEBS J.* PMID 35015337
30. Laberge RM et al. 2015. MTOR regulates the pro-tumorigenic senescence-associated secretory phenotype by promoting IL1A translation. *Nat Cell Biol.* PMID 26147250
31. Herranz N et al. 2015. mTOR regulates MAPKAPK2 translation to control the senescence-associated secretory phenotype. *Nat Cell Biol.* PMID 26280535
32. Gupta AA et al. 2024. Addition of temsirolimus to chemotherapy in children, adolescents, and young adults with intermediate-risk rhabdomyosarcoma (ARST1431). *Lancet Oncol.* PMID 38936378
33. Moiseeva O et al. 2013. Metformin inhibits the senescence-associated secretory phenotype by interfering with IKK/NF-kB activation. *Aging Cell.* PMID 23521863
34. Xu M et al. 2015. JAK inhibition alleviates the cellular senescence-associated secretory phenotype and frailty in old age. *Proc Natl Acad Sci U S A.* PMID 26578790
35. Richards S et al. 2015. Standards and guidelines for the interpretation of sequence variants: a joint consensus recommendation of the ACMG and AMP. *Genet Med.* PMID 25741868
