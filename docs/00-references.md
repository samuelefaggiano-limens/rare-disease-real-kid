# Literature Research & References

This document maintains foundational literature on Mosaic Variegated Aneuploidy (MVA) and tracks specific research useful for both Variant Prediction (Track 1) and Drug Repurposing (Track 2), following the core directives of the hackathon.

---

# MVA Track 1 Reference: What We Know & How to Use It


## 1. What is Mosaic Variegated Aneuploidy (MVA)?

Normally, when a cell divides, it copies its chromosomes exactly and splits them evenly between the two new cells. In MVA, this process breaks down: some cells end up with the wrong number of chromosomes. It's called:
- **"Mosaic"** because only *some* cells in the body are affected, not all of them.
- **"Variegated"** because different cells have different chromosomes affected — it's not always the same one.

This condition is caused by broken "quality control" machinery that normally checks chromosomes are lined up and divided correctly. People with MVA often have developmental problems (slow growth, brain differences) and a higher risk of cancer, because faulty chromosome control can let abnormal cells grow unchecked.

Source: [Mosaic variegated aneuploidy in development, ageing and cancer (review)](https://pubmed.ncbi.nlm.nih.gov/39169218/)


## 2. The Known MVA Gene Panel

MVA is not caused by a single gene — several different genes, when broken, can each independently cause an MVA-like picture. This is the practical value for **Track 1**: instead of scanning the entire genome for the causal variant, you can start by checking whether the child's DNA has damaging variants in **these known genes first**, since they're already proven to cause MVA-type disease in other patients.

| Gene | MVA Subtype | Inheritance pattern seen | Typical features reported | Source |
|---|---|---|---|---|
| **BUB1B** | MVA syndrome 1 | Compound heterozygous (one broken copy from each parent, different mutations) | Structural brain anomalies, developmental delay (seen prenatally in 2nd trimester) | [Genetic analysis of two fetuses with MVA (BUB1B)](https://pubmed.ncbi.nlm.nih.gov/40555658/) |
| **CEP57** | MVA syndrome 2 | Homozygous (same broken copy from both parents; more common when parents are related) | Intrauterine growth restriction, severe short stature, unusual facial features, skeletal abnormalities, congenital heart defects, delayed development | [Novel CEP57 mutation in MVA syndrome 2](https://pubmed.ncbi.nlm.nih.gov/41700350/) |
| **TRIP13** | MVA syndrome 3 | Null (loss-of-function) variant | Strong link to **multiple tumors over a lifetime** (e.g., Wilms' tumor as a young child, followed years later by a tumor behind the eye) | [Wilms' tumor and orbital rhabdomyosarcoma in MVA syndrome 3](https://pubmed.ncbi.nlm.nih.gov/42595739/) |
| **MAD2L1BP** (makes the protein p31comet) | Related aneuploidy/tumor-predisposition condition | Homozygous nonsense variant (biallelic) | Small head size (microcephaly), severe epilepsy, ovarian/testicular tumors in childhood; patient cells showed severe MVA and slowed cell growth in the lab | [Biallelic MAD2L1BP mutation and juvenile granulosa cell tumors](https://pubmed.ncbi.nlm.nih.gov/37796616/) |
| **CEP192** | MVA syndrome with tetraploidy | Biallelic (recessive) | Microcephaly, short stature, limb-extremity dysplasia, tetraploidy, and predisposition to male infertility. | [Mosaic variegated aneuploidy syndrome with tetraploid... mutant CEP192](https://pubmed.ncbi.nlm.nih.gov/37981762/) |
| **SMC5** | Atelis Syndrome / MVA with near-tetraploidy | Compound heterozygous / recessive | Near-tetraploidy, mild developmental delay, hypospadias (milder presentation) or severe growth failure and multi-system abnormalities. | [Recessive SMC5 Variants in a Family with Near-Tetraploidy/Mosaic Variegated Aneuploidy](https://pubmed.ncbi.nlm.nih.gov/41374403/) |
| **BUB1** | Primary Microcephaly 30 (MCPH30) / MVA syndrome | Compound heterozygous | Microcephaly, severe neurodevelopmental disorders, chromosome segregation defects. | [Generation and characterization of... cell lines from a young female with microcephaly carrying a compound heterozygous mutation in BUB1 gene](https://pubmed.ncbi.nlm.nih.gov/39490209/) |

**How to use this table for Track 1:** Match the child's clinical features (growth, brain findings, tumor history, seizures, etc.) against the "Typical features" column. Whichever gene's known phenotype overlaps best with the child's actual symptoms is the strongest candidate to search first in his genome data — this is the essence of *phenotype-guided* prioritization rather than blindly ranking all variants genome-wide.

## 3. Rationale for Prioritizing Coding and Splicing Variants in WGS Analysis
**High Diagnostic Yield in Coding Regions:** Although protein-coding exomes represent merely ~1–2% of the human genome, approximately 85% of documented pathogenic variants responsible for monogenic, Mendelian disorders are located within exons and canonical splice junctions.

**Clear Functional Predictability:** Coding variants—such as frameshifts, nonsense (premature stop codons), canonical splice site disruptions (+/- 1 bp, +/-2 bp), and severe missense changes have interpretable, direct consequences on transcript architecture and protein stability under established clinical guidelines.

**Signal-to-Noise Ratio & Computational Limitations:** Whole Genome Sequencing (WGS) yields millions of non-coding variants per individual. Current in silico tools have limited specificity in deep intronic and intergenic spaces, leading to high false-positive rates and uninterpretable Variants of Uncertain Significance (VUS) without transcriptomic validation.

**Filtering Strategy:** Prioritizing exonic and flanking splice regions is a standard Tier-1 clinical genomics workflow. 

***Non-coding intronic sequences are routinely excluded in early filtration steps unless:*** No causal candidate emerges in primary coding analysis.

***The variant occurs within canonical/extended splice-donor/acceptor regions.***

***A well-characterized deep-intronic cryptic splice site or regulatory variant matches the phenotype and can be validated (e.g., via RNA-Seq or functional reporter assays).***

## 4. A Critical Filtering Caution (Don't Over-Filter by Population Frequency)

When processing a whole-genome VCF file, a standard step is to remove "common" variants — ones seen too frequently in the general population — on the assumption that a variant causing a rare disease should itself be rare. A typical cutoff is **1% frequency** in reference population databases.

However, one study found a disease-causing variant (in a different gene, *MOCS2*, unrelated to MVA) that reaches **~3.6% frequency** within one specific subpopulation (the Roma population), even though it is very rare worldwide. This happens because of a **founder effect** — a mutation that arose in one ancestral group and became relatively common within that group through inheritance, while staying rare everywhere else.

**Practical implication for Track 1:** If your pipeline applies a strict, single global frequency cutoff (e.g., "remove anything above 1%"), you risk accidentally throwing away the real causal variant — especially if the family belongs to a smaller or historically isolated population group. A safer approach is to use a lower and/or population-specific frequency threshold rather than one blanket global cutoff.

Source: [A prevalent MOCS2 variant in the Roma population and mild molybdenum cofactor deficiency](https://pubmed.ncbi.nlm.nih.gov/40707723/)


## Sources at a Glance

1. Review — MVA overview: https://pubmed.ncbi.nlm.nih.gov/39169218/
2. BUB1B compound heterozygous variants in two fetuses: https://pubmed.ncbi.nlm.nih.gov/40555658/
3. Novel CEP57 mutation (MVA syndrome 2): https://pubmed.ncbi.nlm.nih.gov/41700350/
4. TRIP13 variant and multiple tumors (MVA syndrome 3): https://pubmed.ncbi.nlm.nih.gov/42595739/
5. MAD2L1BP / p31comet mutation and tumors: https://pubmed.ncbi.nlm.nih.gov/37796616/
6. MOCS2 founder variant and filtering lesson: https://pubmed.ncbi.nlm.nih.gov/40707723/
7. Genetic diagnosis by whole exome capture: https://pubmed.ncbi.nlm.nih.gov/19861545/
8. ACMG/AMP sequence variant interpretation standards: https://doi.org/10.1038/gim.2015.30
9. The Human Gene Mutation Database (HGMD): https://pubmed.ncbi.nlm.nih.gov/28349240/
10. Clinical interpretation of non-coding genome variants: https://pmc.ncbi.nlm.nih.gov/articles/PMC12192682/
11. CEP192 and MVA with tetraploidy: https://pubmed.ncbi.nlm.nih.gov/37981762/
12. SMC5 variants and Atelis / MVA: https://pubmed.ncbi.nlm.nih.gov/41374403/
13. BUB1 and MCPH30 / MVA: https://pubmed.ncbi.nlm.nih.gov/39490209/
---


# MVA Track 2 Reference: Disease Mechanisms & Potential Drug Targets


## 1. Why This Matters for Track 2

Track 1 is about finding *which gene is broken*. Track 2 is about a different question: **once we know a gene is broken, what actually goes wrong inside the cells, and can an existing (already-approved) drug counteract that damage?**

The three sources below don't identify new causal genes — instead, they dig into the downstream biology of already-known MVA genes (mainly **BUB1B**, also known as **BubR1**) and point to specific cellular processes that could realistically be targeted with medicines that already exist. That's exactly what "drug repurposing" needs: a mechanism plus a drug class that already touches it.


## 2. The Three Mechanisms, Explained Simply

### A. Brain stem cells run out of energy and can't keep up (proteostasis + mitochondria)

Cells with too many/too few chromosomes (aneuploid cells) have to work much harder just to stay alive — they have too much or too little protein to manage, and this constant stress overloads their "protein quality control" system (called proteostasis) and their energy-producing mitochondria. This is especially damaging to **neural stem cells (NSCs)**, the cells that build the brain, because staying a stem cell and dividing repeatedly is an energy-expensive job. When these stressed stem cells can't keep up, fewer of them survive and divide, which is a likely reason MVA causes **microcephaly** (a smaller brain/head).

**What could help (druggable angle):** The study found that boosting three things reduced the damage in lab models:
- **ROS scavengers** — substances that mop up harmful reactive oxygen species ("cell rust") produced by struggling mitochondria.
- **Mitochondrial chaperones** — helper proteins/drugs that assist mitochondria in folding proteins correctly under stress.
- **Apoptosis inhibitors** — drugs that block the "self-destruct" signal, giving stressed-but-still-viable stem cells a chance to survive instead of dying.

All three are existing drug classes (antioxidants, mitochondria-protective compounds, anti-apoptotic drugs), making this a strong candidate area to search for repurposable medicines.

Source: [Proteostasis failure and mitochondrial dysfunction contribute to chromosomal instability-induced microcephaly](https://pubmed.ncbi.nlm.nih.gov/41820377/)


### B. Blocking one death pathway (P53) is not enough on its own

When the BUB1B gene (which makes the protein **BubR1**) is broken, brain cells accumulate chromosome damage. Normally, damaged cells are cleared out through a well-known "cell suicide" pathway controlled by the protein **P53**. It would be tempting to think: *just block P53, and the cells will survive.*

This study tested that idea directly in mice (removing both BubR1 and P53 at the same time) and found it **doesn't fully work** — a lot of brain cell death still happens even without P53, meaning there's a second, separate, still-not-fully-understood cell-death pathway also killing off cells.

**Why this matters for Track 2:** It's an important negative result. A drug strategy built purely around inhibiting P53 (or the classic apoptosis pathway it controls) will likely be insufficient on its own to protect the brain in MVA. Any proposed treatment approach should account for this second pathway, or be combined with something else (like the mitochondrial/ROS-targeting strategies from source A above) rather than relying on P53-blockade alone.

Source: [P53-independent pathogenic mechanisms contribute to BubR1 microcephaly](https://pubmed.ncbi.nlm.nih.gov/37900274/)


### C. Heart defects are linked to a hyperactive Wnt signaling pathway

Many children with MVA (caused by BUB1B/BubR1 mutations) are also born with heart defects, such as holes in the heart wall (septal defects). This study used single-cell analysis of BubR1-deficient mouse hearts and found that the cells destined to become heart muscle cells (cardiomyocytes) don't develop properly — key genes needed for normal heart formation are suppressed.

The likely cause: the **Wnt signaling pathway** (a well-known cell-communication system that tells cells when to grow, specialize, and organize into tissues) becomes **overactive** in the BubR1-deficient heart, disrupting the normal balance needed to build heart tissue correctly.

**Why this matters for Track 2:** Wnt signaling is one of the most heavily studied, most "druggable" pathways in biology — there are already many existing Wnt-pathway inhibitor compounds (some already tested in cancer trials, since Wnt overactivation is also common in tumors). This makes Wnt inhibition a promising, plausible repurposing angle specifically for the **heart defect** side of MVA, separate from the brain/microcephaly angle above.

Source: [BubR1 controls heart development by promoting expression of cardiogenesis regulators](https://pubmed.ncbi.nlm.nih.gov/40055864/)


## 3. Putting It Together: Three Angles for Track 2

| Affected system | What goes wrong | Druggable strategy pointed to | Source |
|---|---|---|---|
| Brain (neural stem cells) → microcephaly | Protein-handling and mitochondrial stress overload energy-hungry stem cells, causing them to die or stop dividing | ROS scavengers, mitochondrial chaperones, apoptosis inhibitors | [Link](https://pubmed.ncbi.nlm.nih.gov/41820377/) |
| Brain (cortex) → microcephaly | Chromosome-damaged cells die via P53 **and** a separate, unidentified P53-independent pathway | Don't rely on P53-pathway drugs alone; likely needs a combination approach | [Link](https://pubmed.ncbi.nlm.nih.gov/37900274/) |
| Heart → congenital defects (e.g., septal defects) | Wnt signaling becomes hyperactive, disrupting normal heart-cell development | Existing Wnt-pathway inhibitor drugs (already explored in other diseases, e.g., cancer) | [Link](https://pubmed.ncbi.nlm.nih.gov/40055864/) |

**Suggested framing for a Track 2 proposal:** Since the child's own specific symptoms (from Track 1's genetic work) will point to which organ systems are most affected, this table lets you match symptoms to mechanism and then to an existing drug class:
- If microcephaly/brain involvement is prominent → look at repurposing antioxidants, mitochondria-protective agents, or anti-apoptotic compounds (but consider combination therapy, given the P53-independent pathway finding).
- If congenital heart defects are prominent → look at repurposing existing Wnt-pathway-modulating drugs.


## Sources at a Glance

1. Proteostasis failure and mitochondrial dysfunction in chromosomal-instability microcephaly: https://pubmed.ncbi.nlm.nih.gov/41820377/
2. P53-independent pathogenic mechanisms in BubR1 microcephaly: https://pubmed.ncbi.nlm.nih.gov/37900274/
3. BubR1 controls heart development via cardiogenesis regulators (Wnt signaling): https://pubmed.ncbi.nlm.nih.gov/40055864/


## 4. BubR1 Dosage: the Evidence That a Repurposing Target Exists at All

The three sources above describe what goes wrong *downstream* of a broken BUB1B. A separate body of work asks a more useful question for Track 2: **does the amount of BubR1 protein itself matter, and can anything change it?** If BubR1 levels were irrelevant, or unreachable by any drug, the only remaining options would be damage control. They are not.

### A. Less BubR1 causes the disease phenotype; more BubR1 protects against it

Mice engineered to make only a small fraction of the normal BubR1 develop progressive aneuploidy together with short lifespan, dwarfism, muscle wasting and cataracts — a premature-ageing picture ([Baker et al. 2004, *Nat Genet*](https://pubmed.ncbi.nlm.nih.gov/15208629/)). The converse experiment is the more informative one: mice engineered to make *extra* BubR1 were protected against aneuploidy and cancer, kept better muscle fibres, performed better on treadmill tests, and lived longer ([Baker et al. 2013, *Nat Cell Biol*](https://pubmed.ncbi.nlm.nih.gov/23242215/)).

**Why this matters for Track 2:** it establishes the direction of travel. BubR1 behaves as a dosage-sensitive protein across its whole range, not as a switch that is simply on or off. Raising it is beneficial in the models where it has been tested. That makes "increase the BubR1 that remains" a coherent therapeutic goal rather than a wish.

### B. There is a known molecular handle on BubR1 levels, and it is reachable by a supplement

BubR1 abundance declines with age because NAD+ declines. The link is the enzyme SIRT2, which keeps BubR1 deacetylated at a specific residue (lysine 668) and thereby protects it from degradation; the acetyltransferase CBP pushes the other way. Treating mice with the NAD+ precursor NMN raised BubR1 *in vivo*, and overexpressing SIRT2 extended the lifespan of the BubR1-deficient mice from source A ([North et al. 2014, *EMBO J*](https://pubmed.ncbi.nlm.nih.gov/24825348/)).

**Why this matters for Track 2:** this is the only route we found that acts on the *quantity* of BubR1 rather than on its consequences. It also has an unusual property for a rare-disease target — the intervention is a nutrient precursor with existing human exposure, not a compound that would have to be developed from scratch.

**Caveat recorded deliberately:** BUB1B is also a dependency of rhabdomyosarcoma cells, which this child has had ([Wan et al. 2012, *Cancer Res*](https://pubmed.ncbi.nlm.nih.gov/23002205/)). Raising BubR1 systemically in a cancer-predisposition background is therefore a mechanistically unresolved question, not an obviously good idea, and is written up as such.

### C. Prior human exposure to NAD+ precursors, including in children

Nicotinamide riboside was given to 24 patients with ataxia-telangiectasia — another genome-instability syndrome with cancer predisposition — 17 of them under 18 years old. It was tolerated, with improved ataxia scores and immunoglobulin levels ([Veenhuis et al. 2021, *Mov Disord*](https://pubmed.ncbi.nlm.nih.gov/34515380/)). Separately, niacin corrected systemic NAD+ deficiency and improved muscle performance in adult-onset mitochondrial myopathy ([Pirinen et al. 2020, *Cell Metab*](https://pubmed.ncbi.nlm.nih.gov/32386566/)).

**Why this matters for Track 2:** it is the difference between a hypothesis that could be tested in people one day and one that could not. Note what it is *not*: neither study is in MVA, and neither tested the BubR1 mechanism.


## 5. The Senescence Arm: a Second, Independent Consequence of Low BubR1

Low BubR1 does not only cause missegregated chromosomes. The same mice accumulate p16-positive **senescent cells** — cells that have stopped dividing but stay alive and secrete inflammatory signals. When those senescent cells were genetically removed, the onset of muscle loss, fat loss and cataracts was delayed, and removing them late in life still slowed progression of disorders that had already started ([Baker et al. 2011, *Nature*](https://pubmed.ncbi.nlm.nih.gov/22048312/)).

This is an unusually direct piece of evidence: the experiment was done *in the BubR1 model itself*, and the phenotype it rescued (muscle wasting) is one this child presents.

Compounds that kill senescent cells selectively — "senolytics" — were then derived by mapping the pro-survival network those cells depend on ([Zhu et al. 2015, *Aging Cell*](https://pubmed.ncbi.nlm.nih.gov/25754370/)), and have since been given to humans in early-phase trials ([Hickson et al. 2019, *EBioMedicine*](https://pubmed.ncbi.nlm.nih.gov/31542391/)). A trial in adult survivors of childhood cancer is ongoing ([NCT04733534](https://clinicaltrials.gov/study/NCT04733534)).

**The gap to be honest about:** the rescue was achieved by *genetic* clearance in mice, not by a drug, and not in MVA. And senescence is also a tumour-suppressive program — removing senescent cells in a cancer-predisposition background cuts both ways.


## 6. Cancer Predisposition and Treatment Sensitivity in MVA

*BUB1B* was the gene originally identified in MVA families, including two with embryonal rhabdomyosarcoma ([Hanks et al. 2004, *Nat Genet*](https://pubmed.ncbi.nlm.nih.gov/15475955/)) — the same tumour type that triggered investigation in this child.

Two further sources matter for how the mechanism interacts with cancer treatment rather than with repurposing:

- Genomic-instability disorders as a group, MVA included, are recognised to have altered sensitivity to radiotherapy and chemotherapy, which shapes surveillance and treatment choices ([Nakano et al. 2024, *Clin Cancer Res*](https://pubmed.ncbi.nlm.nih.gov/39264246/)).
- A published MVA case with rhabdomyosarcoma was managed with **reduced-intensity chemotherapy** ([Nishitani-Isa et al. 2019, *Pediatr Int*](https://pubmed.ncbi.nlm.nih.gov/31184400/)).

Mechanistically this is not a coincidence. Vinca alkaloids and taxanes kill by triggering a mitotic arrest that *depends on a working spindle assembly checkpoint*; where the checkpoint is attenuated, cells slip through mitosis instead of arresting, and resistance has been reported in cancer models. In BUB1B-related MVA the checkpoint is precisely what is broken.

**Scope note:** this is recorded as a research question the mechanism itself raises, not as a clinical recommendation. No comparative data exist in constitutional biallelic *BUB1B* deficiency, and treatment decisions belong to the treating team.


## 7. Routes Examined and Set Aside

Recorded so the reasoning is not repeated. See `docs/02-findings.md` (Track 2, Step 3) for the full map.

| Route | Why it was set aside | Source |
|---|---|---|
| Blocking P53 to stop cell death | Already tested: deleting P53 together with BubR1 did **not** rescue microcephaly. A second, P53-independent death pathway keeps operating. A published negative result. | [PMID 37900274](https://pubmed.ncbi.nlm.nih.gov/37900274/) |
| Growth hormone for short stature | Reported ineffective in a CEP57-related MVA case despite high IGF-1. The defect is mitotic, not endocrine. | [PMID 35434947](https://pubmed.ncbi.nlm.nih.gov/35434947/) |
| KIF18A inhibition (CIN-selective) | Requires a *competent* checkpoint to kill; SAC-deficient cells survive it. MPS1 inhibition can restore the dependency, but neither agent is approved. | [PMID 42647139](https://pubmed.ncbi.nlm.nih.gov/42647139/) |
| Targeting the aneuploid state (proteasome, metabolic) | Approved drugs exist, but this eliminates aneuploid cells rather than correcting the defect — selective toxicity, not rescue, in a body where aneuploidy is constitutional. | [PMID 21315436](https://pubmed.ncbi.nlm.nih.gov/21315436/) |
| Wnt inhibition; neural stem cell protection | Sound evidence, but for congenital heart defects and microcephaly respectively — **neither is reported for this child** in the clinical document. | §Track 2 above |


## Sources at a Glance (Track 2, Steps 3-4)

4. BubR1 insufficiency causes early ageing phenotypes (Baker 2004): https://pubmed.ncbi.nlm.nih.gov/15208629/
5. Increased BubR1 protects against aneuploidy and cancer (Baker 2013): https://pubmed.ncbi.nlm.nih.gov/23242215/
6. SIRT2 induces BubR1 to increase lifespan (North 2014): https://pubmed.ncbi.nlm.nih.gov/24825348/
7. Nicotinamide riboside in ataxia-telangiectasia, incl. children (Veenhuis 2021): https://pubmed.ncbi.nlm.nih.gov/34515380/
8. Niacin, NAD+ deficiency and muscle performance (Pirinen 2020): https://pubmed.ncbi.nlm.nih.gov/32386566/
9. Clearance of p16-positive senescent cells in the BubR1 model (Baker 2011): https://pubmed.ncbi.nlm.nih.gov/22048312/
10. Senolytic drugs from the senescent-cell survival network (Zhu 2015): https://pubmed.ncbi.nlm.nih.gov/25754370/
11. First-in-human senolytic trial (Hickson 2019): https://pubmed.ncbi.nlm.nih.gov/31542391/
12. SEN-SURVIVORS trial in childhood-cancer survivors: https://clinicaltrials.gov/study/NCT04733534
13. Biallelic BUB1B mutations, MVA and cancer predisposition (Hanks 2004): https://pubmed.ncbi.nlm.nih.gov/15475955/
14. FoxM1/Bub1b required for rhabdomyosarcoma growth (Wan 2012): https://pubmed.ncbi.nlm.nih.gov/23002205/
15. Reduced-intensity chemotherapy in MVA rhabdomyosarcoma (Nishitani-Isa 2019): https://pubmed.ncbi.nlm.nih.gov/31184400/
16. Cancer surveillance in genomic instability disorders (Nakano 2024): https://pubmed.ncbi.nlm.nih.gov/39264246/
17. Aneuploidy-selective antiproliferation compounds (Tang 2011): https://pubmed.ncbi.nlm.nih.gov/21315436/
18. SAC competency determines KIF18A sensitivity (2026): https://pubmed.ncbi.nlm.nih.gov/42647139/
19. CEP57 MVA case, growth hormone ineffective (Feng 2022): https://pubmed.ncbi.nlm.nih.gov/35434947/
20. BubR1 pseudokinase domain recruits PP2A-B56 (2020): https://pubmed.ncbi.nlm.nih.gov/33207204/

### Databases queried (Step 3)

| Resource | Role | Endpoint |
|---|---|---|
| **ChEMBL** | Primary source of curated drug mechanisms of action and approval phase | `https://www.ebi.ac.uk/chembl/api/data` |
| **DGIdb** | Cross-check of approved drug-gene interactions and their direction | `https://dgidb.org/api/graphql` |
| **STRING** | Network distance from each target to BUB1B (supporting evidence only) | `https://string-db.org/api/json` |
| **UniProt** | Gene symbol to reviewed human accession | `https://rest.uniprot.org/uniprotkb/search` |
| **Europe PMC** | Citation verification for every PMID recorded above | `https://www.ebi.ac.uk/europepmc/webservices/rest` |

Open Targets was evaluated and **not used**: its GraphQL schema no longer exposes `knownDrugs` on `Target`, and the replacement field was unstable at the time of writing. ChEMBL covers the same curated mechanism data directly.


## Sources at a Glance (Track 2, follow-up investigations and report)

Cited in `02-findings.md` (Track 2 follow-up) and/or `output/track2/lostanos_track2_report.md`; every PMID verified against Europe PMC.

21. Monoallelic BUB1B mutations and PCS syndrome; truncating alleles reduced by NMD, >50 % BubR1 decrease sufficient (Matsuura 2006): https://pubmed.ncbi.nlm.nih.gov/16411201/
22. Oral NMN in healthy Japanese men (Irie 2020): https://pubmed.ncbi.nlm.nih.gov/31685720/
23. Teriflunomide in paediatric multiple sclerosis, TERIKIDS (Chitnis 2021): https://pubmed.ncbi.nlm.nih.gov/34800398/
24. Effects of reversine and proTAME on chromosome segregation in mouse oocytes (Yonekura 2025): https://pubmed.ncbi.nlm.nih.gov/40887771/
25. NAD pool as an antitumour target against cancer stem cells (Navas 2023): https://pubmed.ncbi.nlm.nih.gov/36864434/
26. Long-term nicotinamide riboside in ataxia-telangiectasia (Presterud 2024): https://pubmed.ncbi.nlm.nih.gov/37899683/
27. Nicotinamide riboside in Werner syndrome, randomised crossover trial (Shoji 2025): https://pubmed.ncbi.nlm.nih.gov/40459998/
28. Cellular senescence in ageing (Di Micco 2021): https://pubmed.ncbi.nlm.nih.gov/33328614/
29. Senescence in cancer (Colucci 2025): https://pubmed.ncbi.nlm.nih.gov/40513577/
30. Senolytics and senomorphics (Zhang 2023): https://pubmed.ncbi.nlm.nih.gov/35015337/
31. MTOR regulates the pro-tumorigenic SASP via IL1A translation (Laberge 2015): https://pubmed.ncbi.nlm.nih.gov/26147250/
32. mTOR regulates MAPKAPK2 translation to control the SASP (Herranz 2015): https://pubmed.ncbi.nlm.nih.gov/26280535/
33. Temsirolimus added to chemotherapy in intermediate-risk rhabdomyosarcoma, ARST1431 (Gupta 2024): https://pubmed.ncbi.nlm.nih.gov/38936378/
34. Metformin inhibits the SASP via IKK/NF-κB (Moiseeva 2013): https://pubmed.ncbi.nlm.nih.gov/23521863/
35. Metformin and the hallmarks of ageing (Kulkarni 2020): https://pubmed.ncbi.nlm.nih.gov/32333835/
36. JAK inhibition alleviates the SASP and frailty (Xu 2015): https://pubmed.ncbi.nlm.nih.gov/26578790/
37. ACMG/AMP standards for the interpretation of sequence variants (Richards 2015): https://pubmed.ncbi.nlm.nih.gov/25741868/
