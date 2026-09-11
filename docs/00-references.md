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
