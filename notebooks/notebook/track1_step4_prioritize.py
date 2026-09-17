"""Track 1, Paso 4: priorizacion de genes por fenotipo (HPO) y modelo de herencia.

Entrada: el VCF anotado del Paso 3A (o el del 3B, que es un superconjunto) y
los terminos HPO del paciente (Paso 1). Salida: un ranking de genes y una
tabla de variantes candidatas con todo lo que el Paso 5 necesita para
calcular el epcr.

Que hace, en orden:
  1. Lee los terminos HPO del paciente. Por defecto los saca de la tabla del
     .docx clinico (igual que el notebook del Paso 1); se pueden pasar con
     --hpo para reproducir sin el documento.
  2. Carga la ontologia (hp.obo) y las anotaciones gen->fenotipo de HPO
     (genes_to_phenotype.txt) y calcula la similitud semantica de cada gen
     con el paciente (ver src/hpo.py). El solapamiento exacto de terminos
     da 0 para casi todo porque las enfermedades se anotan con terminos
     hermanos o padres de los del paciente; la ontologia da credito parcial.
  3. Agrupa las variantes por gen, deduce la zigosidad y si el gen tiene
     evidencia bialelica (hom-alt o >= 2 het, candidatas a heterocigosis
     compuesta), que es lo esperable en MVA (autosomica recesiva).
  4. Chequeo explicito del panel MVA (BUB1B, CEP57, TRIP13, MAD2L1BP) y de
     genes del checkpoint del huso / centrosoma: se imprimen sus variantes
     aunque su score fenotipico sea bajo.

Salidas:
  output/track1/track1_step4_gene_ranking.csv   ranking por gen (sin coordenadas)
  data/track1_step4_candidates.tsv              variantes candidatas (contiene
                                                 coordenadas del paciente: data/
                                                 esta gitignoreado a proposito)

Uso:
    python track1_step4_prioritize.py
    python track1_step4_prioritize.py --input track1_step3b_regulatory_annotated.vcf.gz
    python track1_step4_prioritize.py --hpo HP:0002859 HP:0000121 ...
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

from src.hpo import Ontology, PhenotypeMatcher, load_gene_annotations  # noqa: E402

DEFAULT_INPUT = REPO_ROOT / "track1_step3a_snpeff_filter.vcf.gz"
DEFAULT_DOCX = REPO_ROOT / "data" / "Challenge_Clinical_Phenotype_1.docx"
DEFAULT_OBO = REPO_ROOT / "data" / "hp.obo"
DEFAULT_G2P = REPO_ROOT / "data" / "genes_to_phenotype.txt"
DEFAULT_RANKING = REPO_ROOT / "output" / "track1" / "track1_step4_gene_ranking.csv"
DEFAULT_CANDIDATES = REPO_ROOT / "data" / "track1_step4_candidates.tsv"

# Terminos extraidos en track1_step1_exp01-hpo_extraction.ipynb (docs/02-findings.md).
# Solo se usan si no hay .docx y no se pasa --hpo.
FALLBACK_HPO = ["HP:0002859", "HP:0000121", "HP:0004322", "HP:0001508",
                "HP:0003202", "HP:0001622", "HP:0001518", "HP:0200067"]

# Panel MVA: docs/00-references.md seccion 2.
MVA_PANEL = {"BUB1B": "MVA1 (OMIM:257300)", "CEP57": "MVA2 (OMIM:614114)",
             "TRIP13": "MVA3 (OMIM:617598)", "MAD2L1BP": "p31comet, aneuploidia biallelica"}
MVA_DISEASE_IDS = {"OMIM:257300", "OMIM:614114", "OMIM:617598", "ORPHA:1052"}

# Checkpoint del huso mitotico (SAC), cinetocoro y centrosoma: la maquinaria
# cuya perdida produce aneuploidia en mosaico. Lista curada, no exhaustiva.
SAC_CENTROSOME = {
    "BUB1", "BUB3", "MAD1L1", "MAD2L1", "TTK", "CDC20", "CENPE", "KNL1", "ZW10",
    "ZWILCH", "ZWINT", "SPDL1", "AURKB", "PLK1", "CENPF", "PCNT", "CDK5RAP2",
    "CEP152", "CEP135", "CEP63", "NDE1", "CENPJ", "KIF11", "ESPL1", "PTTG1",
    "SGO1", "CDCA8", "INCENP", "BIRC5", "NUF2", "NDC80", "SPC24", "SPC25",
    "MIS12", "DSN1", "ANAPC1", "CDC27", "STAG2", "SMC1A", "SMC3", "RAD21",
}

FLOAT_TAGS = ("GNOMAD_AF", "GNOMAD_EXOME_AF", "CADD_PHRED", "REVEL", "AB_FILTER",
              "AVI_SCORE", "AVI_PHRED", "AVI_QUANTILE")
STR_TAGS = ("GENE_SYMBOL", "CONSEQUENCE", "PUTATIVE_IMPACT", "CLINVAR_SIG",
            "CLINVAR_REVIEW_STATUS", "OMIM_PHENOTYPE", "REGULATORY_FEATURE", "AVI_STATUS")
IMPACT_RANK = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}


def log(msg):
    print(msg, flush=True)


def read_hpo_from_docx(path: Path) -> list[str]:
    import docx  # python-docx
    doc = docx.Document(str(path))
    terms = []
    for table in doc.tables:
        headers = None
        for i, row in enumerate(table.rows):
            cells = [c.text.strip().replace("\n", " ") for c in row.cells]
            if i == 0:
                headers = cells
                continue
            if headers and "HPO ID" in headers:
                val = cells[headers.index("HPO ID")]
                if val.startswith("HP:"):
                    terms.append(val)
    return terms


def parse_info(info: str) -> dict:
    out = {}
    for item in info.split(";"):
        if not item:
            continue
        k, sep, v = item.partition("=")
        out[k] = v if sep else True
    return out


def parse_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def zygosity_from_gt(gt: str) -> str:
    alleles = gt.replace("|", "/").split("/")
    if any(a == "." for a in alleles) or not alleles:
        return "UNKNOWN"
    if all(a == "0" for a in alleles):
        return "HOM_REF"
    if len(set(alleles)) == 1:
        return "HOM_ALT"
    return "HET"


def read_variants(vcf_path: Path):
    """Yield one dict per variant from the Step 3A/3B text VCF."""
    open_fn = gzip.open if str(vcf_path).endswith(".gz") else open
    with open_fn(vcf_path, "rt", encoding="utf-8") as fh:
        fmt_keys = None
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, _, ref, alt, _, _, info, fmt, sample = parts[:10]
            info_d = parse_info(info)
            fmt_keys = fmt.split(":")
            sample_d = dict(zip(fmt_keys, sample.split(":")))
            ad = sample_d.get("AD", "").split(",")
            rec = {
                "chrom": chrom[3:] if chrom.startswith("chr") else chrom,
                "pos": int(pos), "ref": ref, "alt": alt.split(",")[0],
                "gene": info_d.get("GENE_SYMBOL", ""),
                "zygosity": zygosity_from_gt(sample_d.get("GT", "./.")),
                "gt": sample_d.get("GT", ""),
                "dp": sample_d.get("DP", ""), "gq": sample_d.get("GQ", ""),
                "ad_ref": ad[0] if ad else "", "ad_alt": ad[1] if len(ad) > 1 else "",
            }
            for tag in STR_TAGS:
                rec[tag] = info_d.get(tag, "")
            for tag in FLOAT_TAGS:
                rec[tag] = parse_float(info_d.get(tag))
            yield rec


def main():
    ap = argparse.ArgumentParser(description="Paso 4: priorizacion de genes por HPO")
    ap.add_argument("--input", default=str(DEFAULT_INPUT), help="VCF del Paso 3A o 3B")
    ap.add_argument("--docx", default=str(DEFAULT_DOCX))
    ap.add_argument("--hpo", nargs="*", help="terminos HPO del paciente (anula el .docx)")
    ap.add_argument("--obo", default=str(DEFAULT_OBO))
    ap.add_argument("--g2p", default=str(DEFAULT_G2P))
    ap.add_argument("--ranking-out", default=str(DEFAULT_RANKING))
    ap.add_argument("--candidates-out", default=str(DEFAULT_CANDIDATES))
    ap.add_argument("--top", type=int, default=25, help="genes a mostrar en consola")
    args = ap.parse_args()

    for p, label in ((args.input, "VCF anotado"), (args.obo, "hp.obo"), (args.g2p, "genes_to_phenotype.txt")):
        if not Path(p).exists():
            sys.exit(f"[-] No se encontro {label}: {p}")

    # 1. Terminos HPO del paciente
    if args.hpo:
        hpo_terms, src = list(args.hpo), "--hpo"
    elif Path(args.docx).exists():
        hpo_terms, src = read_hpo_from_docx(Path(args.docx)), Path(args.docx).name
    else:
        hpo_terms, src = FALLBACK_HPO, "FALLBACK (docs/02-findings.md); falta el .docx"
    if not hpo_terms:
        sys.exit("[-] No se obtuvo ningun termino HPO")
    log(f"[1/4] {len(hpo_terms)} terminos HPO del paciente ({src}): {' '.join(hpo_terms)}")

    # 2. Ontologia + anotaciones
    onto = Ontology.from_obo(args.obo)
    gene_terms, gene_diseases = load_gene_annotations(args.g2p, onto)
    matcher = PhenotypeMatcher(onto, gene_terms)
    unknown = [t for t in hpo_terms if onto.canonical(t) is None]
    if unknown:
        log(f"      AVISO: terminos no encontrados en hp.obo y descartados: {unknown}")
    log(f"[2/4] Ontologia: {len(onto.parents):,} terminos | genes anotados: {len(gene_terms):,}")

    # 3. Variantes por gen
    by_gene: dict[str, list[dict]] = defaultdict(list)
    n_var = 0
    for rec in read_variants(Path(args.input)):
        n_var += 1
        by_gene[rec["gene"] or "(sin gen)"].append(rec)
    log(f"[3/4] {n_var:,} variantes en {len(by_gene):,} genes leidas de {args.input}")

    gene_rows = []
    for gene, variants in by_gene.items():
        n_het = sum(v["zygosity"] == "HET" for v in variants)
        n_hom = sum(v["zygosity"] == "HOM_ALT" for v in variants)
        biallelic = n_hom >= 1 or n_het >= 2
        pheno = matcher.score_gene(hpo_terms, gene)
        max_impact = max(variants, key=lambda v: IMPACT_RANK.get(v["PUTATIVE_IMPACT"], -1))["PUTATIVE_IMPACT"]
        diseases = gene_diseases.get(gene, set())
        matched_terms = ";".join(f"{t}={s}" for t, s, _ in pheno["matches"] if s > 0)
        gene_rows.append({
            "gene": gene,
            "pheno_score": round(pheno["score"], 4),
            "n_exact_hpo": pheno["n_exact"],
            "has_hpo_annotations": bool(pheno.get("has_annotations")),
            "n_variants": len(variants), "n_het": n_het, "n_hom_alt": n_hom,
            "biallelic_evidence": biallelic,
            "max_impact": max_impact,
            "mva_panel": gene in MVA_PANEL,
            "mva_disease_annotated": bool(diseases & MVA_DISEASE_IDS),
            "sac_centrosome": gene in SAC_CENTROSOME,
            "hpo_matches": matched_terms,
        })
        for v in variants:
            v["pheno_score"] = round(pheno["score"], 4)
            v["n_het_in_gene"] = n_het
            v["biallelic_evidence"] = biallelic
            v["mva_panel"] = gene in MVA_PANEL
            v["mva_disease_annotated"] = bool(diseases & MVA_DISEASE_IDS)
            v["sac_centrosome"] = gene in SAC_CENTROSOME

    gene_rows.sort(key=lambda r: (-r["pheno_score"], -int(r["biallelic_evidence"]),
                                  -IMPACT_RANK.get(r["max_impact"], -1), r["gene"]))

    # 4. Salidas
    Path(args.ranking_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.ranking_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(gene_rows[0].keys()) if gene_rows else ["gene"])
        w.writeheader()
        w.writerows(gene_rows)

    cand_fields = ["chrom", "pos", "ref", "alt", "gene", "zygosity", "gt", "dp", "gq", "ad_ref", "ad_alt",
                   *STR_TAGS, *FLOAT_TAGS, "pheno_score", "n_het_in_gene", "biallelic_evidence",
                   "mva_panel", "mva_disease_annotated", "sac_centrosome"]
    Path(args.candidates_out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.candidates_out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cand_fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for gene, variants in by_gene.items():
            for v in sorted(variants, key=lambda x: (x["chrom"], x["pos"])):
                w.writerow({k: ("" if v.get(k) is None else v.get(k)) for k in cand_fields})

    log(f"[4/4] Ranking de genes -> {args.ranking_out}")
    log(f"      Variantes candidatas -> {args.candidates_out}")

    log(f"\n--- TOP {args.top} GENES POR FENOTIPO ---")
    log(f"{'gene':<14}{'pheno':>7}{'exact':>6}{'nvar':>5}{'het':>4}{'hom':>4}{'biall':>6}  {'impact':<9}{'flags'}")
    for r in gene_rows[:args.top]:
        flags = ",".join(f for f, on in (("MVA", r["mva_panel"]), ("MVAdis", r["mva_disease_annotated"]),
                                          ("SAC", r["sac_centrosome"])) if on)
        log(f"{r['gene']:<14}{r['pheno_score']:>7.3f}{r['n_exact_hpo']:>6}{r['n_variants']:>5}"
            f"{r['n_het']:>4}{r['n_hom_alt']:>4}{'yes' if r['biallelic_evidence'] else 'no':>6}  "
            f"{r['max_impact']:<9}{flags}")

    log("\n--- CHEQUEO EXPLICITO DEL PANEL MVA ---")
    for gene, label in MVA_PANEL.items():
        variants = by_gene.get(gene, [])
        if not variants:
            log(f"{gene:<10} {label:<40} sin variantes HIGH/MODERATE raras")
            continue
        for v in variants:
            log(f"{gene:<10} {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']} {v['zygosity']:<8} "
                f"{v['PUTATIVE_IMPACT']:<9} {v['CONSEQUENCE']:<28} gnomAD={v['GNOMAD_AF']} "
                f"CADD={v['CADD_PHRED']} REVEL={v['REVEL']} AVI_PHRED={v['AVI_PHRED']} "
                f"ClinVar={v['CLINVAR_SIG'] or '-'}")

    sac_hits = [r for r in gene_rows if r["sac_centrosome"] and not r["mva_panel"]]
    if sac_hits:
        log("\n--- GENES SAC / CENTROSOMA CON VARIANTES ---")
        for r in sac_hits:
            log(f"{r['gene']:<10} pheno={r['pheno_score']:.3f} nvar={r['n_variants']} "
                f"biallelic={'yes' if r['biallelic_evidence'] else 'no'} impact={r['max_impact']}")


if __name__ == "__main__":
    main()
