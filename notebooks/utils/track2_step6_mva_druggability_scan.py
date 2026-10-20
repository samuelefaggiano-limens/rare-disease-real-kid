"""Track 2, Step 6: druggability scan across every known MVA gene.

Scalability check for the Track 2 method. Step 3 found that none of the
proteins at STRING distance 1 from BUB1B has an approved, direction-concordant
drug. This script asks whether that druggability gap is specific to BUB1B or
general to mosaic variegated aneuploidy, by running the same lookup
automatically for every gene with a reported MVA / MVA-like phenotype.

For each MVA gene:
  1. take its high-confidence first-shell interactors from STRING
     (combined score >= 0.7, top N by score), plus the gene itself;
  2. resolve each protein to UniProt and query ChEMBL, through every target
     that contains it (single protein, complex, family), for drug mechanisms;
  3. keep molecules with max_phase 4 (approved), not withdrawn, collapsing
     salts to their parent.

Unlike Step 3, there is no curated panel and no desired direction of
modulation here: the scan is fully automatic and counts approved drugs in any
direction. It is therefore an upper bound on druggability, and the step that
needs expert input (choosing the direction of modulation) is kept out of it.

No patient data is used. Inputs are the gene list below and public APIs.

Usage:
    python notebooks/utils/track2_step6_mva_druggability_scan.py [--top 15]

Outputs (output/track2/):
  track2_step6_mva_scan_proteins.csv   one row per (MVA gene, protein)
  track2_step6_mva_scan_summary.csv    one row per MVA gene
"""
from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd

from track2_step3_drug_mapping import (STRING, STRING_MIN_SCORE, _session, chembl_mechanisms,
                                       chembl_molecule, chembl_targets, uniprot_accession)

# Genes with a reported MVA or MVA-like phenotype (docs/00-references.md, Track 1 panel).
MVA_GENES = ["BUB1B", "CEP57", "TRIP13", "MAD2L1BP", "CEP192", "SMC5", "BUB1", "CENATAC"]


def string_partners(gene: str, top: int) -> list[tuple[str, float]]:
    """High-confidence first-shell interactors of `gene`, best first."""
    for attempt in range(5):
        r = _session.get(f"{STRING}/interaction_partners",
                         params={"identifiers": gene, "species": 9606, "limit": top,
                                 "required_score": int(STRING_MIN_SCORE * 1000),
                                 "caller_identity": "mva_hackathon_track2"}, timeout=120)
        if r.status_code == 200:
            rows = r.json()
            return [(x["preferredName_B"], x["score"]) for x in rows
                    if x["score"] >= STRING_MIN_SCORE]
        time.sleep(2 ** attempt)
    raise RuntimeError(f"STRING failed for {gene}")


_protein_cache: dict[str, list[str]] = {}


def approved_drugs(symbol: str) -> list[str]:
    """Approved (max_phase 4, not withdrawn) parent drugs acting on any target containing `symbol`."""
    if symbol in _protein_cache:
        return _protein_cache[symbol]
    acc = uniprot_accession(symbol)
    drugs: set[str] = set()
    if acc:
        for t in chembl_targets(acc):
            for m in chembl_mechanisms(t["target_chembl_id"]):
                if not m.get("molecule_chembl_id"):
                    continue
                mol = chembl_molecule(m["molecule_chembl_id"])
                if mol["max_phase"] >= 4 and not mol["withdrawn"]:
                    drugs.add(mol["drug"])
    _protein_cache[symbol] = sorted(drugs)
    return _protein_cache[symbol]


def scan(top: int, out_dir: Path) -> pd.DataFrame:
    rows = []
    for gene in MVA_GENES:
        partners = string_partners(gene, top)
        print(f"{gene}: {len(partners)} first-shell partners", flush=True)
        for protein, score in [(gene, 1.0)] + partners:
            drugs = approved_drugs(protein)
            rows.append({"mva_gene": gene, "protein": protein,
                         "relation": "MVA gene" if protein == gene else "first-shell interactor",
                         "string_score": round(score, 3), "n_approved_drugs": len(drugs),
                         "approved_drugs": "; ".join(drugs)})
    prot = pd.DataFrame(rows)

    summ = []
    for gene, g in prot.groupby("mva_gene", sort=False):
        shell = g[g.relation == "first-shell interactor"]
        druggable = shell[shell.n_approved_drugs > 0]
        summ.append({
            "mva_gene": gene,
            "gene_itself_has_approved_drug": bool(g[g.relation == "MVA gene"].n_approved_drugs.iloc[0]),
            "first_shell_proteins": len(shell),
            "first_shell_with_approved_drug": len(druggable),
            "fraction_druggable": round(len(druggable) / len(shell), 2) if len(shell) else None,
            "druggable_proteins": "; ".join(f"{p} ({n})" for p, n in
                                            zip(druggable.protein, druggable.n_approved_drugs)),
        })
    summary = pd.DataFrame(summ)
    out_dir.mkdir(parents=True, exist_ok=True)
    prot.to_csv(out_dir / "track2_step6_mva_scan_proteins.csv", index=False)
    summary.to_csv(out_dir / "track2_step6_mva_scan_summary.csv", index=False)
    return summary


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--top", type=int, default=15, help="first-shell partners per MVA gene")
    ap.add_argument("--out", type=Path, default=here.parents[1] / "output" / "track2")
    a = ap.parse_args()
    t0 = time.time()
    s = scan(a.top, a.out)
    print(s.to_string(index=False))
    print(f"done in {(time.time() - t0) / 60:.1f} min")
