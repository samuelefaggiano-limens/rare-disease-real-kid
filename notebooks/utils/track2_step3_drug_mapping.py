"""Track 2, Step 3: map prioritised therapeutic targets to approved drugs.

Input: the target panel from Step 2 (`notebooks/utils/track2_target_panel.tsv`):
one row per gene with its relationship to BUB1B, the desired direction of
modulation and the therapeutic concept. No patient data is used in Track 2,
so every output here is safe to commit.

What it does, in order:
  1. Resolves each gene symbol to its reviewed human UniProt accession.
  2. ChEMBL (primary source): finds every ChEMBL target that contains that
     protein -- single proteins, but also complexes and families, because
     e.g. bortezomib is annotated against "26S proteasome", not PSMB5 -- and
     pulls the curated mechanisms of action (action type, max phase). Each
     molecule is then enriched with approval year, black-box warning, route
     and whether it is a parent compound.
  3. DGIdb (cross-check): approved-drug interactions per gene, with their
     interaction direction. Used to flag drugs ChEMBL does not list as a
     curated mechanism, not to override ChEMBL.
  4. STRING (supporting evidence only): combined interaction score of each
     target with BUB1B, and shortest-path distance to BUB1B in the network
     induced by the panel (edges with score >= 0.7).
  5. Direction check: a drug's action is "concordant" when it pushes the
     target in the desired direction (INHIBITOR/ANTAGONIST/... for
     `decrease`, ACTIVATOR/AGONIST/... for `increase`). Targets marked
     `context` have no desired direction; `caution` marks drugs the child may
     receive whose effect could be modified by the BUB1B defect.

The rule from the procedure applies downstream (Step 4): a drug is never
prioritised just because it hits a protein in the pathway. This step only
builds the evidence trail.

Outputs (output/track2/):
  track2_step3_drug_target_table.csv   one row per drug-target-mechanism
  track2_step3_target_summary.csv      one row per target (druggability gap)
  track2_step3_string_network.csv      edges among panel genes + BUB1B
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import networkx as nx
import pandas as pd
import requests

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT = "https://rest.uniprot.org/uniprotkb/search"
DGIDB = "https://dgidb.org/api/graphql"
STRING = "https://string-db.org/api/json"

DECREASE = {"INHIBITOR", "ANTAGONIST", "NEGATIVE ALLOSTERIC MODULATOR",
            "INVERSE AGONIST", "BLOCKER", "DEGRADER", "NEGATIVE MODULATOR",
            "ANTISENSE INHIBITOR", "RNAI INHIBITOR", "DISRUPTING AGENT",
            "CROSS-LINKING AGENT", "STABILISER"}
INCREASE = {"ACTIVATOR", "AGONIST", "POSITIVE ALLOSTERIC MODULATOR",
            "POSITIVE MODULATOR", "PARTIAL AGONIST", "OPENER", "RELEASING AGENT"}

STRING_MIN_SCORE = 0.7

_session = requests.Session()


def _get(url: str, params: dict | None = None, retries: int = 6) -> dict:
    """GET with exponential backoff; the EBI APIs throttle bursts."""
    for attempt in range(retries):
        try:
            r = _session.get(url, params=params, timeout=60)
            if r.status_code == 200:
                return r.json()
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    raise RuntimeError(f"GET failed: {url} {params}")


_ACCESSION_CACHE = Path(__file__).resolve().parent / "track2_uniprot_accessions.json"


def uniprot_accession(symbol: str) -> str | None:
    """Reviewed human accession for a gene symbol, cached on disk.

    UniProt throttles bursts and intermittently drops requests, which would
    otherwise abort a whole run partway through the panel. The cache makes
    re-runs offline-stable and reproducible.
    """
    cache = json.loads(_ACCESSION_CACHE.read_text()) if _ACCESSION_CACHE.exists() else {}
    if symbol in cache:
        return cache[symbol]
    q = f"gene_exact:{symbol} AND organism_id:9606 AND reviewed:true"
    res = _get(UNIPROT, {"query": q, "fields": "accession", "format": "json", "size": 1})
    hits = res.get("results", [])
    cache[symbol] = hits[0]["primaryAccession"] if hits else None
    _ACCESSION_CACHE.write_text(json.dumps(cache, indent=2, sort_keys=True))
    return cache[symbol]


def chembl_targets(accession: str) -> list[dict]:
    """All ChEMBL targets (single protein, complex, family) containing the protein."""
    res = _get(f"{CHEMBL}/target.json",
               {"target_components__accession": accession,
                "organism": "Homo sapiens", "limit": 100})
    return [{"target_chembl_id": t["target_chembl_id"],
             "target_pref_name": t["pref_name"],
             "target_type": t["target_type"]} for t in res["targets"]]


def chembl_mechanisms(target_chembl_id: str) -> list[dict]:
    res = _get(f"{CHEMBL}/mechanism.json",
               {"target_chembl_id": target_chembl_id, "limit": 1000})
    return res["mechanisms"]


_molecule_cache: dict[str, dict] = {}


def chembl_molecule(chembl_id: str) -> dict:
    if chembl_id not in _molecule_cache:
        m = _get(f"{CHEMBL}/molecule/{chembl_id}.json")
        parent = (m.get("molecule_hierarchy") or {}).get("parent_chembl_id") or chembl_id
        if parent != chembl_id:
            # Salts/hydrates (e.g. IXAZOMIB CITRATE) are reported under their parent.
            p = chembl_molecule(parent)
            m["max_phase"] = max(float(m.get("max_phase") or 0), p["max_phase"])
            m["first_approval"] = m.get("first_approval") or p["first_approval"]
            m["pref_name"] = p["drug"]
        _molecule_cache[chembl_id] = {
            "drug": m.get("pref_name") or chembl_id,
            "parent_chembl_id": parent,
            "max_phase": float(m["max_phase"]) if m.get("max_phase") else 0.0,
            "first_approval": m.get("first_approval"),
            "black_box_warning": bool(m.get("black_box_warning")),
            "oral": bool(m.get("oral")),
            "withdrawn": bool(m.get("withdrawn_flag")),
        }
    return _molecule_cache[chembl_id]


def dgidb_approved(symbols: list[str]) -> pd.DataFrame:
    q = """query($names: [String!]) { genes(names: $names) { nodes { name
             interactions { drug { name approved } interactionScore
               interactionTypes { type directionality } } } } }"""
    res = _session.post(DGIDB, json={"query": q, "variables": {"names": symbols}},
                        timeout=120).json()
    rows = []
    for g in res["data"]["genes"]["nodes"]:
        for i in g["interactions"]:
            if not i["drug"]["approved"]:
                continue
            dirs = {t["directionality"] for t in i["interactionTypes"] if t["directionality"]}
            types = sorted({t["type"] for t in i["interactionTypes"]})
            rows.append({"gene": g["name"], "drug": i["drug"]["name"].upper(),
                         "dgidb_types": ";".join(types),
                         "dgidb_direction": ";".join(sorted(dirs)) or "unknown",
                         "dgidb_score": round(i["interactionScore"], 3)})
    return pd.DataFrame(rows, columns=["gene", "drug", "dgidb_types",
                                       "dgidb_direction", "dgidb_score"])


def string_network(symbols: list[str], anchor: str = "BUB1B") -> tuple[pd.DataFrame, dict, dict]:
    """Edges among the panel; per-gene score with and distance to the anchor."""
    ids = sorted(set(symbols) | {anchor})
    res = _session.post(f"{STRING}/network",
                        data={"identifiers": "\r".join(ids), "species": 9606,
                              "caller_identity": "mva_hackathon_track2"},
                        timeout=120).json()
    edges = pd.DataFrame([{"a": e["preferredName_A"], "b": e["preferredName_B"],
                           "score": e["score"]} for e in res]).drop_duplicates()
    g = nx.Graph()
    g.add_nodes_from(ids)
    strong = edges[edges["score"] >= STRING_MIN_SCORE]
    g.add_edges_from(zip(strong["a"], strong["b"]))
    dist = nx.single_source_shortest_path_length(g, anchor)
    direct = {}
    for _, e in edges.iterrows():
        if anchor in (e["a"], e["b"]):
            other = e["b"] if e["a"] == anchor else e["a"]
            direct[other] = max(direct.get(other, 0), e["score"])
    return edges, direct, dist


def direction_concordance(desired: str, action: str) -> str:
    action = (action or "").upper()
    if desired == "decrease":
        return "concordant" if action in DECREASE else "discordant"
    if desired == "increase":
        return "concordant" if action in INCREASE else "discordant"
    if desired == "caution":
        return "caution"
    return "no_desired_direction"


def build(panel_path: Path, out_dir: Path, verbose: bool = True) -> dict[str, pd.DataFrame]:
    panel = pd.read_csv(panel_path, sep="\t")
    rows, summary = [], []
    for _, t in panel.iterrows():
        acc = uniprot_accession(t.gene)
        targets = chembl_targets(acc) if acc else []
        n_mech = 0
        for tg in targets:
            for m in chembl_mechanisms(tg["target_chembl_id"]):
                if not m.get("molecule_chembl_id"):
                    continue
                mol = chembl_molecule(m["molecule_chembl_id"])
                n_mech += 1
                rows.append({
                    "gene": t.gene, "uniprot": acc,
                    "relationship_class": t.relationship_class,
                    "desired_direction": t.desired_direction,
                    "therapeutic_concept": t.therapeutic_concept,
                    "child_phenotype": t.child_phenotype,
                    **tg,
                    "molecule_chembl_id": m["molecule_chembl_id"],
                    **mol,
                    "action_type": m["action_type"],
                    "mechanism_of_action": m["mechanism_of_action"],
                    "direction": direction_concordance(t.desired_direction, m["action_type"]),
                    "source": "ChEMBL",
                })
        summary.append({"gene": t.gene, "uniprot": acc, "chembl_targets": len(targets),
                        "chembl_mechanisms": n_mech})
        if verbose:
            print(f"{t.gene:9s} {acc or '-':8s} targets={len(targets):2d} mechanisms={n_mech}")

    table = pd.DataFrame(rows)
    # One row per (gene, parent drug, target, action): salts and parents collapse.
    table = (table.sort_values("max_phase", ascending=False)
                  .drop_duplicates(["gene", "parent_chembl_id", "target_chembl_id", "action_type"]))
    table["approved"] = table["max_phase"] >= 4

    dg = dgidb_approved(panel.gene.tolist())
    table = table.merge(dg, on=["gene", "drug"], how="left")
    table["in_dgidb_approved"] = table["dgidb_types"].notna()

    edges, direct, dist = string_network(panel.gene.tolist())
    for df in (table,):
        df["string_score_bub1b"] = df["gene"].map(direct).fillna(0.0)
        df["string_distance_bub1b"] = df["gene"].map(dist)

    summ = pd.DataFrame(summary).merge(panel[["gene", "relationship_class", "desired_direction",
                                              "therapeutic_concept", "child_phenotype"]], on="gene")
    appr = table[table.approved]
    summ["approved_drugs_any_direction"] = summ.gene.map(
        appr.groupby("gene")["drug"].nunique()).fillna(0).astype(int)
    summ["approved_drugs_concordant"] = summ.gene.map(
        appr[appr.direction == "concordant"].groupby("gene")["drug"].nunique()).fillna(0).astype(int)
    summ["approved_concordant_examples"] = summ.gene.map(
        appr[appr.direction == "concordant"].groupby("gene")["drug"]
        .agg(lambda s: "; ".join(sorted(set(s))[:6]))).fillna("")
    summ["dgidb_approved_drugs"] = summ.gene.map(dg.groupby("gene")["drug"].nunique()).fillna(0).astype(int)
    summ["string_score_bub1b"] = summ.gene.map(direct).fillna(0.0)
    summ["string_distance_bub1b"] = summ.gene.map(dist)

    out_dir.mkdir(parents=True, exist_ok=True)
    cols = ["gene", "relationship_class", "desired_direction", "therapeutic_concept",
            "child_phenotype", "drug", "molecule_chembl_id", "parent_chembl_id", "approved",
            "max_phase", "first_approval", "action_type", "direction", "mechanism_of_action",
            "target_chembl_id", "target_pref_name", "target_type", "black_box_warning", "oral",
            "withdrawn", "in_dgidb_approved", "dgidb_direction", "dgidb_score",
            "string_score_bub1b", "string_distance_bub1b", "uniprot", "source"]
    table = table[cols].sort_values(["gene", "approved", "max_phase"], ascending=[True, False, False])
    table.to_csv(out_dir / "track2_step3_drug_target_table.csv", index=False)
    summ.to_csv(out_dir / "track2_step3_target_summary.csv", index=False)
    edges.to_csv(out_dir / "track2_step3_string_network.csv", index=False)
    return {"table": table, "summary": summ, "edges": edges, "dgidb": dg}


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--panel", type=Path, default=here / "track2_target_panel.tsv")
    ap.add_argument("--out", type=Path, default=here.parents[1] / "output" / "track2")
    a = ap.parse_args()
    build(a.panel, a.out)
