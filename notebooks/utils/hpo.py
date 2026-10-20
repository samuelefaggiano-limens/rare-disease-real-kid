"""Ontology-aware phenotype matching between a patient's HPO terms and genes.

Used by Track 1, Step 4 (gene prioritisation). Pure functions over plain
dicts/sets so they can be unit-tested without the real ontology.

Method (Resnik-style semantic similarity, as in Phenolyzer/Exomiser):

  * Every gene's HPO annotations (``genes_to_phenotype.txt``) are propagated
    to all their ancestors in ``hp.obo`` (the "true path rule").
  * The information content of a term is ``IC(t) = -ln(n_t / N)``, where
    ``n_t`` is the number of genes annotated to ``t`` or any descendant and
    ``N`` is the number of annotated genes. Rare, specific terms score high;
    generic ones ("Abnormality of the skeletal system") score near zero.
  * For a patient term ``t`` and a gene ``G``, the match is the IC of the
    most informative common ancestor between ``t`` and any of ``G``'s terms,
    normalised by ``IC(t)`` so it lies in ``[0, 1]``: 1.0 means the gene is
    annotated to ``t`` itself (or something more specific), 0.0 means the
    only shared ancestor is the root.
  * The gene score is the mean over the patient's terms. It is therefore a
    *recall*-flavoured number: how much of the patient's presentation does
    this gene's known phenotype explain.

Why not exact term overlap: the patient's terms are specific (e.g.
HP:0001508 Failure to thrive) while disease annotations often use a sibling
or parent (HP:0001510 Growth delay). Exact overlap gives 0 for these; the
ontology gives partial credit proportional to how specific the shared
ancestor is.
"""

from __future__ import annotations

import math
from collections import defaultdict
from pathlib import Path
from typing import Iterable

HPO_ROOT = "HP:0000001"
PHENOTYPIC_ABNORMALITY = "HP:0000118"


class Ontology:
    """Minimal HPO ontology: parents, alt_id mapping, obsolete replacement."""

    def __init__(self, parents: dict[str, set[str]], alt_ids: dict[str, str] | None = None,
                 names: dict[str, str] | None = None):
        self.parents = parents
        self.alt_ids = alt_ids or {}
        self.names = names or {}
        self._anc_cache: dict[str, frozenset[str]] = {}

    @classmethod
    def from_obo(cls, path: str | Path) -> "Ontology":
        parents: dict[str, set[str]] = {}
        alt_ids: dict[str, str] = {}
        names: dict[str, str] = {}
        replaced_by: dict[str, str] = {}
        cur = None
        in_term = False
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line == "[Term]":
                    in_term, cur = True, None
                    continue
                if line.startswith("["):          # [Typedef] etc.
                    in_term = False
                    continue
                if not in_term or not line:
                    continue
                key, _, val = line.partition(": ")
                if key == "id":
                    cur = val.strip()
                    parents.setdefault(cur, set())
                elif cur is None:
                    continue
                elif key == "name":
                    names[cur] = val.strip()
                elif key == "alt_id":
                    alt_ids[val.strip()] = cur
                elif key == "is_a":
                    parents[cur].add(val.split("!")[0].strip())
                elif key == "replaced_by":
                    replaced_by[cur] = val.strip()
                elif key == "is_obsolete" and val.strip() == "true":
                    parents.pop(cur, None)
        # An obsolete id that names a replacement behaves like an alt_id.
        for old, new in replaced_by.items():
            alt_ids.setdefault(old, new)
        return cls(parents, alt_ids, names)

    def canonical(self, term: str) -> str | None:
        """Primary id for ``term`` (follows alt_id/replaced_by); None if unknown."""
        seen = 0
        while term not in self.parents and term in self.alt_ids and seen < 5:
            term = self.alt_ids[term]
            seen += 1
        return term if term in self.parents else None

    def ancestors(self, term: str) -> frozenset[str]:
        """``term`` plus every ancestor (inclusive closure)."""
        cached = self._anc_cache.get(term)
        if cached is not None:
            return cached
        out = {term}
        stack = list(self.parents.get(term, ()))
        while stack:
            p = stack.pop()
            if p in out:
                continue
            out.add(p)
            stack.extend(self.parents.get(p, ()))
        res = frozenset(out)
        self._anc_cache[term] = res
        return res

    def name(self, term: str) -> str:
        return self.names.get(term, term)


def load_gene_annotations(path: str | Path, ontology: Ontology
                          ) -> tuple[dict[str, set[str]], dict[str, set[str]]]:
    """Parse ``genes_to_phenotype.txt``.

    Returns ``(gene -> direct HPO terms, gene -> disease ids)``. Terms are
    canonicalised through the ontology; unknown ids are dropped.
    Header (2026 format): ncbi_gene_id, gene_symbol, hpo_id, hpo_name,
    frequency, disease_id.
    """
    gene_terms: dict[str, set[str]] = defaultdict(set)
    gene_diseases: dict[str, set[str]] = defaultdict(set)
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        i_gene, i_hpo = col["gene_symbol"], col["hpo_id"]
        i_dis = col.get("disease_id")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(i_gene, i_hpo):
                continue
            term = ontology.canonical(parts[i_hpo])
            if term is None:
                continue
            gene = parts[i_gene]
            gene_terms[gene].add(term)
            if i_dis is not None and len(parts) > i_dis and parts[i_dis]:
                gene_diseases[gene].add(parts[i_dis])
    return dict(gene_terms), dict(gene_diseases)


class PhenotypeMatcher:
    """Scores genes against a patient's HPO terms with Resnik similarity."""

    def __init__(self, ontology: Ontology, gene_terms: dict[str, set[str]]):
        self.onto = ontology
        # Propagate each gene's annotations to ancestors once.
        self.gene_closure: dict[str, frozenset[str]] = {}
        term_gene_count: dict[str, int] = defaultdict(int)
        for gene, terms in gene_terms.items():
            closure: set[str] = set()
            for t in terms:
                closure |= ontology.ancestors(t)
            fz = frozenset(closure)
            self.gene_closure[gene] = fz
            for t in fz:
                term_gene_count[t] += 1
        self.n_genes = max(len(gene_terms), 1)
        self._ic: dict[str, float] = {}
        self._ic_max = math.log(self.n_genes) if self.n_genes > 1 else 1.0
        for t, n in term_gene_count.items():
            self._ic[t] = -math.log(n / self.n_genes) if n else self._ic_max

    def ic(self, term: str) -> float:
        """Information content; a term no gene is annotated to gets the max."""
        return self._ic.get(term, self._ic_max)

    def term_gene_similarity(self, term: str, gene: str) -> tuple[float, str | None]:
        """Normalised Resnik similarity in [0, 1] and the best common ancestor."""
        closure = self.gene_closure.get(gene)
        if not closure:
            return 0.0, None
        ic_t = self.ic(term)
        if ic_t <= 0:
            return 0.0, None
        best, best_anc = 0.0, None
        for anc in self.onto.ancestors(term) & closure:
            v = self.ic(anc)
            if v > best:
                best, best_anc = v, anc
        return min(best / ic_t, 1.0), best_anc

    def score_gene(self, patient_terms: Iterable[str], gene: str) -> dict:
        """Mean normalised similarity over the patient's terms plus a breakdown."""
        terms = [t for t in (self.onto.canonical(x) for x in patient_terms) if t]
        if not terms:
            return {"score": 0.0, "n_terms": 0, "n_exact": 0, "matches": []}
        matches = []
        total = 0.0
        n_exact = 0
        for t in terms:
            sim, anc = self.term_gene_similarity(t, gene)
            total += sim
            if sim >= 0.999:
                n_exact += 1
            matches.append((t, round(sim, 3), anc))
        return {
            "score": total / len(terms),
            "n_terms": len(terms),
            "n_exact": n_exact,
            "matches": matches,
            "has_annotations": gene in self.gene_closure,
        }
