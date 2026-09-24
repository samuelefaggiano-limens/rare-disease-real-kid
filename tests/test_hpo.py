import math

from src.hpo import Ontology, PhenotypeMatcher

# Tiny ontology:
#   ROOT
#    └─ ABN (phenotypic abnormality)
#        ├─ GROWTH
#        │    ├─ SHORT
#        │    └─ FTT  (failure to thrive)
#        └─ TUMOR
#             └─ RMS (rhabdomyosarcoma)
PARENTS = {
    "HP:ROOT": set(), "HP:ABN": {"HP:ROOT"}, "HP:GROWTH": {"HP:ABN"},
    "HP:SHORT": {"HP:GROWTH"}, "HP:FTT": {"HP:GROWTH"}, "HP:TUMOR": {"HP:ABN"},
    "HP:RMS": {"HP:TUMOR"},
}
ALT = {"HP:OLD_RMS": "HP:RMS"}
GENES = {
    "GENE_A": {"HP:SHORT", "HP:RMS"},   # matches patient exactly
    "GENE_B": {"HP:FTT"},               # sibling of SHORT under GROWTH
    "GENE_C": {"HP:TUMOR"},             # parent of RMS
    "GENE_D": {"HP:ABN"},               # generic only
}
PATIENT = ["HP:SHORT", "HP:OLD_RMS"]    # one alt_id on purpose


def make():
    onto = Ontology(PARENTS, ALT)
    return onto, PhenotypeMatcher(onto, GENES)


def test_ancestors_inclusive():
    onto, _ = make()
    assert onto.ancestors("HP:RMS") == frozenset({"HP:RMS", "HP:TUMOR", "HP:ABN", "HP:ROOT"})


def test_alt_id_canonicalised():
    onto, _ = make()
    assert onto.canonical("HP:OLD_RMS") == "HP:RMS"
    assert onto.canonical("HP:NOPE") is None


def test_ic_monotone_with_specificity():
    _, m = make()
    # every gene reaches ROOT and ABN => IC 0; SHORT only GENE_A => max IC
    assert m.ic("HP:ROOT") == 0.0
    assert m.ic("HP:ABN") == 0.0
    assert math.isclose(m.ic("HP:SHORT"), -math.log(1 / 4))
    assert m.ic("HP:GROWTH") < m.ic("HP:SHORT")


def test_exact_match_scores_one():
    _, m = make()
    res = m.score_gene(PATIENT, "GENE_A")
    assert math.isclose(res["score"], 1.0)
    assert res["n_exact"] == 2


def test_sibling_gets_partial_credit_and_parent_more_than_root():
    _, m = make()
    b = m.score_gene(PATIENT, "GENE_B")["score"]   # FTT ~ SHORT share GROWTH
    c = m.score_gene(PATIENT, "GENE_C")["score"]   # TUMOR is parent of RMS
    d = m.score_gene(PATIENT, "GENE_D")["score"]   # only ABN (IC 0)
    assert 0 < b < 1
    assert 0 < c < 1
    assert d == 0.0


def test_unknown_gene_scores_zero_and_is_flagged():
    _, m = make()
    res = m.score_gene(PATIENT, "NOT_A_GENE")
    assert res["score"] == 0.0
    assert res["has_annotations"] is False
