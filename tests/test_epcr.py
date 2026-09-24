from src.epcr import (CandidateRow, assign_epcr, combine, frequency_score, inheritance_score,
                      pathogenicity_score, score_against, what_if_table)


def test_pathogenicity_ordering():
    high = pathogenicity_score("HIGH", None, None, None)
    mod = pathogenicity_score("MODERATE", None, None, None)
    mod_cadd = pathogenicity_score("MODERATE", 32.0, None, None)
    benign = pathogenicity_score("MODERATE", 32.0, 0.9, "Benign")
    patho = pathogenicity_score("LOW", None, None, "Pathogenic")
    mod_avi = pathogenicity_score("MODERATE", None, None, None, avi_phred=31.0)
    assert high > mod
    assert mod_cadd > mod
    assert mod_avi > mod
    assert benign <= 0.10
    assert patho >= 0.95
    for s in (high, mod, mod_cadd, benign, patho):
        assert 0.0 <= s <= 1.0


def test_frequency_score_bounds():
    assert frequency_score(None) == 1.0
    assert frequency_score(0.0) == 1.0
    assert frequency_score(1e-6) == 1.0
    assert frequency_score(0.01) == 0.0
    assert 0.0 < frequency_score(1e-3) < frequency_score(1e-4) < 1.0


def test_inheritance_prefers_biallelic():
    assert inheritance_score("HOM_ALT", 0) == 1.0
    assert inheritance_score("HET", 2) == 1.0
    assert inheritance_score("HET", 1) < inheritance_score("HET", 2)


def test_combine_is_geometric_and_penalises_zero_axis():
    full = combine({"phenotype": 1, "pathogenicity": 1, "inheritance": 1, "frequency": 1})
    no_pheno = combine({"phenotype": 0, "pathogenicity": 1, "inheritance": 1, "frequency": 1})
    assert full == 1.0
    assert no_pheno < 0.3          # a zero axis costs multiplicatively, not additively
    assert combine({"phenotype": 1, "pathogenicity": 1, "inheritance": 1, "frequency": 1}, bonus=1.5) == 1.0


def _rows(scores):
    return [CandidateRow(((f"chr{i}", 1000 + i, "A", "G"),), s, f"G{i}") for i, s in enumerate(scores, 1)]


def test_assign_epcr_enforces_geometric_separation():
    ranked = assign_epcr(_rows([0.9, 0.89, 0.88, 0.87, 0.5, 0.4]), max_rows=5, tail_ratio=0.6)
    assert len(ranked) == 5
    epcrs = [r.epcr for r in ranked]
    assert epcrs == sorted(epcrs, reverse=True)
    for prev, cur in zip(epcrs, epcrs[1:]):
        assert cur <= prev * 0.6 + 1e-9
    assert all(0 < e <= 1 for e in epcrs)


def test_near_ties_would_kill_fmax_but_separation_saves_it():
    truth = frozenset({("chr2", 1002, "A", "G")})     # the answer is row 2
    tied = _rows([0.9, 0.9, 0.9])
    for r in tied:
        r.epcr = r.raw_score                            # no separation applied
    tied_result = score_against(tied, truth)
    separated = assign_epcr(_rows([0.9, 0.9, 0.9]), max_rows=3, tail_ratio=0.6)
    sep_result = score_against(separated, truth)
    assert tied_result["f_max"] < sep_result["f_max"]
    assert sep_result["rank_points"] == 50.0            # rank 2 => 50 points


def test_compound_het_partial_credit():
    pair = CandidateRow((("chr15", 100, "A", "G"), ("chr15", 200, "C", "T")), 0.9, "BUB1B")
    pair.epcr = 0.9
    truth = frozenset({("chr15", 100, "A", "G"), ("chr15", 999, "G", "A")})
    res = score_against([pair], truth)
    assert res["full_match_rank"] is None
    assert res["partial_match_rank"] == 1
    assert res["rank_points"] == 50.0


def test_what_if_table_row_one_is_perfect():
    ranked = assign_epcr(_rows([0.9, 0.7, 0.5]), max_rows=3)
    table = what_if_table(ranked)
    assert table[0]["rank_points"] == 100.0
    assert table[0]["f_max"] == 1.0
