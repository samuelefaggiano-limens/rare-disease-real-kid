import pytest
from src.step5_submission import fnum, fbool, variant_key, note_for

def test_fnum():
    assert fnum("1.23") == 1.23
    assert fnum("") is None
    assert fnum(None) is None
    assert fnum("None") is None
    assert fnum("NaN") is None or str(fnum("NaN")) == "nan"

def test_fbool():
    assert fbool("true") is True
    assert fbool("1") is True
    assert fbool("yes") is True
    assert fbool("True") is True
    
    assert fbool("false") is False
    assert fbool("0") is False
    assert fbool("no") is False
    assert fbool(None) is False
    assert fbool("") is False

def test_variant_key():
    v = {"chrom": "1", "pos": 100, "ref": "A", "alt": "T"}
    # With chr_prefix
    assert variant_key(v, True) == ("chr1", 100, "A", "T")
    # Without chr_prefix
    assert variant_key(v, False) == ("1", 100, "A", "T")

    # Already has chr
    v2 = {"chrom": "chr2", "pos": 200, "ref": "C", "alt": "G"}
    assert variant_key(v2, True) == ("chr2", 200, "C", "G")
    assert variant_key(v2, False) == ("chr2", 200, "C", "G")

def test_note_for():
    v = {
        "gene": "BUB1B",
        "CONSEQUENCE": "stop_gained",
        "PUTATIVE_IMPACT": "HIGH",
        "zygosity": "HET",
        "AVI_PHRED": "30.5"
    }
    s = {
        "components": {
            "phenotype": 0.8,
            "pathogenicity": 0.95,
            "inheritance": 1.0,
            "frequency": 0.9,
        },
        "bonus": 1.2
    }
    
    note = note_for(v, s)
    assert "BUB1B" in note
    assert "stop_gained" in note
    assert "HIGH" in note
    assert "pheno=0.80" in note
    assert "AVI_PHRED=30.5" in note
    assert "MVA-panel" in note

