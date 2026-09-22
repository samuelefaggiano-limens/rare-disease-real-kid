import pytest
from src.step3a_annotate import extract_ann, format_hgvs_id, parse_metadata

def test_extract_ann():
    # Only ANN is present
    assert extract_ann("ANN=T|missense_variant|MODERATE|BUB1B") == "T|missense_variant|MODERATE|BUB1B"
    # Other info before ANN
    assert extract_ann("DP=10;ANN=T|missense_variant|MODERATE|BUB1B;GQ=20") == "T|missense_variant|MODERATE|BUB1B"
    # ANN at the end
    assert extract_ann("DP=10;ANN=T|missense_variant|MODERATE|BUB1B") == "T|missense_variant|MODERATE|BUB1B"
    # No ANN
    assert extract_ann("DP=10;GQ=20") is None

def test_format_hgvs_id():
    # SNV
    assert format_hgvs_id("1", "100", "A", "T") == "chr1:g.100A>T"
    assert format_hgvs_id("chr1", "100", "A", "T") == "chr1:g.100A>T"
    # Deletion (1 base)
    assert format_hgvs_id("1", "100", "AT", "A") == "chr1:g.101del"
    # Deletion (multiple bases)
    assert format_hgvs_id("1", "100", "ATC", "A") == "chr1:g.101_102del"
    # Insertion
    assert format_hgvs_id("1", "100", "A", "AT") == "chr1:g.100_101insT"
    assert format_hgvs_id("1", "100", "A", "ATC") == "chr1:g.100_101insTC"
    # Delins with anchor
    assert format_hgvs_id("1", "100", "ATC", "AGC") == "chr1:g.101_102delinsGC"
    # Delins without anchor
    assert format_hgvs_id("1", "100", "TC", "GC") == "chr1:g.100_101delinsGC"

def test_parse_metadata_empty():
    assert parse_metadata(None) == {}
    assert parse_metadata({"notfound": True}) == {}

def test_parse_metadata_gnomad():
    data = {"gnomad_genome": {"af": {"af": 0.05}}}
    res = parse_metadata(data)
    assert res.get("GNOMAD_AF") == "0.05"

def test_parse_metadata_clinvar():
    data = {
        "clinvar": {
            "rcv": {
                "clinical_significance": "Pathogenic",
                "review_status": "criteria provided, single submitter",
                "conditions": {"identifiers": {"omim": "123456"}}
            }
        }
    }
    res = parse_metadata(data)
    assert res.get("CLINVAR_SIG") == "Pathogenic"
    assert res.get("CLINVAR_REVIEW_STATUS") == "criteria_provided_single_submitter"
    assert res.get("OMIM_PHENOTYPE") == "123456"

    # Conflicting / multiple
    data_multi = {
        "clinvar": {
            "rcv": [
                {"clinical_significance": "Benign"},
                {"clinical_significance": "Pathogenic", "review_status": "reviewed by expert panel"}
            ]
        }
    }
    res_multi = parse_metadata(data_multi)
    assert res_multi.get("CLINVAR_SIG") == "Benign|Pathogenic"
    assert res_multi.get("CLINVAR_REVIEW_STATUS") == "reviewed_by_expert_panel"
