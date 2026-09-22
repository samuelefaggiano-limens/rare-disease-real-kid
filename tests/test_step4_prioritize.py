import pytest
import math
from src.step4_prioritize import parse_info, parse_float, zygosity_from_gt

def test_parse_info():
    # Regular dict matching
    assert parse_info("GENE_SYMBOL=BUB1B;PUTATIVE_IMPACT=HIGH") == {"GENE_SYMBOL": "BUB1B", "PUTATIVE_IMPACT": "HIGH"}
    # Flag fields
    assert parse_info("GENE_SYMBOL=BUB1B;FLAGS") == {"GENE_SYMBOL": "BUB1B", "FLAGS": True}
    # Empty string
    assert parse_info("") == {}

def test_parse_float():
    assert parse_float("1.23") == 1.23
    assert parse_float("0.05") == 0.05
    assert math.isnan(parse_float("NaN"))
    assert parse_float("None") is None
    assert parse_float("") is None

def test_zygosity_from_gt():
    # Homozygous Reference
    assert zygosity_from_gt("0/0") == "HOM_REF"
    assert zygosity_from_gt("0|0") == "HOM_REF"
    # Heterozygous
    assert zygosity_from_gt("0/1") == "HET"
    assert zygosity_from_gt("1|0") == "HET"
    assert zygosity_from_gt("1/2") == "HET"
    # Homozygous Alternate
    assert zygosity_from_gt("1/1") == "HOM_ALT"
    assert zygosity_from_gt("2|2") == "HOM_ALT"
    # Unknown
    assert zygosity_from_gt("./.") == "UNKNOWN"
    assert zygosity_from_gt(".|.") == "UNKNOWN"
