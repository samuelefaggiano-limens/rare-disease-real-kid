import pytest
from src.step2_prefilter import variant_keys

def test_variant_keys():
    # Simple SNV
    assert variant_keys("1", 100, "A", ["T"]) == ["1:100:A:T"]
    # With chr prefix
    assert variant_keys("chr1", 100, "A", ["T"]) == ["1:100:A:T"]
    # Multiple ALTs
    assert variant_keys("chr1", 100, "A", ["T", "C"]) == ["1:100:A:T", "1:100:A:C"]
    
