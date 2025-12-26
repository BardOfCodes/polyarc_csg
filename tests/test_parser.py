"""
Tests for polyset_parser.

Note: parse_csg_to_valid_polyset_csg requires a sketcher object from geolipi.
These tests are placeholders and require geolipi's Sketcher to be properly configured.
"""
import pytest

# Mark all tests as requiring geolipi sketcher setup
pytestmark = pytest.mark.skip(reason="Parser tests require geolipi Sketcher setup")


# TODO: Add random CSG generation. Then 
# Test that occupancy matches. 
# Test that SDF is within error tolerance. For creating GT use Scipy. 


