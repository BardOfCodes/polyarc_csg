"""
Tests for polyset_parser.

Note: parse_csg_to_valid_polyset_csg requires a sketcher object from geolipi.
These tests are placeholders and require geolipi's Sketcher to be properly configured.
"""
import pytest

# Mark all tests as requiring geolipi sketcher setup
pytestmark = pytest.mark.skip(reason="Parser tests require geolipi Sketcher setup")


# TODO: Add parser tests when sketcher can be mocked or properly initialized
# The original tests in scripts/parse_test.py used an older API that didn't require sketcher.
# Current API: parse_csg_to_valid_polyset_csg(expression, sketcher, *args, **kwargs)

