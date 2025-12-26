"""Pytest configuration for polyarc_csg tests."""
import pytest


def pytest_configure(config):
    """Verify dependencies are importable before running tests."""
    try:
        import polyarc_rs
    except ImportError:
        pytest.exit("polyarc_rs not installed. Run: cd ../polyarc_rs && maturin develop --release")
    
    try:
        import geolipi.symbolic
    except ImportError:
        pytest.exit("geolipi not installed.")


@pytest.fixture
def polyarc_rs_module():
    """Fixture providing polyarc_rs module."""
    import polyarc_rs
    return polyarc_rs


@pytest.fixture
def geolipi_symbolic():
    """Fixture providing geolipi.symbolic module."""
    import geolipi.symbolic
    return geolipi.symbolic
