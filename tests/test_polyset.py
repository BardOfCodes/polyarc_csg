"""Tests for polyset validation and conversion."""
import pytest
import polyarc_rs as prs
import geolipi.symbolic as gls
from polyarc_csg.polyset import is_valid_polyset, polyset_to_csg


def create_square(x: float, y: float, size: float, mode: int = 1) -> prs.PolyArc:
    """Creates a square PolyArc centered at (x, y)."""
    h = size / 2
    return prs.PolyArc(
        polyarc=((x - h, y - h, 0.0), (x + h, y - h, 0.0), (x + h, y + h, 0.0), (x - h, y + h, 0.0)),
        is_closed=True,
        mode=mode
    )


class TestValidPolyset:
    def test_single_square(self):
        square = create_square(0, 0, 10, 1)
        assert is_valid_polyset([square]) is True

    def test_two_disjoint_squares(self):
        sq1 = create_square(0, 0, 10, 1)
        sq2 = create_square(20, 20, 10, 1)
        assert is_valid_polyset([sq1, sq2]) is True

    def test_square_with_hole(self):
        outer = create_square(0, 0, 10, 1)
        inner = create_square(0, 0, 5, -1)
        assert is_valid_polyset([outer, inner]) is True

    def test_nested_with_two_holes(self):
        outer = create_square(0, 0, 50, 1)
        hole1 = create_square(-10, -10, 10, -1)
        hole2 = create_square(10, 10, 10, -1)
        assert is_valid_polyset([outer, hole1, hole2]) is True


class TestInvalidPolyset:
    def test_wrong_mode_alternation(self):
        outer = create_square(0, 0, 10, 1)
        inner = create_square(0, 0, 5, 1)  # Should be -1
        assert is_valid_polyset([outer, inner]) is False

    def test_overlapping_squares(self):
        sq1 = create_square(0, 0, 10, 1)
        sq2 = create_square(5, 5, 10, 1)
        assert is_valid_polyset([sq1, sq2]) is False

    def test_nested_wrong_alternation(self):
        outer = create_square(0, 0, 50, 1)
        hole1 = create_square(-10, -10, 10, -1)
        wrong = create_square(10, 10, 10, 1)  # Should be -1
        assert is_valid_polyset([outer, hole1, wrong]) is False


class TestPolysetToCsg:
    def test_single_square(self):
        square = create_square(0, 0, 10)
        csg = polyset_to_csg([square])
        assert isinstance(csg, gls.PolyArc2D)

    def test_square_with_hole(self):
        outer = create_square(0, 0, 10, mode=1)
        inner = create_square(0, 0, 5, mode=-1)
        csg = polyset_to_csg([outer, inner])
        assert isinstance(csg, gls.Difference)

    def test_two_disjoint_squares(self):
        sq1 = create_square(-20, 0, 10)
        sq2 = create_square(20, 0, 10)
        csg = polyset_to_csg([sq1, sq2])
        assert isinstance(csg, gls.Union)

    def test_nested_alternating(self):
        outer = create_square(0, 0, 20, mode=1)
        middle = create_square(0, 0, 15, mode=-1)
        inner = create_square(0, 0, 5, mode=1)
        csg = polyset_to_csg([outer, middle, inner])
        assert isinstance(csg, gls.Difference)

    def test_negative_outer(self):
        outer = create_square(0, 0, 20, mode=-1)
        hole = create_square(0, 0, 10, mode=1)
        csg = polyset_to_csg([outer, hole])
        assert isinstance(csg, gls.Complement)

    def test_complex_nesting(self):
        root = create_square(0, 0, 70, mode=1)
        hole_1 = create_square(-15, -15, 10, mode=-1)
        island_1 = create_square(-15, -15, 5, mode=1)
        hole_2 = create_square(20, 20, 10, mode=-1)
        csg = polyset_to_csg([root, hole_1, island_1, hole_2])
        assert isinstance(csg, gls.Difference)
