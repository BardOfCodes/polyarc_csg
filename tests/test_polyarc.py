"""Tests for polyarc boolean operations and utilities."""
import pytest
import polyarc_rs as prs
from polyarc_csg.polyarc import (
    union_multiple, 
    intersection_multiple, 
    difference_multiple,
    deduplicate_union, 
    deduplicate_intersection,
)


def create_square(x: float, y: float, size: float, mode: int = 1) -> prs.PolyArc:
    """Creates a square PolyArc centered at (x, y)."""
    h = size / 2
    return prs.PolyArc(
        polyarc=((x - h, y - h, 0.0), (x + h, y - h, 0.0), (x + h, y + h, 0.0), (x - h, y + h, 0.0)),
        is_closed=True,
        mode=mode
    )


class TestPolyArcCreation:
    def test_polyarc_creation(self):
        """Test if PolyArc is correctly created and stores its properties."""
        points = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0))
        poly = prs.PolyArc(points, is_closed=True, mode=1)
        
        assert poly.polyarc == points
        assert poly.is_closed is True
        assert poly.mode == 1


class TestBooleanOperations:
    def test_union(self):
        poly1 = prs.PolyArc(((0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0)), is_closed=True, mode=1)
        poly2 = prs.PolyArc(((1, 1, 0), (3, 1, 0), (3, 3, 0), (1, 3, 0)), is_closed=True, mode=1)
        result = prs.boolean_operation(poly1, poly2, "union")
        assert isinstance(result, list) and len(result) > 0

    def test_intersection(self):
        poly1 = prs.PolyArc(((0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0)), is_closed=True, mode=1)
        poly2 = prs.PolyArc(((1, 1, 0), (3, 1, 0), (3, 3, 0), (1, 3, 0)), is_closed=True, mode=1)
        result = prs.boolean_operation(poly1, poly2, "intersection")
        assert isinstance(result, list)

    def test_difference(self):
        poly1 = prs.PolyArc(((0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0)), is_closed=True, mode=1)
        poly2 = prs.PolyArc(((1, 1, 0), (3, 1, 0), (3, 3, 0), (1, 3, 0)), is_closed=True, mode=1)
        result = prs.boolean_operation(poly1, poly2, "difference")
        assert isinstance(result, list)

    def test_xor(self):
        poly1 = prs.PolyArc(((0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0)), is_closed=True, mode=1)
        poly2 = prs.PolyArc(((1, 1, 0), (3, 1, 0), (3, 3, 0), (1, 3, 0)), is_closed=True, mode=1)
        result = prs.boolean_operation(poly1, poly2, "xor")
        assert isinstance(result, list)


class TestSpatialRelationships:
    def test_inner_inside_outer(self):
        outer = prs.PolyArc(((0, 0, 0), (5, 0, 0), (5, 5, 0), (0, 5, 0)), is_closed=True, mode=1)
        inner = prs.PolyArc(((1, 1, 0), (2, 1, 0), (2, 2, 0), (1, 2, 0)), is_closed=True, mode=1)
        assert prs.is_1_inside_2(inner, outer) is True
        assert prs.is_1_inside_2(outer, inner) is False

    def test_disjoint(self):
        outer = prs.PolyArc(((0, 0, 0), (5, 0, 0), (5, 5, 0), (0, 5, 0)), is_closed=True, mode=1)
        disjoint = prs.PolyArc(((10, 10, 0), (12, 10, 0), (12, 12, 0), (10, 12, 0)), is_closed=True, mode=1)
        assert prs.is_disjoint(outer, disjoint) is True

    def test_not_intersecting(self):
        outer = prs.PolyArc(((0, 0, 0), (5, 0, 0), (5, 5, 0), (0, 5, 0)), is_closed=True, mode=1)
        inner = prs.PolyArc(((1, 1, 0), (2, 1, 0), (2, 2, 0), (1, 2, 0)), is_closed=True, mode=1)
        disjoint = prs.PolyArc(((10, 10, 0), (12, 10, 0), (12, 12, 0), (10, 12, 0)), is_closed=True, mode=1)
        assert prs.is_intersected(outer, inner) is False
        assert prs.is_intersected(outer, disjoint) is False


class TestUnionMultiple:
    def test_non_overlapping(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(20, 20, 10)
        result = union_multiple([sq1, sq2])
        assert len(result) == 2

    def test_overlapping_merge(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(5, 5, 10)
        result = union_multiple([sq1, sq2])
        assert len(result) == 1

    def test_nested_merge(self):
        sq1 = create_square(0, 0, 10)
        inner = create_square(0, 0, 5)
        result = union_multiple([sq1, inner])
        assert len(result) == 1

    def test_chain_of_overlaps(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(5, 5, 10)
        sq3 = create_square(10, 10, 10)
        result = union_multiple([sq1, sq2, sq3])
        assert len(result) == 1


class TestIntersectionMultiple:
    def test_non_overlapping_empty(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(20, 20, 10)
        result = intersection_multiple([sq1, sq2])
        assert len(result) == 0

    def test_fully_overlapping(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(0, 0, 10)
        result = intersection_multiple([sq1, sq2])
        assert len(result) == 1

    def test_partial_overlap(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(5, 5, 10)
        result = intersection_multiple([sq1, sq2])
        assert len(result) == 1

    def test_chain_overlap(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(5, 5, 10)
        sq3 = create_square(2, 2, 4)
        result = intersection_multiple([sq1, sq2, sq3])
        assert len(result) == 1


class TestDeduplicateUnion:
    def test_exact_duplicates(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(0, 0, 10)
        result = deduplicate_union([sq1, sq2])
        assert len(result) == 1

    def test_no_duplicates(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(20, 20, 10)
        result = deduplicate_union([sq1, sq2])
        assert len(result) == 2

    def test_enclosed_removed(self):
        outer = create_square(0, 0, 10)
        inner = create_square(0, 0, 5)
        result = deduplicate_union([outer, inner])
        assert len(result) == 1


class TestDeduplicateIntersection:
    def test_exact_duplicates(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(0, 0, 10)
        result = deduplicate_intersection([sq1, sq2])
        assert len(result) == 1

    def test_no_duplicates(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(20, 20, 10)
        result = deduplicate_intersection([sq1, sq2])
        assert len(result) == 2

    def test_overlapping_not_removed(self):
        sq1 = create_square(0, 0, 10)
        sq2 = create_square(5, 5, 10)
        result = deduplicate_intersection([sq1, sq2])
        assert len(result) == 2


class TestDifferenceMultiple:
    def test_no_overlap(self):
        a = create_square(0, 0, 10)
        b = create_square(20, 20, 5)
        result = difference_multiple([a], [b])
        assert len(result) == 1

    def test_complete_coverage(self):
        a = create_square(0, 0, 10)
        b = create_square(0, 0, 20)
        result = difference_multiple([a], [b])
        assert len(result) == 0

    def test_partial_overlap_splits(self):
        a = create_square(0, 0, 10)
        b = create_square(0, 0, 5)
        result = difference_multiple([a], [b])
        assert len(result) >= 1

    def test_multiple_bs(self):
        a = create_square(0, 0, 10)
        b1 = create_square(-2, 0, 4)
        b2 = create_square(2, 0, 4)
        result = difference_multiple([a], [b1, b2])
        assert len(result) >= 1

    def test_multiple_as_and_bs(self):
        a1 = create_square(-15, 0, 10)
        a2 = create_square(15, 0, 10)
        b1 = create_square(-15, 0, 5)
        b2 = create_square(15, 0, 5)
        result = difference_multiple([a1, a2], [b1, b2])
        assert len(result) >= 2

