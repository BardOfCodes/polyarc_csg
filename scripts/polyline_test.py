import polyline_rs as prs
from polyline_csg.polyline import union_multiple, intersection_multiple, deduplicate_intersection, deduplicate_union, difference_multiple

def create_square(x, y, size, mode=1):
    """Creates a square PolyStruct centered at (x, y) with the given size."""
    half_size = size / 2
    polyline = [
        (x - half_size, y - half_size, 0.0),
        (x + half_size, y - half_size, 0.0),
        (x + half_size, y + half_size, 0.0),
        (x - half_size, y + half_size, 0.0),
    ]
    return prs.PolyStruct(polyline=polyline, is_closed=True, mode=mode)


def test_polystruct_creation():
    """Test if PolyStruct is correctly created and stores its properties."""
    points = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)]
    poly = prs.PolyStruct(points, is_closed=True, mode=1)
    
    assert poly.polyline == points
    assert poly.is_closed is True
    assert poly.mode == 1


def test_boolean_operations():
    """Test union, intersection, difference, and XOR between polylines."""
    poly1 = prs.PolyStruct([(0, 0, 0), (2, 0, 0), (2, 2, 0), (0, 2, 0)], is_closed=True, mode=1)
    poly2 = prs.PolyStruct([(1, 1, 0), (3, 1, 0), (3, 3, 0), (1, 3, 0)], is_closed=True, mode=1)

    # Union
    result_union = prs.boolean_operation(poly1, poly2, "union")
    assert isinstance(result_union, list) and len(result_union) > 0, "Union failed"

    # Intersection
    result_intersection = prs.boolean_operation(poly1, poly2, "intersection")
    assert isinstance(result_intersection, list), "Intersection failed"

    # Difference
    result_difference = prs.boolean_operation(poly1, poly2, "difference")
    assert isinstance(result_difference, list), "Difference failed"

    # XOR
    result_xor = prs.boolean_operation(poly1, poly2, "xor")
    assert isinstance(result_xor, list), "XOR failed"

def test_spatial_relationships():
    """Test spatial relationships between polylines."""
    outer = prs.PolyStruct([(0, 0, 0), (5, 0, 0), (5, 5, 0), (0, 5, 0)], is_closed=True, mode=1)
    inner = prs.PolyStruct([(1, 1, 0), (2, 1, 0), (2, 2, 0), (1, 2, 0)], is_closed=True, mode=1)
    disjoint = prs.PolyStruct([(10, 10, 0), (12, 10, 0), (12, 12, 0), (10, 12, 0)], is_closed=True, mode=1)

    assert prs.is_1_inside_2(inner, outer) is True, "Inner polyline should be inside outer"
    assert prs.is_1_inside_2(outer, inner) is False, "Outer should not be inside inner"
    assert prs.is_disjoint(outer, disjoint) is True, "Polylines should be disjoint"
    assert prs.is_intersected(outer, inner) is False, "Should not be intersecting"
    assert prs.is_intersected(outer, disjoint) is False, "Should not be intersecting"

# ====================== TESTING union_multiple ======================


# ====================== Test for union_multiple ======================

def test_union_multiple():
    print("Running union_multiple tests...")

    # Test 1: Non-overlapping squares
    square1 = create_square(0, 0, 10)
    square2 = create_square(20, 20, 10)
    result = union_multiple([square1, square2])
    assert len(result) == 2, "Test 1 Failed: Should result in 2 non-overlapping polylines"

    # Test 2: Overlapping squares (should merge into 1)
    square3 = create_square(5, 5, 10)
    result = union_multiple([square1, square3])
    assert len(result) == 1, "Test 2 Failed: Overlapping squares should merge into 1 polyline"

    # Test 3: Nested squares
    inner_square = create_square(0, 0, 5)
    result = union_multiple([square1, inner_square])
    assert len(result) == 1, "Test 3 Failed: Nested squares should merge into 1 polyline"

    # Test 4: Chain of overlapping squares
    square4 = create_square(10, 10, 10)
    result = union_multiple([square1, square4, square3])
    assert len(result) == 1, "Test 4 Failed: Chain of overlaps should result in 1 polyline"

    print("All union_multiple tests passed!\n")


# ====================== Test for intersection_multiple ======================

def test_intersection_multiple():
    print("Running intersection_multiple tests...")

    # Test 1: Non-overlapping squares
    square1 = create_square(0, 0, 10)
    square2 = create_square(20, 20, 10)
    result = intersection_multiple([square1, square2])
    assert len(result) == 0, "Test 1 Failed: Non-overlapping squares should have no intersection"

    # Test 2: Fully overlapping squares
    square3 = create_square(0, 0, 10)
    result = intersection_multiple([square1, square3])
    assert len(result) == 1, "Test 2 Failed: Fully overlapping squares should have 1 intersection"

    # Test 3: Partial overlap
    square4 = create_square(5, 5, 10)
    result = intersection_multiple([square1, square4])
    assert len(result) == 1, "Test 3 Failed: Partial overlap should result in 1 intersection region"

    # Test 4: Chain of overlapping squares (intersection reduces over time)
    square5 = create_square(2, 2, 4)
    result = intersection_multiple([square1, square4, square5])
    assert len(result) == 1, "Test 4 Failed: Chain overlap should result in 1 small intersected region"

    print("All intersection_multiple tests passed!\n")


# ====================== Test for deduplicate_union ======================

def test_deduplicate_union():
    print("Running deduplicate_union tests...")

    # Test 1: Exact duplicates
    square1 = create_square(0, 0, 10)
    square_dup = create_square(0, 0, 10)
    result = deduplicate_union([square1, square_dup])
    assert len(result) == 1, "Test 1 Failed: Duplicate squares should be deduplicated"

    # Test 2: No duplicates
    square2 = create_square(20, 20, 10)
    result = deduplicate_union([square1, square2])
    assert len(result) == 2, "Test 2 Failed: No duplicates should remain both"

    # Test 3: Enclosed polylines
    inner_square = create_square(0, 0, 5)
    result = deduplicate_union([square1, inner_square])
    assert len(result) == 1, "Test 3 Failed: Enclosed polyline should be removed"

    print("All deduplicate_union tests passed!\n")


# ====================== Test for deduplicate_intersection ======================

def test_deduplicate_intersection():
    print("Running deduplicate_intersection tests...")

    # Test 1: Exact duplicates
    square1 = create_square(0, 0, 10)
    square_dup = create_square(0, 0, 10)
    result = deduplicate_intersection([square1, square_dup])
    assert len(result) == 1, "Test 1 Failed: Duplicate squares should be deduplicated"

    # Test 2: No duplicates
    square2 = create_square(20, 20, 10)
    result = deduplicate_intersection([square1, square2])
    assert len(result) == 2, "Test 2 Failed: No duplicates should remain both"

    # Test 3: Overlapping but not identical
    square3 = create_square(5, 5, 10)
    result = deduplicate_intersection([square1, square3])
    assert len(result) == 2, "Test 3 Failed: Overlapping (non-identical) should not be removed"

    print("All deduplicate_intersection tests passed!\n")


def test_difference_multiple():
    # Case 1: No overlap (A remains unchanged)
    a = create_square(0, 0, 10)
    b = create_square(20, 20, 5)
    result = difference_multiple([a], [b])
    assert len(result) == 1, "No overlap - A should remain intact"

    # Case 2: Complete coverage (A disappears)
    a = create_square(0, 0, 10)
    b = create_square(0, 0, 20)
    result = difference_multiple([a], [b])
    assert len(result) == 0, "Complete coverage - A should be removed"

    # Case 3: Partial overlap (splitting A)
    a = create_square(0, 0, 10)
    b = create_square(0, 0, 5)  # Smaller square inside A
    result = difference_multiple([a], [b])
    assert len(result) > 1, "Partial overlap - A should split into multiple parts"

    # Case 4: Multiple Bs subtracting from A
    a = create_square(0, 0, 10)
    b1 = create_square(-2, 0, 4)
    b2 = create_square(2, 0, 4)
    result = difference_multiple([a], [b1, b2])
    assert len(result) >= 2, "Multiple Bs should carve out two holes in A"

    # Case 5: Multiple As and Bs
    a1 = create_square(-15, 0, 10)
    a2 = create_square(15, 0, 10)
    b1 = create_square(-15, 0, 5)
    b2 = create_square(15, 0, 5)
    result = difference_multiple([a1, a2], [b1, b2])
    assert len(result) >= 2, "Both A1 and A2 should have subtractions"

    print("All test cases passed!")

if __name__ == "__main__":
    test_polystruct_creation()
    test_boolean_operations()
    test_spatial_relationships()
    test_union_multiple()
    test_intersection_multiple()
    test_deduplicate_union()
    test_deduplicate_intersection()
    test_difference_multiple()
    print("All polyline_rs tests passed successfully!")