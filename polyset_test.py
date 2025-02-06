import polyline_rs as prs
from polyline_csg.polyset import is_valid_polyset, polyset_to_csg
import geolipi.symbolic as gls
import woodie.symbolic as ws

def create_square(x, y, size, mode=1):
    """Creates a square PolyStruct centered at (x, y) with the given size and mode."""
    half_size = size / 2
    polyline = [
        (x - half_size, y - half_size, 0.0),
        (x + half_size, y - half_size, 0.0),
        (x + half_size, y + half_size, 0.0),
        (x - half_size, y + half_size, 0.0),
    ]
    return prs.PolyStruct(polyline=polyline, is_closed=True, mode=mode)

def test_valid_polysets():
    """Test cases for valid PolySets."""

    # Case 1: A single square (always valid)
    square_1 = create_square(0, 0, 10, 1)
    assert is_valid_polyset([square_1]) is True, "Test Case 1 Failed"

    # Case 2: Two non-intersecting squares (valid)
    square_2 = create_square(20, 20, 10, 1)
    assert is_valid_polyset([square_1, square_2]) is True, "Test Case 2 Failed"

    # Case 3: One square inside another (alternating signs, valid)
    inner_square = create_square(0, 0, 5, -1)
    assert is_valid_polyset([square_1, inner_square]) is True, "Test Case 3 Failed"

    # Case 6: Two separate nested sequences inside a bounding box (valid)
    outer_box = create_square(0, 0, 50, 1)
    nested_1 = create_square(-10, -10, 10, -1)
    nested_2 = create_square(10, 10, 10, -1)
    assert is_valid_polyset([outer_box, nested_1, nested_2]) is True, "Test Case 6 Failed"

def test_invalid_polysets():
    """Test cases for invalid PolySets."""

    square_1 = create_square(0, 0, 10, 1)

    # Case 4: Incorrect alternation of signs (invalid)
    wrong_inner_square = create_square(0, 0, 5, 1)  # Should be -1 inside
    assert is_valid_polyset([square_1, wrong_inner_square]) is False, "Test Case 4 Failed"

    # Case 5: Overlapping squares (invalid)
    overlapping_square = create_square(5, 5, 10, 1)
    assert is_valid_polyset([square_1, overlapping_square]) is False, "Test Case 5 Failed"

    # Case 7: Two nested sequences, but one breaks alternation (invalid)
    outer_box = create_square(0, 0, 50, 1)
    nested_1 = create_square(-10, -10, 10, -1)
    wrong_nested = create_square(10, 10, 10, 1)  # Should be -1 inside
    assert is_valid_polyset([outer_box, nested_1, wrong_nested]) is False, "Test Case 7 Failed"

def test_polyset_to_csg():
    # Test Case 1: Single square
    square = create_square(0, 0, 10)
    csg = polyset_to_csg([square])
    assert isinstance(csg, gls.PolyLine2D), "Test Case 1 Failed"

    # Test Case 2: Square with a hole
    outer_square = create_square(0, 0, 10, mode=1)
    inner_square = create_square(0, 0, 5, mode=-1)
    csg = polyset_to_csg([outer_square, inner_square])
    assert isinstance(csg, gls.Difference), "Test Case 2 Failed"

    # Test Case 3: Two disjoint squares
    square_1 = create_square(-20, 0, 10)
    square_2 = create_square(20, 0, 10)
    csg = polyset_to_csg([square_1, square_2])
    assert isinstance(csg, gls.Union), "Test Case 3 Failed"

    # Test Case 4: Nested squares (alternating signs)
    outer_square = create_square(0, 0, 20, mode=1)
    middle_square = create_square(0, 0, 15, mode=-1)
    inner_square = create_square(0, 0, 5, mode=1)
    csg = polyset_to_csg([outer_square, middle_square, inner_square])
    assert isinstance(csg, gls.Difference), "Test Case 4 Failed"

    # Test Case 5: Negative top-level loop
    negative_outer = create_square(0, 0, 20, mode=-1)
    hole = create_square(0, 0, 10, mode=1)
    csg = polyset_to_csg([negative_outer, hole])
    assert isinstance(csg, gls.Complement), "Test Case 5 Failed"

    # Test Case 6: Complex nesting
    root = create_square(0, 0, 70, mode=1)
    hole_1 = create_square(-15, -15, 10, mode=-1)
    island_1 = create_square(-15, -15, 5, mode=1)
    hole_2 = create_square(20, 20, 10, mode=-1)
    csg = polyset_to_csg([root, hole_1, island_1, hole_2])
    assert isinstance(csg, gls.Difference), "Test Case 6 Failed"

    print("All test cases passed!")
# ------------------- FINAL VERIFICATION -----------------



if __name__ == "__main__":
    test_valid_polysets()
    test_invalid_polysets()
    test_polyset_to_csg()
    print("All polyset tests passed successfully!")
