import polyline_rs as prs
import geolipi.symbolic as gls
import woodie.symbolic as ws
from polyline_csg.parser import parse_csg_to_valid_polyset
from polyline_csg.polyset import is_valid_polyset

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

def test_single_primitive():
    """Test if a single primitive ws.PolyLine2D is correctly converted to a PolySet."""
    expr = ws.PolyLine2D([
        (0, 0, 0.0), (10, 0, 0.0), (10, 10, 0.0), (0, 10, 0.0)
    ])
    polyset = parse_csg_to_valid_polyset(expr)

    assert len(polyset) == 1, "Test Failed: Expected a single polyline"
    assert polyset[0].mode == 1, "Test Failed: Expected mode 1 for single primitive"
    assert is_valid_polyset(polyset), "Test Failed: Generated PolySet is invalid"

def test_simple_boolean_operations():
    """Test basic boolean operations with two primitives."""
    
    square_1 = ws.PolyLine2D([(0, 0, 0.0), (10, 0, 0.0), (10, 10, 0.0), (0, 10, 0.0)])
    square_2 = ws.PolyLine2D([(5, 5, 0.0), (15, 5, 0.0), (15, 15, 0.0), (5, 15, 0.0)])

    # Test Union
    union_expr = gls.Union(square_1, square_2)
    union_polyset = parse_csg_to_valid_polyset(union_expr)
    assert is_valid_polyset(union_polyset), "Union result is not a valid PolySet"

    # Test Intersection
    intersection_expr = gls.Intersection(square_1, square_2)
    intersection_polyset = parse_csg_to_valid_polyset(intersection_expr)
    assert is_valid_polyset(intersection_polyset), "Intersection result is not a valid PolySet"

    # Test Difference
    difference_expr = gls.Difference(square_1, square_2)
    difference_polyset = parse_csg_to_valid_polyset(difference_expr)
    assert is_valid_polyset(difference_polyset), "Difference result is not a valid PolySet"

def test_large_expression():
    """Test a complex expression combining multiple operations."""
    
    outer_box = ws.PolyLine2D([(0, 0, 0.0), (50, 0, 0.0), (50, 50, 0.0), (0, 50, 0.0)])
    cutout_1 = ws.PolyLine2D([(10, 10, 0.0), (20, 10, 0.0), (20, 20, 0.0), (10, 20, 0.0)])
    cutout_2 = ws.PolyLine2D([(30, 30, 0.0), (40, 30, 0.0), (40, 40, 0.0), (30, 40, 0.0)])

    expr = gls.Difference(outer_box, gls.Union(cutout_1, cutout_2))
    polyset = parse_csg_to_valid_polyset(expr)
    
    assert is_valid_polyset(polyset), "Large expression does not generate a valid PolySet"

def test_union_enclosed_polylines():
    """Union of enclosed polylines should result in a single outer polyline."""
    
    outer_square = create_square(0, 0, 50, 1)
    inner_hole = create_square(0, 0, 20, -1)

    union_expr = gls.Union(ws.PolyLine2D(outer_square.polyline), ws.PolyLine2D(inner_hole.polyline))
    result_polyset = parse_csg_to_valid_polyset(union_expr)

    assert len(result_polyset) == 1, "Test Failed: Expected a single resulting polyline"
    assert prs.is_overlapping(result_polyset[0], outer_square), "Test Failed: Result should be the outer polyline"

def test_union_intersecting_holes():
    """Union between two concentric circles each with a hole that intersects."""
    square_1 = create_square(0, 0, 10, 1)
    square_2 = create_square(0, 0, 20, 1)
    square_3 = create_square(0, 0, 30, 1)
    square_4 = create_square(0, 0, 40, 1)
    squares = [square_1, square_2, square_3, square_4]
    exprs = [ws.PolyLine2D(tuple(square.polyline)) for square in squares]
    expr_1, expr_2, expr_3, expr_4 = exprs
    expr = gls.Difference(expr_4, gls.Difference(expr_3, gls.Difference(expr_2, expr_1)))

    result_polyset = parse_csg_to_valid_polyset(expr)

    assert is_valid_polyset(result_polyset), "Test Failed: Union of intersecting holes does not create a valid PolySet"


    outer_circle_1 = ws.PolyLine2D(
        tuple([(0, 0, 0.5), (1, 0, 0.5), (1, 1, 0.5), (0, 1, 0.5)])
        )
    hole_1 = ws.PolyLine2D(
        tuple([
        (.2, .2, 0.5), (.8, .2, 0.5), (.8, .8, 0.5), (.2, .8, 0.5)
    ])
    )
    outer_circle_2 = ws.PolyLine2D(
        tuple([
        (0.5, 0, 0.5), (1.5, 0, 0.5), (1.5, 1, 0.5), (0.5, 1, 0.5)
    ]))
    hole_2 = ws.PolyLine2D(
        tuple([
        (0.7, .2, 0.5), (1.3, .2, 0.5), (1.3, .8, 0.5), (0.7, .8, 0.5)
    ]))


    expr = gls.Union(
        gls.Difference(outer_circle_1, hole_1),
        gls.Difference(outer_circle_2, hole_2)
    )
    result_polyset = parse_csg_to_valid_polyset(expr)

    assert is_valid_polyset(result_polyset), "Test Failed: Union of intersecting holes does not create a valid PolySet"

if __name__ == "__main__":
    # test_single_primitive()
    # test_simple_boolean_operations()
    # test_large_expression()
    # test_union_enclosed_polylines()
    test_union_intersecting_holes()
    print("All parser tests passed successfully!")