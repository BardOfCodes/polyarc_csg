"""
PolySet operations - validation, CSG conversion, and enclosure tree construction.

A PolySet is a list of PolyArc objects representing a 2D region with potential holes.
"""
import logging
from typing import List, Optional, Tuple

import networkx as nx
import polyarc_rs as prs
import geolipi.symbolic as gls
from geolipi.symbolic.base import GLFunction

__all__ = [
    "UPSCALING_FACTOR",
    "set_upscaling_factor",
    "get_upscaling_factor",
    "construct_enclosure_tree",
    "construct_enclosure_sequences",
    "is_valid_polyset",
    "polyset_to_csg",
    "csg_to_polyset",
    "clean_polyset",
    "upscale_polyexpr",
    "downscale_polyset",
    "extract_primitives_with_signs",
]

logger = logging.getLogger(__name__)

# Default upscaling factor for higher numerical precision
# Can be overridden by setting polyarc_csg.polyset.UPSCALING_FACTOR
UPSCALING_FACTOR = 1000


def set_upscaling_factor(factor: int) -> None:
    """
    Set the global upscaling factor for numerical precision.
    
    Args:
        factor: The new upscaling factor (default is 1000).
    """
    global UPSCALING_FACTOR
    UPSCALING_FACTOR = factor


def get_upscaling_factor() -> int:
    """Get the current upscaling factor."""
    return UPSCALING_FACTOR


def construct_enclosure_tree(polyset: List[prs.PolyArc]) -> nx.DiGraph:
    """
    Constructs a directed graph representing immediate enclosure relationships.
    
    Args:
        polyset: List of PolyArcs.
        
    Returns:
        A transitively reduced directed graph where edges represent
        immediate enclosure (parent encloses child).
    """
    G = nx.DiGraph()
    
    for poly in polyset:
        G.add_node(poly)
    
    for i, poly_a in enumerate(polyset):
        for j, poly_b in enumerate(polyset):
            if i != j and prs.is_1_inside_2(poly_b, poly_a):
                G.add_edge(poly_a, poly_b)
    
    return nx.transitive_reduction(G)


def construct_enclosure_sequences(polyset: List[prs.PolyArc]) -> List[List[prs.PolyArc]]:
    """
    Constructs enclosure sequences from a PolySet by identifying enclosure relationships.

    Args:
        polyset: The input list of PolyArcs.

    Returns:
        A list of enclosure sequences (from outermost to innermost).
    """
    TR = construct_enclosure_tree(polyset)
    roots = [n for n in TR.nodes if TR.in_degree(n) == 0]
    
    sequences = []

    def extract_sequences(node, path):
        path.append(node)
        if TR.out_degree(node) == 0:
            sequences.append(path.copy())
        else:
            for child in TR.successors(node):
                extract_sequences(child, path)
        path.pop()

    for root in roots:
        extract_sequences(root, [])

    return sequences


def is_valid_polyset(polyset: List[prs.PolyArc]) -> bool:
    """
    Checks if a given PolySet is valid.
    
    A valid PolySet must:
    1. Have non-intersecting polyarcs.
    2. Have sequences of alternating positive and negative mode values.
    
    Args:
        polyset: The input list of PolyArcs.
    
    Returns:
        True if valid, False otherwise.
    """
    if not polyset:
        return True
    
    # Check for intersections
    for i, poly_a in enumerate(polyset):
        for j, poly_b in enumerate(polyset):
            if i != j:
                if prs.is_intersected(poly_a, poly_b):
                    return False
                if prs.is_overlapping(poly_a, poly_b):
                    return False

    # Construct enclosure sequences and check alternation
    sequences = construct_enclosure_sequences(polyset)
    for seq in sequences:
        for k in range(len(seq) - 1):
            if seq[k].mode == seq[k + 1].mode:
                return False
    
    return True


def polyset_to_csg(polyset: List[prs.PolyArc]) -> GLFunction:
    """
    Converts a PolySet (list of PolyArcs) into a CSG expression.
    
    Args:
        polyset: The input PolySet.
    
    Returns:
        The corresponding CSG expression.
    """
    if not polyset:
        return gls.NullExpression2D()
    
    TR = construct_enclosure_tree(polyset)
    roots = [n for n in TR.nodes if TR.in_degree(n) == 0]

    def construct_csg(node):
        children = list(TR.successors(node))
        
        if not children:
            return gls.PolyArc2D(node.polyarc)

        child_exprs = [construct_csg(child) for child in children]

        if len(child_exprs) == 1:
            child_expr = child_exprs[0]
        else:
            child_expr = gls.Union(*child_exprs)
        
        return gls.Difference(gls.PolyArc2D(node.polyarc), child_expr)

    root_exprs = [construct_csg(root) for root in roots]
    root_exprs = [
        gls.Complement(expr) if roots[ind].mode == -1 else expr 
        for ind, expr in enumerate(root_exprs)
    ]

    if len(root_exprs) > 1:
        return gls.Union(*root_exprs)
    elif len(root_exprs) == 1:
        return root_exprs[0]
    else:
        return gls.NullExpression2D()


def extract_primitives_with_signs(
    expression: GLFunction, 
    current_sign: int = 1
) -> List[Tuple[gls.PolyArc2D, int]]:
    """
    Recursively extracts primitives from the expression tree, tracking their signs.
    
    Args:
        expression: The CSG expression.
        current_sign: The current sign (+1 or -1) based on parent operations.
    
    Returns:
        A list of (primitive, sign) tuples.
    """
    primitives = []

    if isinstance(expression, gls.PolyArc2D):
        primitives.append((expression, current_sign))
    elif isinstance(expression, gls.Complement):
        primitives.extend(extract_primitives_with_signs(expression.args[0], -current_sign))
    elif isinstance(expression, gls.Difference):
        primitives.extend(extract_primitives_with_signs(expression.args[0], current_sign))
        primitives.extend(extract_primitives_with_signs(expression.args[1], -current_sign))
    elif isinstance(expression, (gls.Union, gls.Intersection)):
        for child in expression.args:
            primitives.extend(extract_primitives_with_signs(child, current_sign))

    return primitives


def _tuple_to_tuple(polyarc) -> tuple:
    """Convert a polyarc argument to a tuple of (x, y, bulge) tuples."""
    return tuple((x[0], x[1], x[2]) for x in polyarc)


def csg_to_polyset(expression: GLFunction) -> List[prs.PolyArc]:
    """
    Converts a CSG expression to a PolySet.
    
    Args:
        expression: The CSG expression.
        
    Returns:
        A list of PolyArc objects.
    """
    primitives = extract_primitives_with_signs(expression)
    polyset = []
    for prim, sign in primitives:
        updated_points = _tuple_to_tuple(prim.args[0])
        polyset.append(prs.PolyArc(polyarc=updated_points, is_closed=True, mode=sign))
    return polyset


def clean_polyset(polyset: List[prs.PolyArc]) -> List[prs.PolyArc]:
    """
    Cleans a PolySet by removing redundant vertices and empty shapes.
    
    Args:
        polyset: The input PolySet.
        
    Returns:
        A cleaned PolySet.
    """
    cleaned_polyset = []
    for poly in polyset:
        try:
            cleaned_poly = prs.remove_redundant_vertices(poly, epsilon=1e-5)
            if cleaned_poly is None:
                cleaned_poly = poly
            
            length = prs.compute_path_length(cleaned_poly)
            area = prs.compute_area(cleaned_poly)

            if abs(area) > 0.0 and abs(length) > 0.0:
                cleaned_polyset.append(cleaned_poly)
        except Exception as e:
            logger.warning(f"Failed to clean polyarc: {e}")
    
    return cleaned_polyset


def upscale_polyexpr(polyexpr: GLFunction, factor: int = UPSCALING_FACTOR) -> GLFunction:
    """
    Upscales coordinates in a poly expression for better numerical precision.
    
    Args:
        polyexpr: The expression to upscale.
        factor: The scaling factor.
        
    Returns:
        The upscaled expression.
    """
    if isinstance(polyexpr, gls.PolyArc2D):
        points = polyexpr.args[0]
        upscaled_points = tuple((x[0] * factor, x[1] * factor, x[2]) for x in points)
        return gls.PolyArc2D(upscaled_points)
    else:
        new_args = [upscale_polyexpr(arg, factor) for arg in polyexpr.args]
        return type(polyexpr)(*new_args)


def downscale_polyset(
    polyset: List[prs.PolyArc], 
    factor: int = UPSCALING_FACTOR
) -> List[prs.PolyArc]:
    """
    Downscale a PolySet by a given factor.
    
    Args:
        polyset: The input PolySet.
        factor: The downscaling factor.
    
    Returns:
        The downscaled PolySet.
    """
    def downscale_polyarc(polyarc, factor):
        return tuple((x[0] / factor, x[1] / factor, x[2]) for x in polyarc)

    return [
        prs.PolyArc(downscale_polyarc(poly.polyarc, factor), poly.is_closed, poly.mode) 
        for poly in polyset
    ]
