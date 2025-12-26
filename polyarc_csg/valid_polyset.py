"""
Parser for converting CSG expressions to valid PolySets.

Handles DNF/CNF transformations and iterative refinement to produce
valid PolySets from arbitrary CSG expressions.
"""
from typing import Any, List, Optional, Tuple

from geolipi.symbolic.base import GLFunction
from geolipi.torch_compute.sketcher import Sketcher
import geolipi.symbolic as gls
from geolipi.torch_compute.batch_compile import expr_to_dnf, expr_to_cnf
import polyarc_rs as prs

from .polyset import (
    is_valid_polyset, 
    polyset_to_csg, 
    csg_to_polyset, 
    upscale_polyexpr, 
    downscale_polyset, 
    clean_polyset
)
from .polyarc import union_multiple, intersection_multiple, difference_multiple
from .generic_to_polyset import expr_to_polyarc_expr, resolve_difference

__all__ = [
    "expr_to_valid_polyset_expr",
    "resolve_intersection",
    "resolve_unions",
]



def _csg_to_polyterms(expr_list: List[GLFunction]) -> List[prs.PolyArc]:
    """Convert a list of CSG expressions to PolyArc objects."""
    return [
        prs.PolyArc(
            polyarc=tuple(tuple(x) for x in arg.args[0]), 
            is_closed=True, 
            mode=1
        ) 
        for arg in expr_list
    ]


def _get_pos_and_neg(args: Tuple[GLFunction, ...]) -> Tuple[List[prs.PolyArc], List[prs.PolyArc]]:
    """Separate arguments into positive and negative (complemented) terms."""
    positives = [arg for arg in args if not isinstance(arg, gls.Complement)]
    negatives = [arg.args[0] for arg in args if isinstance(arg, gls.Complement)]
    return _csg_to_polyterms(positives), _csg_to_polyterms(negatives)


def _flip_sign(polyset: List[prs.PolyArc]) -> List[prs.PolyArc]:
    """Flips the sign of all primitives in the PolySet."""
    return [
        prs.PolyArc(poly.polyarc, poly.is_closed, -poly.mode) 
        for poly in polyset
    ]


def resolve_intersection(dnf_expression: GLFunction) -> GLFunction:
    """
    Resolves intersection terms in a DNF (Disjunctive Normal Form) expression.
    
    Processes each intersection term by performing the actual geometric
    intersection and converts the result back to a CSG expression.
    
    Args:
        dnf_expression: The DNF form of the expression (Union of Intersections).
    
    Returns:
        Simplified expression with resolved intersections.
    """
    if isinstance(dnf_expression, gls.Union):
        terms = list(dnf_expression.args)
    else:
        terms = [dnf_expression]

    resolved_terms = []
    intersection_terms = []
    root_level_terms = []
    
    for term in terms:
        if isinstance(term, gls.Intersection):
            positives, negatives = _get_pos_and_neg(term.args)
            
            positives_resolved = intersection_multiple(positives) if positives else None
            negatives_resolved = union_multiple(negatives) if negatives else None

            if negatives_resolved and positives_resolved:
                resolved = difference_multiple(positives_resolved, negatives_resolved)
            elif positives_resolved:
                resolved = positives_resolved
            elif negatives_resolved:
                resolved = _flip_sign(negatives_resolved)
            else:
                resolved = None

            if resolved:
                resolved_expr = polyset_to_csg(resolved)
                if not isinstance(resolved_expr, gls.NullExpression2D):
                    intersection_terms.append(resolved_expr)
        else:
            root_level_terms.append(term)

    if root_level_terms:
        positives, negatives = _get_pos_and_neg(root_level_terms)
        
        positives_resolved = union_multiple(positives) if positives else None
        negatives_resolved = intersection_multiple(negatives) if negatives else None

        if negatives_resolved and positives_resolved:
            resolved = _flip_sign(difference_multiple(negatives_resolved, positives_resolved))
        elif positives_resolved:
            resolved = positives_resolved
        elif negatives_resolved:
            resolved = _flip_sign(negatives_resolved)
        else:
            resolved = None

        if resolved:
            resolved_expr = polyset_to_csg(resolved)
            if not isinstance(resolved_expr, gls.NullExpression2D):
                resolved_terms.append(resolved_expr)

    resolved_terms.extend(intersection_terms)
    
    if len(resolved_terms) > 1:
        return gls.Union(*resolved_terms)
    elif len(resolved_terms) == 1:
        return resolved_terms[0]
    else:
        return gls.NullExpression2D()


def resolve_unions(cnf_expression: GLFunction) -> GLFunction:
    """
    Resolves union terms in a CNF (Conjunctive Normal Form) expression.
    
    Processes each union term by performing the actual geometric
    union and converts the result back to a CSG expression.
    
    Args:
        cnf_expression: The CNF form of the expression (Intersection of Unions).
    
    Returns:
        Simplified expression with resolved unions.
    """
    if isinstance(cnf_expression, gls.Intersection):
        terms = list(cnf_expression.args)
    else:
        terms = [cnf_expression]

    resolved_terms = []
    union_terms = []
    root_level_terms = []
    
    for term in terms:
        if isinstance(term, gls.Union):
            positives, negatives = _get_pos_and_neg(term.args)

            positives_resolved = union_multiple(positives) if positives else None
            negatives_resolved = intersection_multiple(negatives) if negatives else None
            ## FIX - If negative_resolved is "everything", then Pos - Everthing is Null
            if negatives and not negatives_resolved:
                positives_resolved = []
            if negatives_resolved and positives_resolved:
                resolved = _flip_sign(difference_multiple(negatives_resolved, positives_resolved))
            elif positives_resolved:
                resolved = positives_resolved
            elif negatives_resolved:
                resolved = _flip_sign(negatives_resolved)
            else:
                resolved = None

            if resolved:
                resolved_expr = polyset_to_csg(resolved)
                if not isinstance(resolved_expr, gls.NullExpression2D):
                    union_terms.append(resolved_expr)
        else:
            root_level_terms.append(term)

    if root_level_terms:
        positives, negatives = _get_pos_and_neg(root_level_terms)
        
        positives_resolved = intersection_multiple(positives) if positives else None
        negatives_resolved = union_multiple(negatives) if negatives else None

        ## FIX - If negative_resolved is "everything", then Pos - Everthing is Null
        if negatives and not negatives_resolved:
            positives_resolved = []
        if negatives_resolved and positives_resolved:
            resolved = difference_multiple(positives_resolved, negatives_resolved)
        elif positives_resolved:
            resolved = positives_resolved
        elif negatives_resolved:
            resolved = _flip_sign(negatives_resolved)
        else:
            resolved = None

        if resolved:
            resolved_expr = polyset_to_csg(resolved)
            if not isinstance(resolved_expr, gls.NullExpression2D):
                resolved_terms.append(resolved_expr)

    resolved_terms.extend(union_terms)
    
    if len(resolved_terms) > 1:
        return gls.Intersection(*resolved_terms)
    elif len(resolved_terms) == 1:
        return resolved_terms[0]
    else:
        return gls.NullExpression2D()


def expr_to_valid_polyset_expr(
    expression: GLFunction, 
    sketcher: Sketcher, 
    *args: Any, 
    **kwargs: Any
) -> GLFunction:
    """
    Converts an arbitrary CSG expression into a valid PolySet CSG expression.
    
    This function handles transforms, nested operations, and iteratively
    refines the expression until it produces a valid PolySet.
    
    Args:
        expression: The input CSG expression.
        sketcher: The sketcher object for affine transformations.
        *args, **kwargs: Additional arguments passed to transform resolution.
    
    Returns:
        A CSG expression representing a valid PolySet.
        
    Raises:
        RuntimeError: If unable to resolve to a valid PolySet after max iterations.
    """
    max_iterations = 20
    iteration_count = 0
    
    # Resolve transforms and convert to polyarc expressions
    expression = expr_to_polyarc_expr(expression, sketcher, *args, **kwargs)
    expression = upscale_polyexpr(expression)
    expression = resolve_difference(expression)

    polyset = csg_to_polyset(expression)
    polyset = clean_polyset(polyset)

    # Check if already valid
    if is_valid_polyset(polyset):
        polyset = downscale_polyset(polyset)
        return polyset_to_csg(polyset)

    # Iterative refinement using DNF/CNF transformations
    while not is_valid_polyset(polyset) and iteration_count < max_iterations:
        # DNF pass
        expression = resolve_difference(expression)
        dnf_expr = expr_to_dnf(expression)
        expression = resolve_intersection(dnf_expr)

        if isinstance(expression, gls.NullExpression2D):
            return expression
        
        # CNF pass
        expression = resolve_difference(expression)
        cnf_expr = expr_to_cnf(expression)
        expression = resolve_unions(cnf_expr)

        polyset = csg_to_polyset(expression)
        polyset = clean_polyset(polyset)
        
        iteration_count += 1
        
        if isinstance(expression, gls.NullExpression2D):
            return expression

    if iteration_count == max_iterations:
        raise RuntimeError("Failed to resolve to a valid PolySet after maximum attempts.")
    
    polyset = downscale_polyset(polyset)
    return polyset_to_csg(polyset)
