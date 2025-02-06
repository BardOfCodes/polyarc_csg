import networkx as nx
from geolipi.symbolic.types import COMBINATOR_TYPE, PRIM_TYPE
from geolipi.torch_compute.utils import INVERTED_MAP, NORMAL_MAP
from geolipi.symbolic.base_symbolic import GLFunction
import woodie.symbolic as ws
import geolipi.symbolic as gls
from geolipi.torch_compute.compile_expression import expr_to_dnf, expr_to_cnf
import polyline_rs as prs
from . import polyline as pl
from itertools import chain, combinations
from .polyset import is_valid_polyset, polyset_to_csg, csg_to_polyset, upscale_polyexpr, downscale_polyset
from .polyline import union_multiple, intersection_multiple, difference_multiple
from .expression_parser import resolve_transfroms_and_params, resolve_difference


def flip_sign(polyset):
    """
    Flips the sign of all primitives in the PolySet.
    
    Args:
        polyset (List[prs.PolyStruct]): The input PolySet.
    
    Returns:
        List[prs.PolyStruct]: The PolySet with flipped signs.
    """
    for ind, poly in enumerate(polyset):
        polyset[ind] = prs.PolyStruct(poly.polyline, poly.is_closed, -poly.mode)
    return polyset

def resolve_intersection(dnf_expression):
    """
    Resolves intersections in the DNF expression.
    
    Args:
        dnf_expression (GLFunction): The DNF form of the expression.
    
    Returns:
        GLFunction: Simplified expression with resolved intersections.
    """
    if isinstance(dnf_expression, gls.Union):
        terms = dnf_expression.args
    else:
        terms = [dnf_expression]

    resolved_terms = []
    intersection_terms = []
    root_level_terms = []
    for term in terms:
        if isinstance(term, gls.Intersection):
            positives = [arg for arg in term.args if not isinstance(arg, gls.Complement)]
            negatives = [arg.args[0] for arg in term.args if isinstance(arg, gls.Complement)]
            # convert to a list of PolyStructs
            positives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=1) for arg in positives]
            negatives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=1) for arg in negatives]

            if positives:
                positives_resolved = intersection_multiple(positives)
            else:
                positives_resolved = None

            if negatives:
                negatives_resolved = union_multiple(negatives)
            else:
                negatives_resolved = None

            if negatives_resolved and positives_resolved:
                # Here we need to do this considering that there could be multiple positive A and multiple negative B. 
                resolved = difference_multiple(positives_resolved, negatives_resolved)
            else:
                if positives_resolved:
                    # No negatives. 
                    resolved = positives_resolved
                else: 
                    # No positives. 
                    resolved = negatives_resolved
                    resolved = flip_sign(resolved)
            # Now this is a "valid" polyset. 
            # convert to csg expression
            # Consider what should be done here. 
            # -> First only deal with polygons.
            resolved_expr = polyset_to_csg(resolved)
            # The idea is that if 
            intersection_terms.append(resolved_expr)
        else:
            # Root level primitives -> must simply be combined. 
            root_level_terms.append(term)
            # raise ValueError(f"Unexpected term type: {type(term)}")
    if root_level_terms:
        positives = [arg for arg in root_level_terms if not isinstance(arg, gls.Complement)]
        negatives = [arg for arg in root_level_terms if isinstance(arg, gls.Complement)]
        positives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=1) for arg in positives]
        negatives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=1) for arg in negatives]

        if positives:
            positives_resolved = union_multiple(positives)
        else:
            positives_resolved = None

        if negatives:
            negatives_resolved = intersection_multiple(negatives)
        else:
            negatives_resolved = None

        if negatives_resolved and positives_resolved:
            # Here we need to do this considering that there could be multiple positive A and multiple negative B. 
            resolved = difference_multiple(negatives_resolved, positives_resolved)
            resolved = flip_sign(resolved)
        else:
            if positives_resolved:
                # No negatives. 
                resolved = positives_resolved
            else: 
                # No positives. 
                resolved = negatives_resolved
                resolved = flip_sign(resolved)

        
        resolved_expr = polyset_to_csg(resolved)
        resolved_terms.append(resolved_expr)

    if intersection_terms:
        # check if there are some that can be resolved here directly. 
        # TODO: How to combine intersection terms along with the resolved term?
        resolved_terms.extend(intersection_terms) 
        # Can the union be done here itself to resolve the situation?
    new_expr = gls.Union(*resolved_terms) if len(resolved_terms) > 1 else resolved_terms[0]
    return new_expr


def resolve_unions(cnf_expression: GLFunction):
    """
    Resolves intersections in the CNF expression.
    
    Args:
        cnf_expression (GLFunction): The CNF form of the expression.
    
    Returns:
        GLFunction: Simplified expression with resolved intersections.
    """
    if isinstance(cnf_expression, gls.Intersection):
        terms = cnf_expression.args
    else:
        terms = [cnf_expression]

    resolved_terms = []
    intersection_terms = []
    root_level_terms = []
    for term in terms:
        if isinstance(term, gls.Union):
            positives = [arg for arg in term.args if not isinstance(arg, gls.Complement)]
            negatives = [arg.args[0] for arg in term.args if isinstance(arg, gls.Complement)]
            # convert to a list of PolyStructs
            positives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=1) for arg in positives]
            negatives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=1) for arg in negatives]

            if positives:
                positives_resolved = union_multiple(positives)
            else:
                positives_resolved = None

            if negatives:
                negatives_resolved = intersection_multiple(negatives)
            else:
                negatives_resolved = None

            if negatives_resolved and positives_resolved:
                # Here we need to do this considering that there could be multiple positive A and multiple negative B. 
                resolved = difference_multiple(negatives_resolved, positives_resolved)
                resolved = flip_sign(resolved)
            else:
                if positives_resolved:
                    # No negatives. 
                    resolved = positives_resolved
                else: 
                    # No positives. 
                    resolved = negatives_resolved
                    resolved = flip_sign(resolved)
            # Now this is a "valid" polyset. 
            # convert to csg expression
            # Consider what should be done here. 
            # -> First only deal with polygons.
            resolved_expr = polyset_to_csg(resolved)
            # The idea is that if 
            intersection_terms.append(resolved_expr)
        else:
            # Root level primitives -> must simply be combined. 
            root_level_terms.append(term)
            # raise ValueError(f"Unexpected term type: {type(term)}")
    if root_level_terms:
        positives = [arg for arg in root_level_terms if not isinstance(arg, gls.Complement)]
        negatives = [arg for arg in root_level_terms if isinstance(arg, gls.Complement)]
        positives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=1) for arg in positives]
        negatives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=1) for arg in negatives]

        if positives:
            positives_resolved = intersection_multiple(positives)
        else:
            positives_resolved = None

        if negatives:
            negatives_resolved = union_multiple(negatives)
        else:
            negatives_resolved = None

        if negatives_resolved and positives_resolved:
            # Here we need to do this considering that there could be multiple positive A and multiple negative B. 
            resolved = difference_multiple(positives_resolved, negatives_resolved)
        else:
            if positives_resolved:
                # No negatives. 
                resolved = positives_resolved
            else: 
                # No positives. 
                resolved = negatives_resolved
                resolved = flip_sign(resolved)
        
        resolved_expr = polyset_to_csg(resolved)
        resolved_terms.append(resolved_expr)

    if intersection_terms:
        # check if there are some that can be resolved here directly. 
        # TODO: How to combine intersection terms along with the resolved term?
        resolved_terms.extend(intersection_terms) 
        # Can the union be done here itself to resolve the situation?
    new_expr = gls.Intersection(*resolved_terms) if len(resolved_terms) > 1 else resolved_terms[0]
    return new_expr


def parse_csg_to_valid_polyset_csg(expression, *args, **kwargs):
    """
    Converts an arbitrary CSG expression into a valid PolySet.
    
    Args:
        expression (GLFunction): The input CSG expression.
    
    Returns:
        List[prs.PolyStruct]: A valid PolySet representation.
    """
    # Initial extraction of primitives
    max_iterations = 20
    iteration_count = 0
    # PolyExpr
    expression = resolve_transfroms_and_params(expression, *args, **kwargs)
    expression = upscale_polyexpr(expression)

    expression = resolve_difference(expression)

    polyset = csg_to_polyset(expression)

    # Validate the initial polyset
    if is_valid_polyset(polyset):
        polyset = downscale_polyset(polyset)
        polyset_csg = polyset_to_csg(polyset)
        return polyset_csg

    # Iterative refinement
    while not is_valid_polyset(polyset) and iteration_count < max_iterations:
        # Convert to DNF
        dnf_expr = expr_to_dnf(expression)
        expression = resolve_intersection(dnf_expr)

        cnf_expr = expr_to_cnf(expression)
        expression = resolve_unions(cnf_expr)

        # Extract updated primitives
        polyset = csg_to_polyset(expression)
        # primitives_with_signs = extract_primitives_with_signs(expression)
        # polyset = [prs.PolyStruct(polyline=list(tuple(x) for x in prim.args[0]), is_closed=True, mode=sign)
        #            for prim, sign in primitives_with_signs]
        
        iteration_count += 1
    

    if iteration_count == max_iterations:
        raise RuntimeError("Failed to resolve to a valid PolySet after maximum attempts.")
    
    polyset = downscale_polyset(polyset)
    polyset_csg = polyset_to_csg(polyset)
    return polyset_csg

