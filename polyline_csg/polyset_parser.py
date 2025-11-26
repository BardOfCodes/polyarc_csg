import networkx as nx
from geolipi.symbolic.symbol_types import COMBINATOR_TYPE, PRIM_TYPE
from geolipi.torch_compute.maps import INVERTED_MAP, NORMAL_MAP
from geolipi.symbolic.base import GLFunction
import woodie.symbolic as ws
import geolipi.symbolic as gls
# from geolipi.torch_compute.compile_expression import expr_to_dnf, expr_to_cnf
# import polyline_rs as prs
from itertools import chain, combinations
from .polyset import is_valid_polyset, polyset_to_csg, csg_to_polyset, upscale_polyexpr, downscale_polyset, clean_polyset
from .polyline import union_multiple, intersection_multiple, difference_multiple
from .expression_parser import resolve_to_transform_free_polyline_expr, resolve_difference

def csg_to_polyterms(expr_list):
    poly_list =  [prs.PolyStruct(polyline=tuple([tuple(x) for x in arg.args[0]]), is_closed=True, mode=1) for arg in expr_list]
    return poly_list

def get_pos_and_neg(args):
    positives = [arg for arg in args if not isinstance(arg, gls.Complement)]
    negatives = [arg.args[0] for arg in args if isinstance(arg, gls.Complement)]
    # convert to a list of PolyStructs
    positives = csg_to_polyterms(positives)
    negatives = csg_to_polyterms(negatives)
    return positives, negatives
    
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
            positives, negatives = get_pos_and_neg(term.args)
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
                elif negatives_resolved: 
                    # No positives. 
                    resolved = negatives_resolved
                    resolved = flip_sign(resolved)
                else:
                    resolved = None
            # Now this is a "valid" polyset. 
            # convert to csg expression
            # Consider what should be done here. 
            # -> First only deal with polygons.
            if resolved:
                resolved_expr = polyset_to_csg(resolved)
                # The idea is that if 
                if not isinstance(resolved_expr, gls.NullExpression2D):
                    intersection_terms.append(resolved_expr)
        else:
            # Root level primitives -> must simply be combined. 
            root_level_terms.append(term)
            # raise ValueError(f"Unexpected term type: {type(term)}")
    if root_level_terms:
        positives, negatives = get_pos_and_neg(root_level_terms)
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
            elif negatives_resolved: 
                # No positives. 
                resolved = negatives_resolved
                resolved = flip_sign(resolved)
            else:
                resolved = None
        if resolved:
            resolved_expr = polyset_to_csg(resolved)
            
            if not isinstance(resolved_expr, gls.NullExpression2D):
                resolved_terms.append(resolved_expr)

    resolved_terms.extend(intersection_terms) 
    # Can the union be done here itself to resolve the situation?
    if len(resolved_terms) > 1:
        new_expr = gls.Union(*resolved_terms)
    elif len(resolved_terms) == 1:
        new_expr = resolved_terms[0]
    else:
        new_expr = gls.NullExpression2D()
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
            positives, negatives = get_pos_and_neg(term.args)

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
                elif negatives_resolved: 
                    # No positives. 
                    resolved = negatives_resolved
                    resolved = flip_sign(resolved)
                else:
                    resolved = None
            # Now this is a "valid" polyset. 
            # convert to csg expression
            # Consider what should be done here. 
            # -> First only deal with polygons.
            if resolved:
                resolved_expr = polyset_to_csg(resolved)
                # The idea is that if 
                if not isinstance(resolved_expr, gls.NullExpression2D):
                    intersection_terms.append(resolved_expr)
        else:
            # Root level primitives -> must simply be combined. 
            root_level_terms.append(term)
            # raise ValueError(f"Unexpected term type: {type(term)}")
    if root_level_terms:
        positives, negatives = get_pos_and_neg(root_level_terms)
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
            elif negatives_resolved: 
                # No positives. 
                resolved = negatives_resolved
                resolved = flip_sign(resolved)
            else:
                resolved = None
        if resolved:
            resolved_expr = polyset_to_csg(resolved)
            if not isinstance(resolved_expr, gls.NullExpression2D):
                resolved_terms.append(resolved_expr)

    # check if there are some that can be resolved here directly. 
    # TODO: How to combine intersection terms along with the resolved term?
    resolved_terms.extend(intersection_terms) 
        # Can the union be done here itself to resolve the situation?
    if len(resolved_terms) > 1:
        new_expr = gls.Intersection(*resolved_terms)
    elif len(resolved_terms) == 1:
        new_expr = resolved_terms[0]
    else:
        new_expr = gls.NullExpression2D()
    return new_expr


def parse_csg_to_valid_polyset_csg(expression, sketcher, *args, **kwargs):
    """
    Converts an arbitrary CSG expression into a valid PolySet.
    Input Expression should be in Tuple Form
    
    Args:
        expression (GLFunction): The input CSG expression.
    
    Returns:
        List[prs.PolyStruct]: A valid PolySet representation.
    """
    # Initial extraction of primitives
    max_iterations = 20
    iteration_count = 0
    # PolyExpr
    expression = resolve_to_transform_free_polyline_expr(expression, sketcher, *args, **kwargs)
    expression = upscale_polyexpr(expression)
    expression = resolve_difference(expression)
    

    polyset = csg_to_polyset(expression)
    # Remove nulls. 
    polyset = clean_polyset(polyset)
    

    # Validate the initial polyset
    if is_valid_polyset(polyset):
        polyset = downscale_polyset(polyset)
        polyset_csg = polyset_to_csg(polyset)
        return polyset_csg

    # Iterative refinement
    while not is_valid_polyset(polyset) and iteration_count < max_iterations:
        # Convert to DNF
        expression = resolve_difference(expression)
        dnf_expr = expr_to_dnf(expression)
        expression = resolve_intersection(dnf_expr)

        if isinstance(expression, gls.NullExpression2D):
            return expression
        
        expression = resolve_difference(expression)
        cnf_expr = expr_to_cnf(expression)
        expression = resolve_unions(cnf_expr)

        # Extract updated primitives
        polyset = csg_to_polyset(expression)
        polyset = clean_polyset(polyset)
        # primitives_with_signs = extract_primitives_with_signs(expression)
        # polyset = [prs.PolyStruct(polyline=list(tuple(x) for x in prim.args[0]), is_closed=True, mode=sign)
        #            for prim, sign in primitives_with_signs]
        
        iteration_count += 1
        if isinstance(expression, gls.NullExpression2D):
            return expression
    

    if iteration_count == max_iterations:
        raise RuntimeError("Failed to resolve to a valid PolySet after maximum attempts.")
    
    polyset = downscale_polyset(polyset)
    polyset_csg = polyset_to_csg(polyset)
    return polyset_csg

