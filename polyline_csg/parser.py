import networkx as nx
from geolipi.symbolic.types import COMBINATOR_TYPE, PRIM_TYPE
from geolipi.torch_compute.utils import INVERTED_MAP, NORMAL_MAP
from geolipi.symbolic.base_symbolic import GLFunction
import woodie.symbolic as ws
import geolipi.symbolic as gls
from geolipi.torch_compute.compile_expression import expr_to_dnf, expr_to_cnf
import polyline_rs as prs
from .polyset import is_valid_polyset, polyset_to_csg
from . import polyline as pl
from itertools import chain, combinations
from .polyline import union_multiple, intersection_multiple, difference_multiple

def extract_primitives_with_signs(expression, current_sign=1):
    """
    Recursively extracts primitives (ws.PolyLine2D) from the expression tree, tracking their signs.
    
    Args:
        expression (GLFunction): The CSG expression.
        current_sign (int): The current sign (+1 or -1) based on parent operations.
    
    Returns:
        List[Tuple[ws.PolyLine2D, int]]: A list of primitives with their corresponding signs.
    """
    primitives = []

    if isinstance(expression, ws.PolyLine2D):
        primitives.append((expression, current_sign))
    
    elif isinstance(expression, gls.Complement):
        primitives.extend(extract_primitives_with_signs(expression.args[0], -current_sign))

    elif isinstance(expression, gls.Difference):
        # Left child retains the sign, right child flips the sign
        primitives.extend(extract_primitives_with_signs(expression.args[0], current_sign))
        primitives.extend(extract_primitives_with_signs(expression.args[1], -current_sign))
    
    elif isinstance(expression, (gls.Union, gls.Intersection)):
        for child in expression.args:
            primitives.extend(extract_primitives_with_signs(child, current_sign))

    return primitives

def csg_to_polyset(expression: GLFunction):

    primitives = extract_primitives_with_signs(expression)
    polyset = []
    for prim in primitives:
        polyset.append(prs.PolyStruct(polyline=list(tuple(x) for x in prim[0].args[0]), is_closed=True, mode=prim[1]))
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
            negatives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=-1) for arg in negatives]

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
        positives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=1) for arg in root_level_terms]
        negatives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=-1) for arg in root_level_terms]

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
            negatives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=-1) for arg in negatives]

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
            else:
                if positives_resolved:
                    # No negatives. 
                    resolved = positives_resolved
                else: 
                    # No positives. 
                    resolved = negatives_resolved
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
        positives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=1) for arg in root_level_terms]
        negatives = [prs.PolyStruct(polyline=list(tuple(x) for x in arg.args[0]), is_closed=True, mode=-1) for arg in root_level_terms]

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
        else:
            if positives_resolved:
                # No negatives. 
                resolved = positives_resolved
            else: 
                # No positives. 
                resolved = negatives_resolved
        
        resolved_expr = polyset_to_csg(resolved)
        resolved_terms.append(resolved_expr)

    if intersection_terms:
        # check if there are some that can be resolved here directly. 
        # TODO: How to combine intersection terms along with the resolved term?
        resolved_terms.extend(intersection_terms) 
        # Can the union be done here itself to resolve the situation?
    new_expr = gls.Complement(gls.Union(*resolved_terms)) if len(resolved_terms) > 1 else gls.Complement(resolved_terms[0])
    return new_expr

def resolve_difference(expression: GLFunction, ):
    """
    Compiles a GL expression into a format suitable for batch evaluation, gathering transformations,
    primitive parameters, and handling difference and complement operations.

    This function traverses the expression tree to extract and organize necessary data for
    rendering the expression. It accounts for transformations, inversion modes, and primitive
    parameters, and resolves any complex operations like Difference and Complement.

    Parameters:
        expression (GLExpr): The GL expression to be compiled.
        sketcher (Sketcher, optional): The sketcher object used for affine transformations.
            Defaults to None.
        rectify_transform (bool, optional): Flag to determine if transformations should be rectified.
            Defaults to RECTIFY_TRANSFORM.

    Returns:
        tuple: A tuple containing the compiled expression, primitive transformations, primitive
        inversions, and primitive parameters. These are organized to facilitate batch evaluation
        of the expression.
    """

    # inversions
    inversion_mode = False
    inversion_stack = [inversion_mode]

    execution_stack = []
    execution_pointer_index = []
    operator_stack = []
    operator_nargs_stack = []
    operator_params_stack = []

    parser_list = [expression]
    while parser_list:
        cur_expr = parser_list.pop()
        inversion_mode = inversion_stack.pop()
        if isinstance(cur_expr, COMBINATOR_TYPE):
            tree_branches, param_list = [], []
            for arg in cur_expr.args:
                if arg in cur_expr.lookup_table:
                    param_list.append(cur_expr.lookup_table[arg])
                else:
                    tree_branches.append(arg)
            n_args = len(tree_branches)

            if type(cur_expr) == gls.Difference:
                inversion_stack.append(not inversion_mode)
            else:
                inversion_stack.append(inversion_mode)

            inversion_stack.extend([inversion_mode for x in range(n_args - 1)])

            if inversion_mode:
                current_symbol = INVERTED_MAP[type(cur_expr)]
            else:
                current_symbol = NORMAL_MAP[type(cur_expr)]
            operator_stack.append(current_symbol)
            operator_nargs_stack.append(n_args)
            operator_params_stack.append(param_list)
            next_to_parse = tree_branches[::-1]
            parser_list.extend(next_to_parse)
            execution_pointer_index.append(len(execution_stack))
        elif isinstance(cur_expr, PRIM_TYPE):

            if inversion_mode:
                prim_expr = gls.Complement(cur_expr)
            else:
                prim_expr = cur_expr
            execution_stack.append(prim_expr)
        else:
            raise ValueError(f"Unknown expression type {type(cur_expr)}")

        while (
            operator_stack
            and len(execution_stack) - execution_pointer_index[-1]
            >= operator_nargs_stack[-1]
        ):
            n_args = operator_nargs_stack.pop()
            operator = operator_stack.pop()
            _ = execution_pointer_index.pop()
            params = operator_params_stack.pop()
            args = execution_stack[-n_args:]
            new_canvas = operator(*args, *params)
            execution_stack = execution_stack[:-n_args] + [new_canvas]
    expression = execution_stack[0]
    return expression



def expr_to_dnf_full(expression):
    # convert all difference operators to intersection and complement
    expression = resolve_difference(expression)
    dnf_expression = expr_to_dnf(expression)
    return dnf_expression

def expr_to_cnf_full(expression):
    # convert all difference operators to intersection and complement
    expression = resolve_difference(expression)
    cnf_expression = expr_to_cnf(expression)
    return cnf_expression

def parse_csg_to_valid_polyset(expression):
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
    polyset = csg_to_polyset(expression)

    # Validate the initial polyset
    if is_valid_polyset(polyset):
        return polyset

    # Iterative refinement
    while not is_valid_polyset(polyset) and iteration_count < max_iterations:
        # Convert to DNF
        dnf_expr = expr_to_dnf_full(expression)

        # Resolve intersections
        expression = resolve_intersection(dnf_expr)

        cnf_expr = expr_to_cnf_full(expression)
        expression = resolve_unions(cnf_expr)

        # Extract updated primitives
        primitives_with_signs = extract_primitives_with_signs(expression)
        polyset = [prs.PolyStruct(polyline=list(tuple(x) for x in prim.args[0]), is_closed=True, mode=sign)
                   for prim, sign in primitives_with_signs]
        
        iteration_count += 1
    

    if iteration_count == max_iterations:
        raise RuntimeError("Failed to resolve to a valid PolySet after maximum attempts.")
    return polyset
