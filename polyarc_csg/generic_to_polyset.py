"""
Expression parser for transforming CSG expressions into PolyArc2D primitives.

Handles transforms, parametric expressions, and primitive conversion.
"""
from typing import Any, List, Optional, Tuple, Type

import torch as th
import sympy as sp
import geolipi.symbolic as gls
from geolipi.symbolic.base import GLFunction
from geolipi.torch_compute.evaluate_expression import singledispatch, _parse_param_from_expr
from geolipi.symbolic.symbol_types import PRIM_TYPE, COMBINATOR_TYPE, MOD_TYPE
from geolipi.torch_compute.sketcher import Sketcher
from geolipi.torch_compute.maps import INVERTED_MAP, NORMAL_MAP, MODIFIER_MAP
from .prim_map import PRIMITIVE_MAP

__all__ = [
    "resolve_difference",
    "expr_to_polyarc_expr",
]


def expr_to_polyarc_expr(
    expression: GLFunction, 
    sketcher: Sketcher, 
    transform: Optional[th.Tensor] = None
) -> GLFunction:
    """
    Convert a CSG expression to a PolyArc2D expression by resolving transforms.
    
    Args:
        expression: The input CSG expression.
        sketcher: The sketcher object providing coordinate transformations.
        transform: Optional initial affine transform. If None, identity is used.
    
    Returns:
        A CSG expression with all primitives converted to PolyArc2D.
    """
    if transform is None:
        transform = sketcher.get_affine_identity()
    resolved_expr = resolve_expression(expression, sketcher, transform)
    
    return resolved_expr


@singledispatch
def resolve_expression(
    expression: GLFunction, 
    sketcher: Sketcher, 
    transform: Optional[th.Tensor] = None
) -> GLFunction:
    """
    Recursively resolve a CSG expression to PolyArc2D primitives.
    
    Args:
        expression: The CSG expression to resolve.
        sketcher: The sketcher object.
        transform: The accumulated affine transform.
    
    Returns:
        Resolved expression with PolyArc2D primitives.
    
    Raises:
        NotImplementedError: If the expression type is not supported.
    """
    raise NotImplementedError(f"Expression type {type(expression)} not implemented")


@resolve_expression.register
def resolve_expression_combinator(
    expression: COMBINATOR_TYPE, 
    sketcher: Sketcher, 
    transform: Optional[th.Tensor] = None
) -> GLFunction:
    """Resolve combinator expressions (Union, Intersection, Difference, etc.)."""
    tree_branches: List[GLFunction] = []
    param_list: List[Any] = []
    
    for arg in expression.args:
        if arg in expression.lookup_table:
            assert isinstance(arg, sp.Symbol), "Argument must be a symbol"
            param_list.append(expression.lookup_table[arg])
        else:
            tree_branches.append(arg)
    
    resolved_expr_list: List[GLFunction] = []
    for child in tree_branches:
        cur_resolved_expr = resolve_expression(child, sketcher, transform.clone())
        resolved_expr_list.append(cur_resolved_expr)
    
    new_expr = type(expression)(*resolved_expr_list, *param_list)
    return new_expr


@resolve_expression.register
def resolve_expression_mod(
    expression: MOD_TYPE, 
    sketcher: Sketcher, 
    transform: Optional[th.Tensor] = None
) -> GLFunction:
    """Resolve modifier expressions (transforms like Translate, Rotate, etc.)."""
    sub_expr = expression.args[0]
    params = expression.args[1:]
    params = _parse_param_from_expr(expression, params, sketcher)
    new_transform = MODIFIER_MAP[type(expression)](transform, *params)
    new_resolved_expr = resolve_expression(sub_expr, sketcher, new_transform)
    return new_resolved_expr


@resolve_expression.register
def resolve_expression_prim(
    expression: PRIM_TYPE, 
    sketcher: Sketcher, 
    transform: Optional[th.Tensor] = None
) -> GLFunction:
    """Resolve primitive expressions to PolyArc2D."""
    params = expression.args
    params = _parse_param_from_expr(expression, params, sketcher)
    
    # Get polyarc points from primitive, shape (k, 3)
    polyarc_points = PRIMITIVE_MAP[type(expression)](*params).to(device=sketcher.device)
    
    # Transform coordinates, shape (k, 2) -> (k, 2)
    polyarc_xy = polyarc_points[:, :2]
    transformed_xy = sketcher.get_coords(transform.inverse(), polyarc_xy)
    
    # Combine transformed xy with original bulge values
    transformed_points = th.cat([transformed_xy, polyarc_points[:, 2:]], dim=1)
    transformed_points = transformed_points.cpu().numpy().tolist()
    
    tuple_form: Tuple[Tuple[float, float, float], ...] = tuple(
        tuple(x) for x in transformed_points
    )
    new_expr = gls.PolyArc2D(tuple_form)
    return new_expr


def resolve_difference(expression: GLFunction) -> GLFunction:
    """
    Resolve Difference operations by converting them to Complement form.
    
    This function transforms the expression tree so that Difference operations
    are represented using Union, Intersection, and Complement. This simplifies
    subsequent DNF/CNF transformations.
    
    Args:
        expression: The input CSG expression.
    
    Returns:
        Equivalent expression with Difference resolved to Complement form.
    
    Raises:
        ValueError: If an unknown expression type is encountered.
    """
    # Inversion tracking
    inversion_mode: bool = False
    inversion_stack: List[bool] = [inversion_mode]

    execution_stack: List[GLFunction] = []
    execution_pointer_index: List[int] = []
    operator_stack: List[Type[GLFunction]] = []
    operator_nargs_stack: List[int] = []
    operator_params_stack: List[List[Any]] = []

    parser_list: List[GLFunction] = [expression]
    
    while parser_list:
        cur_expr = parser_list.pop()
        inversion_mode = inversion_stack.pop()
        
        if isinstance(cur_expr, COMBINATOR_TYPE):
            tree_branches: List[GLFunction] = []
            param_list: List[Any] = []
            
            for arg in cur_expr.args:
                if arg in cur_expr.lookup_table:
                    param_list.append(cur_expr.lookup_table[arg])
                else:
                    tree_branches.append(arg)
            n_args = len(tree_branches)

            # Handle inversion for different expression types
            if isinstance(cur_expr, gls.Difference):
                inversion_stack.append(not inversion_mode)
                inversion_stack.append(inversion_mode)
            elif isinstance(cur_expr, gls.Complement):
                inversion_stack.extend([not inversion_mode for _ in range(n_args)])
            else:
                inversion_stack.extend([inversion_mode for _ in range(n_args)])

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

        # Process completed operator calls
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
            
            if issubclass(operator, gls.Complement):
                new_canvas = args[0]
            else:
                new_canvas = operator(*args, *params)
            execution_stack = execution_stack[:-n_args] + [new_canvas]
    
    return execution_stack[0]
