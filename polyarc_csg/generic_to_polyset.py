"""
Expression parser for transforming CSG expressions into PolyArc2D primitives.

Handles transforms, parametric expressions, and primitive conversion.
"""
import torch as th
import geolipi.symbolic as gls
import sympy as sp
from typing import Dict, List, Any
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


def expr_to_polyarc_expr(expression: GLFunction, sketcher:Sketcher, transform: th.Tensor = None) -> GLFunction:
    if transform is None:
        transform = sketcher.get_affine_identity()
    resolved_expr = resolve_expression(expression, sketcher, transform)
    
    return resolved_expr

@singledispatch
def resolve_expression(expression: GLFunction, sketcher:Sketcher, transform: th.Tensor = None) -> GLFunction:
    raise NotImplementedError(f"Expression type {type(expression)} not implemented")

@resolve_expression.register
def resolve_expression_combinator(expression: COMBINATOR_TYPE, sketcher:Sketcher, transform: th.Tensor = None) -> GLFunction:

    tree_branches, param_list = [], []
    for arg in expression.args:
        if arg in expression.lookup_table:
            assert isinstance(arg, sp.Symbol), "Argument must be a symbol"
            param_list.append(expression.lookup_table[arg])
        else:
            tree_branches.append(arg)
    resolved_expr_list = []
    for child in tree_branches:
        cur_resolved_expr = resolve_expression(child, sketcher, transform.clone())
        resolved_expr_list.append(cur_resolved_expr)
    new_expr = type(expression)(*resolved_expr_list, *param_list)
    return new_expr

@resolve_expression.register
def resolve_expression_mod(expression: MOD_TYPE, sketcher:Sketcher, transform: th.Tensor = None) -> GLFunction:
    sub_expr = expression.args[0]
    params = expression.args[1:]
    params = _parse_param_from_expr(expression, params, sketcher)
    new_transform = MODIFIER_MAP[type(expression)](transform, *params)
    new_resolved_expr = resolve_expression(sub_expr, sketcher, new_transform)
    return new_resolved_expr

@resolve_expression.register
def resolve_expression_prim(expression: PRIM_TYPE, sketcher:Sketcher, transform: th.Tensor = None) -> GLFunction:
    
    params = expression.args
    params = _parse_param_from_expr(expression, params, sketcher)
    # shape k, 3
    polyarc_points = PRIMITIVE_MAP[type(expression)](*params).to(device=sketcher.device)
    # shape 3 x 3
    polyarc_xy = polyarc_points[:, :2]
    transformed_xy = sketcher.get_coords(transform.inverse(), polyarc_xy)
    transformed_points = th.cat([transformed_xy, polyarc_points[:, 2:]], dim=1)
    transformed_points = transformed_points.cpu().numpy().tolist()
    tuple_form  = tuple([tuple(x) for x in transformed_points])
    new_expr = gls.PolyArc2D(tuple_form)
    return new_expr


def resolve_difference(expression: GLFunction, ):
    """
    This can be improved.
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

            
            if isinstance(cur_expr, (gls.Difference)):
                inversion_stack.append(not inversion_mode)
                inversion_stack.append(inversion_mode)
            elif isinstance(cur_expr, (gls.Complement)):
                inversion_stack.extend([not inversion_mode for x in range(n_args)])
            else:
                inversion_stack.extend([inversion_mode for x in range(n_args)])

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
            if issubclass(operator, gls.Complement):
                new_canvas = args[0]
            else:
                new_canvas = operator(*args, *params)
            execution_stack = execution_stack[:-n_args] + [new_canvas]
    expression = execution_stack[0]
    return expression