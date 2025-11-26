
import torch as th
import geolipi.symbolic as gls
import woodie.symbolic as ws
import sympy as sp
from typing import Dict
from geolipi.symbolic.base import GLFunction
from geolipi.symbolic.symbol_types import PRIM_TYPE, COMBINATOR_TYPE, MOD_TYPE
from geolipi.torch_compute.sketcher import Sketcher
from geolipi.torch_compute.maps import INVERTED_MAP, NORMAL_MAP, MODIFIER_MAP

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


def resolve_to_transform_free_polyline_expr(expression, sketcher:Sketcher, uniforms: Dict[str, th.Tensor],):

    transforms_stack = [sketcher.get_affine_identity()]
    execution_stack = []
    operator_stack = []
    operator_nargs_stack = []
    operator_params_stack = []
    execution_pointer_index = []
    parser_list = [expression]
    device = sketcher.device


    while parser_list:
        cur_expr = parser_list.pop()
        if isinstance(cur_expr, COMBINATOR_TYPE):
            operator_stack.append(type(cur_expr))
            # what about parameterized combinators?
            tree_branches, cur_params = [], []
            for arg in cur_expr.args:
                if arg in cur_expr.lookup_table:
                    cur_params.append(cur_expr.lookup_table[arg])
                else:
                    tree_branches.append(arg)
            n_args = len(tree_branches)
            operator_nargs_stack.append(n_args)
            operator_params_stack.append(cur_params)
            transform = transforms_stack.pop()
            transform_chain = [transform.clone() for x in range(n_args)]
            transforms_stack.extend(transform_chain)
            next_to_parse = tree_branches[::-1]
            parser_list.extend(next_to_parse)
            execution_pointer_index.append(len(execution_stack))
            ### What if this also has a param???

        elif isinstance(cur_expr, MOD_TYPE):
            params = cur_expr.args[1:]
            params = recursive_parse_param(cur_expr, params, uniforms)
            # This is a hack unclear how to deal with other types)
            # if isinstance(cur_expr, (gls.EulerRotate2D)):
            #     params = [-params[0]]
            # elif isinstance(cur_expr, gls.Scale2D):
            #     params = [1.0 / params[0]]
            transform = transforms_stack.pop()
            identity_mat = sketcher.get_affine_identity()
            new_transform = MODIFIER_MAP[type(cur_expr)](identity_mat, *params)
            transform = th.matmul(new_transform, transform)
            transforms_stack.append(transform)
            next_to_parse = cur_expr.args[0]
            parser_list.append(next_to_parse)
        elif isinstance(cur_expr, PRIM_TYPE):
            # Here ideally we will also map the parameters to the correct values

            params = cur_expr.args
            params = recursive_parse_param(cur_expr, params, uniforms)
            # shape k, 3
            polyline_points = PRIMITIVE_MAP[type(cur_expr)](*params).to(device=device)
            # shape 3 x 3
            transform = transforms_stack.pop()
            polyline_xy = polyline_points[:, :2]
            transformed_xy = sketcher.get_coords(transform.inverse(), polyline_xy)
            transformed_points = th.cat([transformed_xy, polyline_points[:, 2:]], dim=1)
            transformed_points = transformed_points.cpu().numpy().tolist()
            tuple_form  = tuple([tuple(x) for x in transformed_points])
            new_expr = gls.PolyLine2D(tuple_form)
            # Get points using the primitive mapper. 
            # Apply transform to the points. 
            # Add the new primitive to the execution stack
            execution_stack.append(new_expr)
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

    assert len(execution_stack) == 1
    sdf = execution_stack[0]
    return sdf


def polyline_prim(params):
    params = th.stack(params, dim=0)
    return params

def rect_prim(params):
    size = params / 2.0
    points = th.tensor(
        [   [-size[0], -size[1], 0],
            [size[0], -size[1], 0],
            [size[0], size[1], 0],
            [-size[0], size[1], 0],
        ],
    )
    return points

def trapezoid_prim(r1, r2, he):
    # Construct the trapezoid corners in CCW order
    # bottom-left -> bottom-right -> top-right -> top-left
    return th.tensor([
        [-r2, -he, 0.0],
        [ r2, -he, 0.0],
        [ r1,  he, 0.0],
        [-r1,  he, 0.0],
    ], dtype=th.float32)
    
    return points

def circle_prim(params):
    r_small = th.sqrt(params/2.0)
    points = th.tensor([[-r_small[0], r_small[0], 1],
                         [r_small[0], -r_small[0], 1]],)
    return points

def no_param_rect_prim():
    points = th.tensor(
        [   [-0.5, -0.5, 0],
            [0.5, -0.5, 0],
            [0.5, 0.5, 0],
            [-0.5, 0.5, 0],
        ]
    )
    return points

PRIMITIVE_MAP = {
    gls.PolyLine2D: polyline_prim,
    gls.Rectangle2D: rect_prim,
    gls.Circle2D: circle_prim,
    gls.NoParamRectangle2D: no_param_rect_prim,
    gls.Trapezoid2D: trapezoid_prim,
}

uniform_type_map = {
    ws.UniformFloat: "float",
    ws.UniformVec2: "vec2",
    ws.UniformVec3: "vec3",
}

# A mapper from operation type to shader line template
# Now the same functions but for tuples: 
map_op_map_vec = {
    "ADD": lambda x, y: x + y,
    "SUB": lambda x, y: x - y,
    "MUL": lambda x, y: x * y,
    "DIV": lambda x, y: x / y,
    "POW": lambda x, y: x ** y,
    "NEG": lambda x: -x,
    "SIN": lambda x: th.sin(x),
    "COS": lambda x: th.cos(x),
    "TAN": lambda x: th.tan(x),
    "ASIN": lambda x: th.asin(x),
    "ACOS": lambda x: th.acos(x),
    "ATAN": lambda x: th.atan(x),
    "ATAN2": lambda x, y: th.atan2(x, y),
    "LOG": lambda x: th.log(x),
    "EXP": lambda x: th.exp(x),
    "SQRT": lambda x: th.sqrt(x),
    "ABS": lambda x: th.abs(x),
    "MIN": lambda x, y: th.min(x, y),
    "MAX": lambda x, y: th.max(x, y),
    "FLOOR": lambda x: th.floor(x),
    "CEIL": lambda x: th.ceil(x),
    "ROUND": lambda x: th.round(x),
    "FRAC": lambda x: th.frac(x),
    "SIGN": lambda x: th.sign(x),
    # "STEP": lambda x, y: th.step(x, y),
    "MOD": lambda x, y: th.mod(x, y),
    "NORMALIZE": lambda x: x / (th.norm(x) + 1e-9),
    "NORM": lambda x: th.norm(x),
}

def param_primitive_process(param, dtype=th.float32, device=th.device("cpu")):

    if isinstance(param, str):
        pass
    elif isinstance(param, (int, float, sp.Integer, sp.Float)):
        param = th.tensor(float(param), dtype=dtype, device=device)
    else:
        if len(param) == 0:
            if isinstance(param, sp.Integer):
                param = th.tensor(float(param), dtype=dtype, device=device)
        else:
            if isinstance(param[0], (list, tuple, sp.Tuple)):
                param = [th.tensor([float(y) for y in x], dtype=dtype, device=device) for x in param]
            else:
                param = th.tensor([float(x) for x in param], dtype=dtype, device=device)
    return param


def recursive_parse_param(expression , params, uniforms, dtype=th.float32, device=th.device("cpu")):

    shader_params = []
    for ind, param in enumerate(params):
        if isinstance(param, (tuple, sp.Tuple)):
            processed_param = param_primitive_process(param, dtype, device)
            shader_params.append(processed_param)
        elif isinstance(param, str):
            # Are we sure?
            shader_params.append(param)
        elif isinstance(param, sp.Symbol):
            if param in expression.lookup_table:
                cur_param = expression.lookup_table[param]
                processed_param = param_primitive_process(cur_param, dtype, device)
                shader_params.append(processed_param)
            else:
                shader_params.append(param.name)
        elif isinstance(param, (sp.Integer, sp.Float)):
            if isinstance(param, sp.Integer):
                param = th.tensor(float(param), dtype=dtype, device=device)
            shader_params.append(param)
        elif isinstance(param, (ws.VecList)):
            vector_list, n_vecs = param.args
            under_params = vector_list
            under_expression = param
            cur_params = recursive_parse_param(under_expression, under_params, uniforms)
            shader_params.append(cur_params)
        elif isinstance(param, (ws.UniformVec2, ws.UniformVec3)):
            min_val, default_val, max_val, uniform_name = param.args
            uniform_name = uniform_name.name
            if uniform_name in uniforms:
                param = uniforms[uniform_name]
            else:
                param = default_val
            param = param_primitive_process(param, dtype, device)
            shader_params.append(param)
        elif isinstance(param, (ws.UniformFloat)):
            min_val, default_val, max_val, uniform_name = param.args
            uniform_name = uniform_name.name
            if uniform_name in uniforms:
                param = uniforms[uniform_name]
            else:
                param = default_val
            param = param_primitive_process(param, dtype, device)
            param = param.unsqueeze(0)
            shader_params.append(param)
        elif isinstance(param, (ws.Vec2, ws.Vec3, ws.Vec4)):
            # Now its input can be a math node, or a variable. 
            under_expression = param
            under_params = param.args
            cur_params = recursive_parse_param(under_expression, under_params, uniforms)
            cur_params = th.cat(cur_params, dim=0)
            shader_params.append(cur_params)
        elif isinstance(param, (ws.Float,)):
            # Now its input can be a math node, or a variable. 
            under_expression = param
            under_params = param.args
            cur_params = recursive_parse_param(under_expression, under_params, uniforms)
            shader_params.append(cur_params[0])
        elif isinstance(param, ws.VarSplitter):
            under_expression = param
            under_params = param.args
            cur_params = recursive_parse_param(under_expression, under_params, uniforms)
            selected_ind = int(float(cur_params[1]))
            new_param = cur_params[0]
            shader_params.append(new_param[selected_ind])
            
        elif isinstance(param, (ws.UnaryOperator, ws.VectorOperator)):
            under_expression = param
            under_params = param.args
            cur_params = recursive_parse_param(under_expression, under_params, uniforms)
            new_param = cur_params[0]
            op = cur_params[1]
            op_func = map_op_map_vec[op]
            param = op_func(new_param)
            shader_params.append(param)
        elif isinstance(param, ws.BinaryOperator):
            under_expression = param
            under_params = param.args
            cur_params = recursive_parse_param(under_expression, under_params, uniforms)
            new_param1 = cur_params[0]
            new_param2 = cur_params[1]
            op = cur_params[2]
            op_func = map_op_map_vec[op]
            param = op_func(new_param1, new_param2)
            shader_params.append(param)
        else:
            raise NotImplementedError
            
    return shader_params

