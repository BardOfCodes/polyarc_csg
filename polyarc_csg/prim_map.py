"""
Primitive mapping functions for converting CSG primitives to PolyArc vertices.

Each function takes primitive parameters and returns a tensor of vertices
in the format (x, y, bulge) where bulge=0 for line segments and bulge!=0 for arcs.
"""
from typing import Dict, List, Tuple, Type, Callable, Union

import torch as th
import geolipi.symbolic as gls
from geolipi.symbolic.base import GLFunction

__all__ = [
    "PRIMITIVE_MAP",
]

# Type alias for primitive functions
PrimitiveFunction = Callable[..., th.Tensor]


def polyarc_prim(params: List[th.Tensor]) -> th.Tensor:
    """
    Convert PolyArc2D parameters directly to a tensor.
    
    Args:
        params: List of tensors representing polyarc vertices.
    
    Returns:
        Stacked tensor of shape (N, 3) with (x, y, bulge) per vertex.
    """
    params = th.stack(params, dim=0)
    return params


def rect_prim(params: th.Tensor) -> th.Tensor:
    """
    Create rectangle vertices from size parameters.
    
    Args:
        params: Tensor of shape (2,) with (width, height).
    
    Returns:
        Tensor of shape (4, 3) with rectangle vertices in CCW order.
    """
    size = params / 2.0
    points = th.tensor(
        [   [-size[0], -size[1], 0],
            [size[0], -size[1], 0],
            [size[0], size[1], 0],
            [-size[0], size[1], 0],
        ],
    )
    return points


def trapezoid_prim(r1: th.Tensor, r2: th.Tensor, he: th.Tensor) -> th.Tensor:
    """
    Construct a trapezoid with given parameters.
    
    Args:
        r1: Half-width at the top edge.
        r2: Half-width at the bottom edge.
        he: Half-height of the trapezoid.
    
    Returns:
        Tensor of shape (4, 3) with vertices in CCW order.
    """
    return th.tensor([
        [-r2, -he, 0.0],
        [ r2, -he, 0.0],
        [ r1,  he, 0.0],
        [-r1,  he, 0.0],
    ], dtype=th.float32)


def circle_prim(params: th.Tensor) -> th.Tensor:
    """
    Create circle vertices (two semicircles using bulge=1).
    
    Args:
        params: Tensor of shape (1,) with radius.
    
    Returns:
        Tensor of shape (2, 3) representing a circle via two semicircular arcs.
    """
    r_small = params
    points = th.tensor([[-r_small[0], 0, 1],
                         [r_small[0], 0, 1]],)
    return points


def no_param_rect_prim() -> th.Tensor:
    """
    Create a unit square centered at origin.
    
    Returns:
        Tensor of shape (4, 3) with unit square vertices.
    """
    points = th.tensor(
        [   [-0.5, -0.5, 0],
            [0.5, -0.5, 0],
            [0.5, 0.5, 0],
            [-0.5, 0.5, 0],
        ]
    )
    return points


def triangle_prim_params(params: Tuple[th.Tensor, th.Tensor, th.Tensor]) -> th.Tensor:
    """
    Create triangle vertices from three corner points.
    
    Args:
        params: Tuple of three tensors, each of shape (2,) with (x, y) coordinates.
    
    Returns:
        Tensor of shape (3, 3) with triangle vertices.
    """
    p0, p1, p2 = params
    points = th.tensor(
        [[p0[0], p0[1], 0],
         [p1[0], p1[1], 0],
         [p2[0], p2[1], 0],
        ]
    )
    return points


# Mapping from primitive types to their vertex generation functions
PRIMITIVE_MAP: Dict[Type[GLFunction], PrimitiveFunction] = {
    gls.PolyArc2D: polyarc_prim,
    gls.Rectangle2D: rect_prim,
    gls.Circle2D: circle_prim,
    gls.NoParamRectangle2D: no_param_rect_prim,
    gls.Trapezoid2D: trapezoid_prim,
    gls.Triangle2D: triangle_prim_params,
}

# Potential Primitives: 
# Circle
# Rounded Box
# Chamfer Box
# Box
# Oriented Box
# Segment
# Rhombus
# Isosceles Trapezoid
# Parallelogram
# Equilateral Triangle
# Isosceles Triangle
# Triangle
# Uneven Capsule
# Regular Pentagon
# Regular Hexagon
# Regular Octagon
# Hexagram
# Pentagram
# Regular Star
# Pie
# Cut Disk
# Arc
# Ring
# Vecisa
# Moon
# Heart
# Cross
