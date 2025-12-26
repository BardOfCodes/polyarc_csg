import torch as th
import geolipi.symbolic as gls

__all__ = [
    "PRIMITIVE_MAP",
]


def polyarc_prim(params):
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
    r_small = params # th.sqrt(params/2.0)
    points = th.tensor([[-r_small[0], 0, 1],
                         [r_small[0], 0, 1]],)
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

def triangle_prim_params(params):
    p0, p1, p2 = params
    points = th.tensor(
        [[p0[0], p0[1], 0],
         [p1[0], p1[1], 0],
         [p2[0], p2[1], 0],
        ]
    )
    return points



PRIMITIVE_MAP = {
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