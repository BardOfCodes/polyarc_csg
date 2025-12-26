"""
PolyArc CSG - Exact signed distance field evaluation for 2D CSG expressions.

This package converts arbitrary 2D CSG expressions into intersection-free
boundary representations composed of line segments and circular arcs (PolyArcs).
"""

from .polyset import (
    UPSCALING_FACTOR,
    set_upscaling_factor,
    get_upscaling_factor,
    upscaling_context,
    construct_enclosure_tree,
    construct_enclosure_sequences,
    is_valid_polyset,
    polyset_to_csg,
    csg_to_polyset,
    clean_polyset,
    upscale_polyexpr,
    downscale_polyset,
    extract_primitives_with_signs,
    PolyArcCleaningError,
)

from .polyarc import (
    union,
    intersection,
    difference,
    xor,
    is_disjoint,
    is_overlapping,
    is_intersected,
    is_enclosed,
    is_1_inside_2,
    is_2_inside_1,
    union_multiple,
    intersection_multiple,
    difference_multiple,
    deduplicate_union,
    deduplicate_intersection,
)

from .valid_polyset import (
    expr_to_valid_polyset_expr,
    resolve_intersection,
    resolve_unions,
)

from .generic_to_polyset import (
    resolve_difference,
    expr_to_polyarc_expr,
)

from .offset import (
    get_offset_expr,
    determine_polyarc_orientation,
    make_all_clockwise,
    get_reverse_sequence,
)

__version__ = "0.1.0"

__all__ = [
    # Version
    "__version__",
    # Polyset operations
    "UPSCALING_FACTOR",
    "set_upscaling_factor",
    "get_upscaling_factor",
    "upscaling_context",
    "construct_enclosure_tree",
    "construct_enclosure_sequences",
    "is_valid_polyset",
    "polyset_to_csg",
    "csg_to_polyset",
    "clean_polyset",
    "upscale_polyexpr",
    "downscale_polyset",
    "extract_primitives_with_signs",
    "PolyArcCleaningError",
    # PolyArc boolean operations
    "union",
    "intersection",
    "difference",
    "xor",
    "is_disjoint",
    "is_overlapping",
    "is_intersected",
    "is_enclosed",
    "is_1_inside_2",
    "is_2_inside_1",
    "union_multiple",
    "intersection_multiple",
    "difference_multiple",
    "deduplicate_union",
    "deduplicate_intersection",
    # Valid polyset conversion
    "expr_to_valid_polyset_expr",
    "resolve_intersection",
    "resolve_unions",
    # Expression parsing
    "resolve_difference",
    "expr_to_polyarc_expr",
    # Offset operations
    "get_offset_expr",
    "determine_polyarc_orientation",
    "make_all_clockwise",
    "get_reverse_sequence",
]

