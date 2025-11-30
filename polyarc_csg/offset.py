"""Offset operations for PolySets."""
import logging
from typing import List

import numpy as np
import torch
import polyarc_rs as prs
from geolipi.symbolic.base import GLFunction

from .polyset import (
    csg_to_polyset, 
    polyset_to_csg, 
    upscale_polyexpr, 
    downscale_polyset, 
    UPSCALING_FACTOR
)

__all__ = [
    "get_offset_expr",
    "determine_polyarc_orientation",
    "make_all_clockwise",
    "get_reverse_sequence",
]

logger = logging.getLogger(__name__)

EPSILON = 1e-7


def determine_polyarc_orientation(poly: prs.PolyArc, device: str = "cuda") -> str:
    """
    Determines whether a closed polyarc is clockwise (CW) or counterclockwise (CCW).
    
    Uses the shoelace formula for chords plus the circular-segment area for arcs.
    
    Args:
        poly: The polyarc to check.
        device: Torch device to use for computation.
    
    Returns:
        "CCW" if area > 0, "CW" if area < 0.
    """
    polyarc_tensor = torch.tensor(poly.polyarc, dtype=torch.float64, device=device)
    x, y, bulge = polyarc_tensor[:, 0], polyarc_tensor[:, 1], polyarc_tensor[:, 2]

    # Append first to end for closure
    x = torch.cat([x, x[:1]])
    y = torch.cat([y, y[:1]])
    bulge = torch.cat([bulge, bulge[:1]])

    # Roll to get "next" point for each segment
    x_next = torch.roll(x, shifts=-1)
    y_next = torch.roll(y, shifts=-1)

    # Shoelace area
    shoelace_area = 0.5 * torch.sum(x * y_next - x_next * y)

    # Arc contributions
    arc_mask = bulge.abs() > EPSILON
    arc_areas = torch.zeros_like(bulge)

    if arc_mask.any():
        dx = x_next - x
        dy = y_next - y
        chord_length = torch.sqrt(dx**2 + dy**2)
        theta = 4.0 * torch.atan(bulge)
        R = chord_length / (2.0 * torch.sin(theta / 2.0))
        raw_arc_areas = 0.5 * R**2 * (theta - torch.sin(theta))
        arc_areas = raw_arc_areas
        arc_areas[~arc_mask] = 0.0

    total_area = shoelace_area + torch.sum(arc_areas)
    return "CCW" if total_area > 0 else "CW"


def get_reverse_sequence(polyarc: tuple) -> tuple:
    """
    Reverses the vertex sequence of a polyarc, adjusting bulges accordingly.
    
    Args:
        polyarc: The polyarc vertices as a tuple of (x, y, bulge).
        
    Returns:
        The reversed polyarc.
    """
    xes = [x[0] for x in polyarc]
    yes = [x[1] for x in polyarc]
    bulges = [x[2] for x in polyarc]
    
    reverse_x = xes[::-1]
    reverse_y = yes[::-1]
    reverse_bulges = bulges[::-1]
    
    # Roll and negate bulge
    reverse_bulges = reverse_bulges[1:] + reverse_bulges[:1]
    reverse_bulges = [-x for x in reverse_bulges]
    
    return tuple(zip(reverse_x, reverse_y, reverse_bulges))


def make_all_clockwise(polyset: List[prs.PolyArc]) -> List[prs.PolyArc]:
    """
    Ensures all polyarcs have consistent orientation matching their mode.
    
    Positive mode polyarcs should have positive area (CCW),
    negative mode polyarcs should have negative area (CW).
    
    Args:
        polyset: The input PolySet.
        
    Returns:
        The PolySet with corrected orientations.
    """
    result = []
    for poly in polyset:
        poly_area = prs.compute_area(poly)
        if not np.sign(poly_area) == np.sign(poly.mode):
            reverse_sequence = get_reverse_sequence(poly.polyarc)
            poly = prs.PolyArc(polyarc=reverse_sequence, is_closed=True, mode=poly.mode)
        result.append(poly)
    return result


def get_offset_expr(
    expression: GLFunction, 
    offset: float = 0.05, 
    upscale: bool = True, 
    factor: int = UPSCALING_FACTOR
) -> GLFunction:
    """
    Returns the offset expression of a given expression.
    
    Args:
        expression: The input CSG expression.
        offset: The offset distance (positive = outward, negative = inward).
        upscale: Whether to upscale for numerical precision.
        factor: The upscaling factor.
    
    Returns:
        The offset CSG expression.
    """
    expression = expression.sympy()
    
    if upscale:
        expression = upscale_polyexpr(expression, factor)
    
    polyset = csg_to_polyset(expression)
    polyset = make_all_clockwise(polyset)
    
    offset_polyset = []
    
    for poly in polyset:
        try:
            cleaned_poly = prs.remove_redundant_vertices(poly, epsilon=1e-5)
            if cleaned_poly is None:
                cleaned_poly = poly
            
            actual_offset = offset * factor if upscale else offset
            offset_curves = prs.offset_polyarc(cleaned_poly, actual_offset, True)
            
            offset_polyset.extend(offset_curves)
        except Exception as e:
            logger.warning(f"Failed to offset polyarc: {e}")
    
    if upscale:
        offset_polyset = downscale_polyset(offset_polyset, factor)
    
    return polyset_to_csg(offset_polyset)
