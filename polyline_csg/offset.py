import geolipi.symbolic as gls
from .polyset import csg_to_polyset, is_valid_polyset, polyset_to_csg
import polyline_rs as prs
import numpy as np
from .polyset import upscale_polyexpr, downscale_polyset, UPSCALING_FACTOR

EPSILON = 1e-7
import torch
import math

def determine_polyline_orientation(poly, device="cuda"):
    """
    Determines whether a closed polyline is clockwise (CW) or counterclockwise (CCW).
    by computing a signed area. It uses the shoelace formula for chords + the 
    circular-segment area for arcs (determined by bulge).
    
    Args:
        polyline (tuple or list): Sequence of (x, y, bulge) points.
            - If the first and last points are not the same, it will be auto-closed.
    
    Returns:
        str: "CCW" if area > 0, "CW" if area < 0.
    """
    polyline_tensor = torch.tensor(poly.polyline, dtype=torch.float64, device=device)  # (N, 3)
    x, y, bulge = polyline_tensor[:, 0], polyline_tensor[:, 1], polyline_tensor[:, 2]

    # append first to end
    x = torch.cat([x, x[:1]])
    y = torch.cat([y, y[:1]])
    bulge = torch.cat([bulge, bulge[:1]])

    # Roll to get "next" point for each segment
    x_next = torch.roll(x, shifts=-1)
    y_next = torch.roll(y, shifts=-1)
    bulge_next = torch.roll(bulge, shifts=-1)

    # Shoelace area from the polygon formed by connecting chord endpoints
    shoelace_area = 0.5 * torch.sum(x * y_next - x_next * y)

    # Arc contributions
    # bulge.abs() -> which segments are arcs (vs. nearly straight lines)
    arc_mask = bulge.abs() > EPSILON
    arc_areas = torch.zeros_like(bulge)

    if arc_mask.any():
        # chord length
        dx = x_next - x
        dy = y_next - y
        chord_length = torch.sqrt(dx**2 + dy**2)

        # central angle of arc
        theta = 4.0 * torch.atan(bulge)

        # radius from chord + angle: R = c / (2 sin(theta/2))
        R = chord_length / (2.0 * torch.sin(theta / 2.0))

        # area of the circular segment: (1/2) * R^2 * (theta - sin(theta))
        raw_arc_areas = 0.5 * R**2 * (theta - torch.sin(theta))

        # sign depends on whether bulge is positive (CCW) or negative (CW)
        arc_areas = raw_arc_areas

        # For orientation, the sign is handled by the sign of `theta`, 
        # which is in turn set by bulge. If you'd prefer:
        # arc_areas *= torch.sign(bulge)   # But typically `theta` is negative 
        #                                  # if bulge is negative, so it's 
        #                                  # already accounted for.

        # Only include arc areas where arc_mask is True
        arc_areas[~arc_mask] = 0.0

    # Total area
    total_area = shoelace_area + torch.sum(arc_areas)

    return "CCW" if total_area > 0 else "CW"

def make_all_clockwise(polyset):
    
    for ind, poly in enumerate(polyset):
        poly_area = prs.compute_area(poly)
        if not np.sign(poly_area) == np.sign(poly.mode):
            reverse_sequence = get_reverse_sequence(poly.polyline)
            new_poly = prs.PolyStruct(polyline=reverse_sequence, is_closed=True, mode=poly.mode)
            polyset[ind] = new_poly
    return polyset

def get_reverse_sequence(polyline):
    xes = [x[0] for x in polyline]
    yes = [x[1] for x in polyline]
    bulges = [x[2] for x in polyline]
    reverse_x = xes[::-1]
    reverse_y = yes[::-1]
    reverse_bulges = bulges[::-1]
    # roll bulge
    reverse_bulges = reverse_bulges[1:] + reverse_bulges[:1]
    reverse_bulges = [-x for x in reverse_bulges]
    reverse_sequence = tuple(list(zip(reverse_x, reverse_y, reverse_bulges)))
    return reverse_sequence

def get_offset_expr(expression, offset=0.05, upscale=True, factor=UPSCALING_FACTOR):
    """
    Returns the offset expression of a given expression.
    
    Args:
        expr (GLFunction): The input expression.
        offset (float): The offset value.
    
    Returns:
        GLFunction: The offset expression.
    """
    expression = expression.sympy()
    if upscale:
        expression = upscale_polyexpr(expression, factor)
    polyset = csg_to_polyset(expression)
    # Make all clockwise is wrong. It should be alternative, based on sign of offset.
    polyset = make_all_clockwise(polyset)
    offset_polyset = []
    # NOTE CAN FAILE WHEN NOW OF THEM IS ZERO.
    for poly in polyset:
        try:
            cleaned_poly = prs.remove_redundant_vertices(poly, epsilon=1e-5)
            if cleaned_poly is None:     
                cleaned_poly = poly
            if upscale:
                # make sure its clean
                offset_curves = prs.offset_polyline(cleaned_poly, offset * factor, True)
            else:
                offset_curves = prs.offset_polyline(cleaned_poly, offset, True)
        except:
            print("FAILURE with offseting polyline")
            offset_curves = []
        
        for new_poly in offset_curves:
            offset_polyset.append(new_poly)
    if upscale:
        offset_polyset = downscale_polyset(offset_polyset, factor)
    offset_expr = polyset_to_csg(offset_polyset)
    return offset_expr