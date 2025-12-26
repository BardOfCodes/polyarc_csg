"""Boolean operations on lists of PolyArcs."""
from typing import List, Literal
import polyarc_rs as prs

# Type alias for boolean operation modes
BooleanMode = Literal["union", "intersection", "difference", "xor"]

__all__ = [
    "union", "intersection", "difference", "xor",
    "is_disjoint", "is_overlapping", "is_intersected", "is_enclosed",
    "is_1_inside_2", "is_2_inside_1",
    "union_multiple", "intersection_multiple", "difference_multiple",
    "deduplicate_union", "deduplicate_intersection",
]


def union(a: prs.PolyArc, b: prs.PolyArc) -> List[prs.PolyArc]:
    """Returns the union of two polyarcs."""
    return prs.boolean_operation(a, b, "union")


def intersection(a: prs.PolyArc, b: prs.PolyArc) -> List[prs.PolyArc]:
    """Returns the intersection of two polyarcs."""
    return prs.boolean_operation(a, b, "intersection")


def difference(a: prs.PolyArc, b: prs.PolyArc) -> List[prs.PolyArc]:
    """Returns the difference (A - B) of two polyarcs."""
    return prs.boolean_operation(a, b, "difference")


def xor(a: prs.PolyArc, b: prs.PolyArc) -> List[prs.PolyArc]:
    """Returns the exclusive OR (symmetric difference) of two polyarcs."""
    return prs.boolean_operation(a, b, "xor")


def is_disjoint(a: prs.PolyArc, b: prs.PolyArc) -> bool:
    """Returns True if the two polyarcs are disjoint (no overlap)."""
    return prs.is_disjoint(a, b)


def is_overlapping(a: prs.PolyArc, b: prs.PolyArc) -> bool:
    """Returns True if polyarc A overlaps exactly with polyarc B."""
    return prs.is_overlapping(a, b)


def is_intersected(a: prs.PolyArc, b: prs.PolyArc) -> bool:
    """Returns True if polyarc A intersects with polyarc B."""
    return prs.is_intersected(a, b)


def is_enclosed(a: prs.PolyArc, b: prs.PolyArc) -> bool:
    """Returns True if one polyarc is completely inside the other."""
    return prs.is_enclosed(a, b)


def is_1_inside_2(a: prs.PolyArc, b: prs.PolyArc) -> bool:
    """Returns True if polyarc A is completely inside polyarc B."""
    return prs.is_1_inside_2(a, b)


def is_2_inside_1(a: prs.PolyArc, b: prs.PolyArc) -> bool:
    """Returns True if polyarc B is completely inside polyarc A."""
    return prs.is_2_inside_1(a, b)


def deduplicate_union(polyarcs: List[prs.PolyArc], bool_mode: BooleanMode = "union") -> List[prs.PolyArc]:
    """
    Removes redundant polyarcs from the list by checking for:
    - Exact overlaps
    - Enclosed polyarcs (fully inside another)
    - Intersecting polyarcs (merged together)
    
    Args:
        polyarcs: List of polyarcs to deduplicate.
        bool_mode: Boolean operation mode for merging.

    Returns:
        Deduplicated list of polyarcs.
    """
    unique_polyarcs = []
    
    for poly in polyarcs:
        is_redundant = False

        for unique_poly in unique_polyarcs:
            if prs.is_overlapping(poly, unique_poly):
                is_redundant = True
                break
            if prs.is_1_inside_2(poly, unique_poly):
                is_redundant = True
                break
            if prs.is_intersected(poly, unique_poly):
                merged = prs.boolean_operation(poly, unique_poly, bool_mode)
                unique_polyarcs.remove(unique_poly)
                unique_polyarcs.extend(merged)
                is_redundant = True
                break

        if not is_redundant:
            unique_polyarcs.append(poly)

    return unique_polyarcs


def deduplicate_intersection(polyarcs: List[prs.PolyArc]) -> List[prs.PolyArc]:
    """
    Removes duplicate (overlapping) polyarcs from the list.

    Args:
        polyarcs: List of polyarcs.

    Returns:
        List with duplicates removed.
    """
    unique_polyarcs = []

    for poly in polyarcs:
        is_duplicate = False
        for unique_poly in unique_polyarcs:
            if prs.is_overlapping(poly, unique_poly):
                is_duplicate = True
                break

        if not is_duplicate:
            unique_polyarcs.append(poly)

    return unique_polyarcs


def union_multiple(polyarcs: List[prs.PolyArc]) -> List[prs.PolyArc]:
    """
    Performs union operation on a list of polyarcs.

    Args:
        polyarcs: The polyarcs to union.

    Returns:
        Resulting polyarcs after the union.
    """
    if not polyarcs:
        return []

    result_polyarcs = [polyarcs[0]]

    for poly in polyarcs[1:]:
        new_results = []

        for res_poly in result_polyarcs:
            union_result = prs.boolean_operation(res_poly, poly, "union")
            new_results.extend(union_result)

        if not new_results:
            new_results.append(poly)

        result_polyarcs = deduplicate_union(new_results)

    return result_polyarcs


def intersection_multiple(polyarcs: List[prs.PolyArc]) -> List[prs.PolyArc]:
    """
    Perform intersection across multiple polyarcs.

    Args:
        polyarcs: List of polyarcs to intersect.

    Returns:
        Resulting polyarcs after intersection.
    """
    if not polyarcs:
        return []
    if len(polyarcs) == 1:
        return polyarcs

    result = [polyarcs[0]]

    for poly in polyarcs[1:]:
        new_result = []
        for res_poly in result:
            intersection_result = prs.boolean_operation(res_poly, poly, "intersection")
            new_result.extend(intersection_result)

        result = deduplicate_intersection(new_result)

        if not result:
            break

    return result


def difference_multiple(A: List[prs.PolyArc], B: List[prs.PolyArc]) -> List[prs.PolyArc]:
    """
    Perform the difference between a list of polyarcs (A) and another list (B).

    Args:
        A: Positive polyarcs.
        B: Polyarcs to subtract.

    Returns:
        Resulting polyarcs after difference (A - B).
    """
    result_polyarcs = []

    for a in A:
        working_polyarcs = [a]

        for b in B:
            new_working_polyarcs = []

            for poly in working_polyarcs:
                diff_result = difference(poly, b)
                new_working_polyarcs.extend(diff_result)

            working_polyarcs = new_working_polyarcs

        result_polyarcs.extend(working_polyarcs)

    return result_polyarcs

