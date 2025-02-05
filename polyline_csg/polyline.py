import polyline_rs as prs

def union(a: prs.PolyStruct, b: prs.PolyStruct):
    """Returns the union of two polylines."""
    return prs.boolean_operation(a, b, "union")

def intersection(a: prs.PolyStruct, b: prs.PolyStruct):
    """Returns the intersection of two polylines."""
    return prs.boolean_operation(a, b, "intersection")

def difference(a: prs.PolyStruct, b: prs.PolyStruct):
    """Returns the difference (A - B) of two polylines."""
    return prs.boolean_operation(a, b, "difference")

def xor(a: prs.PolyStruct, b: prs.PolyStruct):
    """Returns the exclusive OR (symmetric difference) of two polylines."""
    return prs.boolean_operation(a, b, "xor")

def is_disjoint(a: prs.PolyStruct, b: prs.PolyStruct) -> bool:
    """Returns True if the two polylines are disjoint (no overlap)."""
    return prs.is_disjoint(a, b)

def is_overlapping(a: prs.PolyStruct, b: prs.PolyStruct) -> bool:
    """Returns True if polyline A overlaps with polyline B."""
    return prs.is_overlap(a, b)

def is_intersected(a: prs.PolyStruct, b: prs.PolyStruct) -> bool:
    """Returns True if polyline A overlaps with polyline B."""
    return prs.is_intersected(a, b)

def is_enclosed(a: prs.PolyStruct, b: prs.PolyStruct) -> bool:
    """Returns True if polyline A is completely inside polyline B."""
    return prs.is_enclosed(a, b)

def is_1_inside_2(a: prs.PolyStruct, b: prs.PolyStruct) -> bool:
    """Returns True if polyline A is completely inside polyline B."""
    return prs.is_1_inside_2(a, b)

def is_2_inside_1(a: prs.PolyStruct, b: prs.PolyStruct) -> bool:
    """Returns True if polyline A is completely inside polyline B."""
    return prs.is_2_inside_1(a, b)

def deduplicate_union(polylines, bool_mode="union"):
    """
    Removes redundant polylines from the list by checking for:
    - Exact overlaps
    - Enclosed polylines (fully inside another)
    
    Args:
        polylines (list of prs.PolyStruct): List of polylines to deduplicate.

    Returns:
        list of prs.PolyStruct: Deduplicated list of polylines.
    """
    unique_polylines = []
    
    for poly in polylines:
        is_redundant = False

        for unique_poly in unique_polylines:
            if prs.is_overlapping(poly, unique_poly):
                is_redundant = True  # Exact overlap
                break
            if prs.is_1_inside_2(poly, unique_poly):
                is_redundant = True  # Fully enclosed inside
                break
            if prs.is_intersected(poly, unique_poly):
                # Merge intersecting polylines
                merged = prs.boolean_operation(poly, unique_poly, bool_mode)
                unique_polylines.remove(unique_poly)
                unique_polylines.extend(merged)
                is_redundant = True
                break

        if not is_redundant:
            unique_polylines.append(poly)

    return unique_polylines

def union_multiple(polylines):
    """
    Performs union operation on a list of polylines.

    Args:
        polylines (list of prs.PolyStruct): The polylines to union.

    Returns:
        list of prs.PolyStruct: Resulting polylines after the union.
    """
    if not polylines:
        return []

    result_polylines = [polylines[0]]

    for poly in polylines[1:]:
        new_results = []

        # Union with each existing polyline in the result
        for res_poly in result_polylines:
            union_result = prs.boolean_operation(res_poly, poly, "union")
            new_results.extend(union_result)

        # Add the current poly if it didn't merge with any existing ones
        if not new_results:
            new_results.append(poly)

        # Deduplicate after each union operation
        result_polylines = deduplicate_union(new_results)

    return result_polylines


def intersection_multiple(polylines: list) -> list:
    """
    Perform intersection across multiple polylines.

    Args:
        polylines (list of prs.PolyStruct): List of polylines to intersect.

    Returns:
        list of prs.PolyStruct: Resulting polylines after intersection.
    """
    if not polylines:
        return []
    if len(polylines) == 1:
        return polylines

    # Start with the first polyline
    result = [polylines[0]]

    for poly in polylines[1:]:
        new_result = []
        for res_poly in result:
            # Perform intersection between current result polyline and next polyline
            intersection_result = prs.boolean_operation(res_poly, poly, "intersection")
            new_result.extend(intersection_result)

        # Remove duplicate polylines
        result = deduplicate_intersection(new_result)

        # Early exit: If at any point the intersection is empty, the final result will be empty
        if not result:
            break

    return result

def deduplicate_intersection(polylines: list) -> list:
    """
    Removes duplicate (overlapping) polylines from the list.

    Args:
        polylines (list of prs.PolyStruct): List of polylines.

    Returns:
        list of prs.PolyStruct: List with duplicates removed.
    """
    unique_polylines = []

    for poly in polylines:
        is_duplicate = False
        for unique_poly in unique_polylines:
            if prs.is_overlapping(poly, unique_poly):
                is_duplicate = True
                break  # No need to check further

        if not is_duplicate:
            unique_polylines.append(poly)

    return unique_polylines


import polyline_rs as prs

def difference(a_list, b_list):
    """Performs difference between two polylines."""
    return prs.boolean_operation(a_list, b_list, "difference")

def difference_multiple(A, B):
    """
    Perform the difference between a list of polylines (A) and another list (B).

    Args:
        A (list of prs.PolyStruct): Positive polylines.
        B (list of prs.PolyStruct): Positive polylines to subtract.

    Returns:
        list of prs.PolyStruct: Resulting polylines after difference.
    """
    result_polylines = []

    for a in A:
        working_polylines = [a]

        for b in B:
            new_working_polylines = []

            for poly in working_polylines:
                # Subtract b from poly
                diff_result = difference(poly, b)
                new_working_polylines.extend(diff_result)

            working_polylines = new_working_polylines  # Update for next b

        # Add remaining parts to result
        result_polylines.extend(working_polylines)

    return result_polylines
