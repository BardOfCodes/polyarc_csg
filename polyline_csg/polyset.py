
import networkx as nx
import polyline_rs as prs
import geolipi.symbolic as gls
import woodie.symbolic as ws
from geolipi.symbolic.base_symbolic import GLFunction

UPSCALING_FACTOR = 1000  # Upscaling factor for higher resolution
def construct_enclosure_sequences(polyset):
    """
    Constructs enclosure sequences from a given PolySet by identifying immediate enclosure relationships.

    Args:
        polyset (list of prs.PolyStruct): The input list of PolyStructs.

    Returns:
        list of list of prs.PolyStruct: A list of enclosure sequences (from outermost to innermost).
    """

    G = nx.DiGraph()

    # Step 1: Add nodes (polylines)
    for poly in polyset:
        G.add_node(poly)

    # Step 2: Add directed edges based on enclosure
    for i, poly_a in enumerate(polyset):
        for j, poly_b in enumerate(polyset):
            if i != j and prs.is_1_inside_2(poly_b, poly_a):  # A encloses B
                G.add_edge(poly_a, poly_b)

    # Step 3: Apply transitive reduction to keep only immediate enclosures
    TR = nx.transitive_reduction(G)

    # Step 4: Identify root nodes (polylines with no incoming edges)
    roots = [n for n in TR.nodes if TR.in_degree(n) == 0]

    # Step 5: Extract enclosure sequences using DFS
    sequences = []

    def extract_sequences(node, path):
        """Recursively extract nested enclosure sequences."""
        path.append(node)
        if TR.out_degree(node) == 0:
            sequences.append(path.copy())  # Store sequence when reaching a leaf
        else:
            for child in TR.successors(node):
                extract_sequences(child, path)
        path.pop()

    for root in roots:
        extract_sequences(root, [])

    return sequences


def is_valid_polyset(polyset):
    """
    Checks if a given PolySet is valid.
    
    A valid PolySet must:
    1. Have non-intersecting polylines.
    2. Have sequences of alternating positive and negative mode values.
    
    Args:
        polyset (list of prs.PolyStruct): The input list of PolyStructs.
    
    Returns:
        bool: True if valid, False otherwise.
    """
    # Check for intersections
    for i, poly_a in enumerate(polyset):
        for j, poly_b in enumerate(polyset):
            if i != j:
                if prs.is_intersected(poly_a, poly_b):
                    return False  # Intersecting polylines are not allowed
                if prs.is_overlapping(poly_a, poly_b):
                    return False

    # Construct enclosure sequences
    sequences = construct_enclosure_sequences(polyset)

    # Ensure alternating signs in each sequence
    for seq in sequences:
        for k in range(len(seq) - 1):
            if seq[k].mode == seq[k + 1].mode:
                return False  # Consecutive elements should not have the same mode
    # Also No two should match
    
    return True  # Passed all checks


def polyset_to_csg(polyset):
    """
    Converts a PolySet (list of PolyStructs) into a CSG expression.
    
    Args:
        polyset (list of prs.PolyStruct]): The input PolySet.
    
    Returns:
        geolipi.symbolic.GLFunction: The corresponding CSG expression.
    """
    # Step 1: Construct the enclosure tree
    G = nx.DiGraph()

    # Step 1: Add nodes (polylines)
    for poly in polyset:
        G.add_node(poly)

    # Step 2: Add directed edges based on enclosure
    for i, poly_a in enumerate(polyset):
        for j, poly_b in enumerate(polyset):
            if i != j and prs.is_1_inside_2(poly_b, poly_a):  # A encloses B
                G.add_edge(poly_a, poly_b)

    # Step 3: Apply transitive reduction to keep only immediate enclosures
    TR = nx.transitive_reduction(G)


    # Identify root nodes (polylines with no enclosing parent)
    roots = [n for n in TR.nodes if TR.in_degree(n) == 0]

    # Step 2: Recursive function to convert tree to CSG expression
    def construct_csg(node):
        children = list(TR.successors(node))
        
        if not children:
            return gls.PolyLine2D(node.polyline)  # Leaf node is just a primitive

        # Recursively construct expressions for children
        child_exprs = [construct_csg(child) for child in children]

        # If the current node is positive: Difference with union of holes
        if len(child_exprs) == 1:
            child_expr = child_exprs[0]
        else:
            child_expr = gls.Union(*child_exprs)
        
        return gls.Difference(gls.PolyLine2D(node.polyline), child_expr)


    # Step 3: Handle multiple roots (disjoint top-level shapes)
    root_exprs = [construct_csg(root) for root in roots]
    root_exprs = [gls.Complement(expr) if roots[ind].mode == -1 else expr for ind, expr in enumerate(root_exprs)]

    # Combine multiple roots using Union
    # what to do in case of empty?
    if len(root_exprs) > 1:
        final_expr = gls.Union(*root_exprs)
    elif len(root_exprs) == 1:
        final_expr = root_exprs[0]
    else:
        final_expr = gls.NullExpression2D()
        
    return final_expr

def sp_tuple_to_tuple(polyline):
    """
    Reverts a PolyLine2D back to its argument form.
    
    Args:
        polyline (gls.PolyLine2D): The input PolyLine2D.
    
    Returns:
        tuple: The argument form of the PolyLine2D.
    """
    arg = tuple([(x[0], x[1], x[2]) for x in polyline])
    return arg

def extract_primitives_with_signs(expression, current_sign=1):
    """
    Recursively extracts primitives (gls.PolyLine2D) from the expression tree, tracking their signs.
    
    Args:
        expression (GLFunction): The CSG expression.
        current_sign (int): The current sign (+1 or -1) based on parent operations.
    
    Returns:
        List[Tuple[gls.PolyLine2D, int]]: A list of primitives with their corresponding signs.
    """
    primitives = []

    if isinstance(expression, gls.PolyLine2D):
        primitives.append((expression, current_sign))
    
    elif isinstance(expression, gls.Complement):
        primitives.extend(extract_primitives_with_signs(expression.args[0], -current_sign))

    elif isinstance(expression, gls.Difference):
        # Left child retains the sign, right child flips the sign
        primitives.extend(extract_primitives_with_signs(expression.args[0], current_sign))
        primitives.extend(extract_primitives_with_signs(expression.args[1], -current_sign))
    
    elif isinstance(expression, (gls.Union, gls.Intersection)):
        for child in expression.args:
            primitives.extend(extract_primitives_with_signs(child, current_sign))

    return primitives

def csg_to_polyset(expression: GLFunction):
    primitives = extract_primitives_with_signs(expression)
    polyset = []
    for prim in primitives:
        updated_points = sp_tuple_to_tuple(prim[0].args[0])
        polyset.append(prs.PolyStruct(polyline=updated_points, is_closed=True, mode=prim[1]))
    return polyset

def clean_polyset(polyset):
    """
    Cleans a PolySet by removing redundant primitives and empty
    """
    cleaned_polyset = []
    for poly in polyset:
        # Remove redundant vertices
        try:
            cleaned_poly = prs.remove_redundant_vertices(poly, epsilon=1e-5)
            if cleaned_poly is None:        
                failed_cleaning = True
                cleaned_poly = poly
            else:
                failed_cleaning = False
            bbox = prs.compute_extents(cleaned_poly)
            
            # Compute path length
            length = prs.compute_path_length(cleaned_poly)
            # Compute area
            area = prs.compute_area(cleaned_poly)

            condition_met = (abs(area) > 0.0) and (abs(length) > 0.0)
            if condition_met:
                cleaned_polyset.append(cleaned_poly)
        except:
            print("FAILURE with cleaning polyline")
            failed_cleaning = True
    return cleaned_polyset
        

def upscale_polyexpr(polyexpr, factor=UPSCALING_FACTOR):
    if isinstance(polyexpr, gls.PolyLine2D):
        points = polyexpr.args[0]
        upscaled_points = tuple([(x[0] * factor, x[1] * factor, x[2]) for x in points])
        new_polyline = gls.PolyLine2D(upscaled_points)
        return new_polyline
    else:    
        new_args = []
        for arg in polyexpr.args:
            out = upscale_polyexpr(arg, factor)
            new_args.append(out)
        return type(polyexpr)(*new_args)

def downscale_polyset(polyset, factor=UPSCALING_FACTOR):
    """
    Downscale a PolySet by a given factor.
    
    Args:
        polyset (list of prs.PolyStruct): The input PolySet.
        factor (int): The downscaling factor.
    
    Returns:
        list of prs.PolyStruct: The downscaled PolySet.
    """
    def downscale_polyline(polyline, factor):
        return tuple([(x[0] / factor, x[1] / factor, x[2]) for x in polyline])

    return [prs.PolyStruct(downscale_polyline(poly.polyline, factor), poly.is_closed, poly.mode) for poly in polyset]