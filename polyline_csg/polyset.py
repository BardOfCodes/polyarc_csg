
import networkx as nx
import polyline_rs as prs
import geolipi.symbolic as gls
import woodie.symbolic as ws

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
            if i != j and prs.is_intersected(poly_a, poly_b):
                return False  # Intersecting polylines are not allowed

    # Construct enclosure sequences
    sequences = construct_enclosure_sequences(polyset)

    # Ensure alternating signs in each sequence
    for seq in sequences:
        for k in range(len(seq) - 1):
            if seq[k].mode == seq[k + 1].mode:
                return False  # Consecutive elements should not have the same mode

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
            return ws.PolyLine2D(tuple(node.polyline))  # Leaf node is just a primitive

        # Recursively construct expressions for children
        child_exprs = [construct_csg(child) for child in children]

        # If the current node is positive: Difference with union of holes
        if len(child_exprs) == 1:
            child_expr = child_exprs[0]
        else:
            child_expr = gls.Union(*child_exprs)
        
        return gls.Difference(ws.PolyLine2D(tuple(node.polyline)), child_expr)


    # Step 3: Handle multiple roots (disjoint top-level shapes)
    root_exprs = [construct_csg(root) for root in roots]
    root_exprs = [gls.Complement(expr) if roots[ind].mode == -1 else expr for ind, expr in enumerate(root_exprs)]

    # Combine multiple roots using Union
    final_expr = gls.Union(*root_exprs) if len(root_exprs) > 1 else root_exprs[0]

    return final_expr