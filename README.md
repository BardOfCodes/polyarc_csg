# polyarc_csg

Higher-level CSG (Constructive Solid Geometry) operations on **PolySets** - collections of closed PolyArcs representing complex 2D shapes with holes. Built on top of [polyarc_rs](../polyarc_rs) and [geolipi](https://github.com/bardofcodes/geolipi) for symbolic expression handling.

## What is a PolySet?

A **PolySet** is a list of `PolyArc` objects that together define a 2D region, potentially with:
- Multiple disjoint outer boundaries
- Nested holes (alternating positive/negative modes)
- Complex enclosure hierarchies

For example, a donut is a PolySet with two PolyArcs: an outer circle (mode=+1) and an inner hole (mode=-1).

## Features

- **PolySet validation**: Check if a PolySet is valid (no intersections, proper mode alternation)
- **CSG ↔ PolySet conversion**: Convert between geolipi symbolic expressions and PolySets
- **Boolean operations on PolySets**: Union, intersection, difference of multiple shapes
- **Expression parsing**: Convert arbitrary CSG expressions (with transforms) into valid PolySets
- **Offset operations**: Parallel offset of complex shapes

## Dependencies

- `polyarc_rs`: Low-level polyarc operations
- `geolipi`: Symbolic CSG expression framework
- `networkx`: Enclosure tree construction
- `numpy`, `torch`: Numerical operations

## Installation

```bash
# Install dependencies
pip install networkx numpy torch

# Install polyarc_rs (see ../polyarc_rs)
cd ../polyarc_rs && maturin develop --release

# Install geolipi (follow geolipi installation instructions)
```

## Usage

### PolySet Validation

```python
import polyarc_rs as prs
from polyarc_csg.polyset import is_valid_polyset, polyset_to_csg

# Create a donut shape
outer = prs.PolyArc(polyarc=((-10,-10,0),(10,-10,0),(10,10,0),(-10,10,0)), is_closed=True, mode=1)
hole = prs.PolyArc(polyarc=((-5,-5,0),(5,-5,0),(5,5,0),(-5,5,0)), is_closed=True, mode=-1)

polyset = [outer, hole]
print(is_valid_polyset(polyset))  # True
```

### PolySet to CSG Expression

```python
from polyarc_csg.polyset import polyset_to_csg

csg_expr = polyset_to_csg(polyset)
# Returns: gls.Difference(gls.PolyArc2D(...), gls.PolyArc2D(...))
```

### CSG Expression to PolySet

```python
import geolipi.symbolic as gls
from polyarc_csg.polyset import csg_to_polyset

expr = gls.Difference(
    gls.PolyArc2D(((-10,-10,0),(10,-10,0),(10,10,0),(-10,10,0))),
    gls.PolyArc2D(((-5,-5,0),(5,-5,0),(5,5,0),(-5,5,0)))
)
polyset = csg_to_polyset(expr)
```

### Complex CSG Parsing

```python
from polyarc_csg.polyset_parser import parse_csg_to_valid_polyset_csg

# Handles transforms, nested operations, and produces valid PolySets
expr = gls.Union(
    gls.Difference(outer_shape, hole1),
    gls.Difference(outer_shape2, hole2)
)
result_csg = parse_csg_to_valid_polyset_csg(expr, sketcher, uniforms)
```

### Boolean Operations on PolySets

```python
from polyarc_csg.polyarc import union_multiple, intersection_multiple, difference_multiple

# Union of multiple polyarcs
result = union_multiple([polyarc1, polyarc2, polyarc3])

# Difference: A - B (each is a list of polyarcs)
result = difference_multiple(positive_list, negative_list)
```

### Offset Operations

```python
from polyarc_csg.offset import get_offset_expr

offset_expr = get_offset_expr(csg_expression, offset=0.1)
```

## Module Structure

| Module | Purpose |
|--------|---------|
| `polyset.py` | PolySet validation, CSG↔PolySet conversion, enclosure tree |
| `polyarc.py` | Boolean operations on lists of PolyArcs |
| `polyset_parser.py` | Parse complex CSG to valid PolySets (handles DNF/CNF) |
| `expression_parser.py` | Transform primitives to PolyArc2D, handle parametric expressions |
| `offset.py` | Offset operations on polysets |

## Key Concepts

### Valid PolySet Rules
1. No two PolyArcs may intersect or overlap
2. Enclosure sequences must have alternating modes (+1, -1, +1, ...)
3. Outer boundaries have mode=+1, holes have mode=-1

### Upscaling
Operations use a `UPSCALING_FACTOR=1000` to improve numerical precision. Coordinates are scaled up before operations and scaled down after.

## Running Tests

```bash
# Install dev dependencies
pip install pytest

# Run tests
pytest tests/ -v
```

## License

MIT

