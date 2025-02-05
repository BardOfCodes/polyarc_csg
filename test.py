import polyline_rs
import polyline_rs

# Define positive and negative PolyStructs

# Define some sample PolyStructs
# poly_a = polyline_rs.PolyStruct(1, [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
# poly_b = polyline_rs.PolyStruct(1, [(0.5, 0.5, 0), (1.5, 0.5, 0), (1.5, 1.5, 0), (0.5, 1.5, 0)])
# bounding_box = polyline_rs.PolyStruct(1, [(-2, -2, 0), (2, -2, 0), (2, 2, 0), (-2, 2, 0)])

# # Test Union
# union_result = polyline_rs.union([poly_a, poly_b])
# print("Union Result:", [(p.poly_type, p.polyline) for p in union_result])

# result = union_result

# Test Intersection
# intersection_result = polyline_rs.intersection([poly_a, poly_b])
# print("Intersection Result:", [(p.poly_type, p.polyline) for p in intersection_result])

# # Test Difference
# difference_result = polyline_rs.difference(poly_a, poly_b)
# print("Difference Result:", [(p.poly_type, p.polyline) for p in difference_result])

# # Test Complement
# complement_result = polyline_rs.complement(poly_b, bounding_box)
# print("Complement Result:", [(p.poly_type, p.polyline) for p in complement_result])

# # # Define sample contours
# square = polyline_rs.PolyStruct(1, [(-0.5, -0.5, 0.0), (-0.5, 0.5, 0.0), (0.5, 0.5, 0.0), (0.5, -0.5, 0.0)])
# circle = polyline_rs.PolyStruct(1, [(-0.75, -0.25, 0.0), (-0.75, 0.25, 0.0), (0.75, 0.25, 0.0), (0.75, -0.25, 0.0)])

# # circle = [(-0.75, -0.25, 0.0), (-0.75, 0.25, 0.0), (0.75, 0.25, 0.0), (0.75, -0.25, 0.0)]
# # circle = [(0.75, 0.0, -0.5), (-0.75, 0.0, -0.5),]

# # result = polyline_rs.difference([square, circle])
# result = polyline_rs.difference(square, circle)

# # Perform difference_multiple
# positive_lines = [
#     polyline_rs.PolyStruct(1, [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0)])
# ]
# negative_lines = [
#     polyline_rs.PolyStruct(-1, [(0.5, 0.5, 0), (1.5, 0.5, 0), (1.5, 1.5, 0), (0.5, 1.5, 0)]),
#     polyline_rs.PolyStruct(-1, [(0.25, 0.75, 0), (1.5, 0.75, 0), (1.5, 1.5, 0), (0.25, 1.5, 0)])
# ]


# result = polyline_rs.difference_multiple(positive_lines, negative_lines)
# for poly in result:
#     print(f"Type: {poly.poly_type}, Polyline: {poly.polyline}")

# print("Union Result:", result)


# # Intersection
# result_intersection = polyline_rs.intersection([square, circle])
# print("Intersection Result:", result_intersection)

# # Difference
# result_difference = polyline_rs.difference(square, circle)
# print("Difference Result:", result_difference)
# for ind, r in enumerate(result):
#     result[ind] = (r[0], r[1], -r[2])

# Now I want to create an expression and save it. 
import geolipi.symbolic as gls
import woodie.symbolic as ws
from splitweaver.symbolic.base import CustomJSONEncoder
import json
from woodie.dag.mapper import convert_to_woodie_dag
import woodie.dag as wdag
from polyline_csg.parser import parse_csg_to_valid_polyset
from polyline_csg.polyset import polyset_to_csg

outer_circle_1 = ws.PolyLine2D(
    tuple([(0, 0, 0.5), (1, 0, 0.5), (1, 1, 0.5), (0, 1, 0.5)])
    )
# increase resolution
RESOLUTION = 1000
outer_circle_1 = ws.PolyLine2D(tuple([(x[0] * RESOLUTION, x[1] * RESOLUTION, x[2]) for x in outer_circle_1.args[0]]))

hole_1 = ws.PolyLine2D(
    tuple([
    (.2, .2, 0.5), (.8, .2, 0.5), (.8, .8, 0.5), (.2, .8, 0.5)
])
)
hole_1 = ws.PolyLine2D(tuple([(x[0] * RESOLUTION, x[1] * RESOLUTION, x[2]) for x in hole_1.args[0]]))
outer_circle_2 = ws.PolyLine2D(
    tuple([
    (0.5, 0, 0.5), (1.5, 0, 0.5), (1.5, 1, 0.5), (0.5, 1, 0.5)
]))
outer_circle_2 = ws.PolyLine2D(tuple([(x[0] * RESOLUTION, x[1] * RESOLUTION, x[2]) for x in outer_circle_2.args[0]]))
hole_2 = ws.PolyLine2D(
    tuple([
    (0.7, .2, 0.5), (1.3, .2, 0.5), (1.3, .8, 0.5), (0.7, .8, 0.5)
]))
hole_2 = ws.PolyLine2D(tuple([(x[0] * RESOLUTION, x[1] * RESOLUTION, x[2]) for x in hole_2.args[0]]))

expr = gls.Union(
    gls.Difference(outer_circle_1, hole_1),
    gls.Difference(outer_circle_2, hole_2)
)

result_polyset = parse_csg_to_valid_polyset(expr)
expr = polyset_to_csg(result_polyset)
print(expr)

def reduce_resolution(expr):
    if isinstance(expr, ws.PolyLine2D):
        return ws.PolyLine2D(tuple([(x[0] / RESOLUTION, x[1] / RESOLUTION, x[2]) for x in expr.args[0]]))
    else:
        return expr.__class__(*(reduce_resolution(arg) for arg in expr.args))
    
expr = reduce_resolution(expr)


plane_origin = (0, 0, 0)
plane_normal = (0, 1, 0)
pos_exprs = []
neg_exprs = []
        

expression = ws.LinkedHeightField3D(ws.Plane3D(plane_origin, plane_normal),
                                    ws.ApplyHeight(expr, (0.75,))
                                    )
expr_node_graph = convert_to_woodie_dag(expression)


pre_final = wdag.SetMaterial(expr=expr_node_graph.output_sockets['expr'], material=(2,),)
final = wdag.RegisterGeometry(expr=pre_final.output_sockets['expr'], name="base", bbox=(1, 1, 1))

main_dict = final.to_dict()

render =  wdag.RegisterState(expr=wdag.NamedGeometry(name="base"), state=(0,))

dict_2 = render.to_dict()
for key, value in dict_2.items():
    main_dict[key].extend(value)

positions = {}
for node_dict in main_dict['nodes']:
    positions[node_dict['id']] = {'x':0, 'y':0}

main_dict['positions'] = positions


graph_data = {'moduleList': {"main": main_dict}}

json_data = json.dumps(graph_data, cls=CustomJSONEncoder)

output_file = "/Users/adityaganeshan/projects/wood/wood_app/frontend/src/assets/projects/project.json"

with open(output_file, "w") as f:
    f.write(json_data)
    print(f"Saved to {output_file}")
