import json
import math
import os

import geopandas as gpd
import numpy as np
import pulp
from shapely.geometry import LineString, Point, Polygon
import momepy
import networkx as nx
import osmnx as ox
from collections import Counter
import gurobipy as gp
from gurobipy import GRB
import random
from deap import base, creator, tools, algorithms
from multiprocessing import Pool, cpu_count
from functools import partial

def edge_bearing(edge):
    tuple1 = ox.projection.project_geometry(Point(edge[0]), crs=32639, to_latlong=True)[0]
    tuple2 = ox.projection.project_geometry(Point(edge[1]), crs=32639, to_latlong=True)[0]
    bearing = ox.bearing.calculate_bearing(tuple1.y, tuple1.x, tuple2.y, tuple2.x)
    return bearing

def all_bearings_sorted(G, Node):
    bearings = []
    if G.degree(Node) > 0:
        Edges = G.edges(Node, data=True)
        for u, v, attr in Edges:
            if u == Node:
                bearing = edge_bearing((u, v))
            elif v == Node:
                bearing = edge_bearing((v, u))
            bearings.append(bearing)
        for i in range(len(bearings)):
            bearings[i] = bearings[i] % 360
        bearings.sort()
    return bearings

def nearest_quarter(ang):
    ang = ang % 360
    return (ang // 90) * 90

def calculate_edge_angles(bearings_sorted):
    num_edges = len(bearings_sorted)
    angles = []
    for i in range(num_edges):
        next_index = (i + 1) % num_edges
        angle = (bearings_sorted[next_index] - bearings_sorted[i] + 360) % 360
        angles.append(angle)
    return angles

def check_if_the_node_is_uncertain(bearings_sorted):
    angles_between_edges = calculate_edge_angles(bearings_sorted)
    for bearingss in bearings_sorted:
        angles_relative = [- bearingss + angle for angle in bearings_sorted]
        angles_relative = [angle % 360 for angle in angles_relative]
        nearest_quarters = [nearest_quarter(angle) for angle in angles_relative]
        if max(Counter(nearest_quarters).values()) >= 2 and len(bearings_sorted) > 1:
            return True
    return False

def check_if_two_outgoing_are_uncertain(out1, out2):
    if nearest_quarter(out1 % 360) == nearest_quarter(out2 % 360):
        return True
    return False

def check_incoming_uncertainty(G, Node, edge):
    bearings = all_bearings_sorted(G, Node)
    incoming_bearing = (edge_bearing(edge) + 180) % 360
    for ic in range(len(bearings)):
        bearings[ic] = (bearings[ic] - incoming_bearing) % 360
    bearings = list(set(bearings))
    bearings.remove(0)
    bearings.sort()
    for i in range(len(bearings)):
        if i < len(bearings)-1:
            if check_if_two_outgoing_are_uncertain(bearings[i], bearings[i+1]):
                return True
        elif i == len(bearings)-1:
            if check_if_two_outgoing_are_uncertain(bearings[i], bearings[0]):
                return True
    return False

def check_edge_uncertainty(G, Node):
    bearings = all_bearings_sorted(G, Node)
    edge_uncertainty = {}
    edge_uncertainty_by_coordinates = {}
    for u, v, attr in G.edges(Node, data=True):
        incoming_bearing = edge_bearing((v, u))
        for u1, v1, attr1 in G.edges(Node, data=True):
            outgoing_bearing1 = edge_bearing((u1, v1))
            outgoing_bearing1_corrected = (outgoing_bearing1 - incoming_bearing) % 360
            for u2, v2, attr2 in G.edges(Node, data=True):
                outgoing_bearing2 = edge_bearing((u2, v2))
                outgoing_bearing2_corrected = (outgoing_bearing2 - incoming_bearing) % 360
                incoming_bearing_rounded = round(incoming_bearing, 2)
                incoming_bearing_rounded = (incoming_bearing_rounded + 180) % 360
                outgoing_bearing1_rounded = round(outgoing_bearing1, 2)
                outgoing_bearing2_rounded = round(outgoing_bearing2, 2)
                u_x = round(u[0], 0)
                u_y = round(u[1], 0)
                v_x = round(v[0], 0)
                v_y = round(v[1], 0)
                u1_x = round(u1[0], 0)
                u1_y = round(u1[1], 0)
                v1_x = round(v1[0], 0)
                v1_y = round(v1[1], 0)
                u2_x = round(u2[0], 0)
                u2_y = round(u2[1], 0)
                v2_x = round(v2[0], 0)
                v2_y = round(v2[1], 0)
                if check_if_two_outgoing_are_uncertain(outgoing_bearing1_corrected, outgoing_bearing2_corrected) and not np.allclose(outgoing_bearing1, outgoing_bearing2, atol=0.1):
                    if incoming_bearing_rounded in edge_uncertainty:
                        if (outgoing_bearing1_rounded, outgoing_bearing2_rounded) not in edge_uncertainty[incoming_bearing_rounded] and (outgoing_bearing2_rounded, outgoing_bearing1_rounded) not in edge_uncertainty[incoming_bearing_rounded]:
                            edge_uncertainty[incoming_bearing_rounded].append((outgoing_bearing1_rounded, outgoing_bearing2_rounded))
                            edge_uncertainty_by_coordinates[((v_x, v_y), (u_x, u_y))].append((((u1_x, u1_y), (v1_x, v1_y)), ((u2_x, u2_y), (v2_x, v2_y))))
                    else:
                        edge_uncertainty[incoming_bearing_rounded] = [(outgoing_bearing1_rounded, outgoing_bearing2_rounded)]
                        edge_uncertainty_by_coordinates[((v_x, v_y), (u_x, u_y))] = [(((u1_x, u1_y), (v1_x, v1_y)), ((u2_x, u2_y), (v2_x, v2_y)))]
    return edge_uncertainty, edge_uncertainty_by_coordinates

def uncertainty_weights_by_ID(edges, uncertain_edges):
    for edge in edges:
        u = edge[0]
        v = edge[1]
        u_x = round(u[0], 0)
        u_y = round(u[1], 0)
        v_x = round(v[0], 0)
        v_y = round(v[1], 0)
        for key, value in uncertain_edges.items():
            for key1, value1 in value.items():
                if ((u_x, u_y), (v_x, v_y)) in value1[0] or ((v_x, v_y), (u_x, u_y)) in value1[0]:
                    if 'Weight' in edge[2]:
                        edge[2]['Weight'] += 1
                    else:
                        edge[2]['Weight'] = 1
    return edges

def find_vertices_within_distance(polygon_A, point_B, distance):
    polygon_A = polygon_A.to_crs(point_B.crs)
    point_B_geometry = point_B.geometry.iloc[0]
    buffered_point_B = point_B_geometry.buffer(distance)
    vertices_within_distance = []
    vertices_within_distance_ID = []
    for index, row in polygon_A.iterrows():
        vertices = list(row['geometry'].exterior.coords)
        for vertex in vertices:
            if Point(vertex).within(buffered_point_B):
                vertices_within_distance.append(vertex)
    return vertices_within_distance

def find_vertices_within_distance_from_points(Corners, geomB, distance):
    Corners = Corners.to_crs(geomB.crs)
    point_B_geometry = geomB.geometry.iloc[0]
    buffered_point_B = point_B_geometry.buffer(distance, cap_style=3).buffer(1)
    vertices_within_distance = []
    vertices_within_distance_ID = []
    for index, row in Corners.iterrows():
        if Point(row['geometry']).within(buffered_point_B):
            vertices_within_distance.append(row['geometry'])
            vertices_within_distance_ID.append(row['ID'])
    return vertices_within_distance, vertices_within_distance_ID

def ILP_solver(edges):
    vertex_vars = {vertex: pulp.LpVariable(f'X_{vertex}', cat=pulp.LpBinary) for edge in edges for vertex in edge['vertices']}
    prob = pulp.LpProblem("Vertex_Selection", pulp.LpMaximize)
    prob += pulp.lpSum([edge['Weight'] * vertex_vars[vertex] for edge in edges for vertex in edge['vertices']]), "Total_Weight"
    prob += pulp.lpSum([vertex_vars[vertex] for edge in edges for vertex in edge['vertices']]), "Total_Vertices"
    for edge in edges:
        prob += pulp.lpSum([vertex_vars[vertex] for vertex in edge['vertices']]) == 1, f"Edge_{edge['ID']}_Selection"
    prob.solve()
    selected_vertices = [v.name.split('_')[1] for v in prob.variables() if v.varValue == 1]
    total_weight = sum(edges[0]['Weight'] for edge in edges for vertex in edge['vertices'] if vertex in selected_vertices)
    return selected_vertices

def ILP_solver_Gorubi(edges):
    model = gp.Model("Vertex_Selection")
    vertex_vars = {}
    for edge in edges:
        for vertex in edge['vertices']:
            vertex_vars[vertex] = model.addVar(vtype=GRB.BINARY, name=f'X_{vertex}')
    model.update()
    total_weight_expr = gp.LinExpr()
    for edge in edges:
        for vertex in edge['vertices']:
            total_weight_expr += edge['Weight'] * vertex_vars[vertex]
    model.setObjective(total_weight_expr, GRB.MAXIMIZE)
    for edge in edges:
        model.addConstr(gp.quicksum(vertex_vars[vertex] for vertex in edge['vertices']) == 1, f"Edge_{edge['ID']}_Selection")
    model.optimize()
    selected_vertices = [var.varName.split('_')[1] for var in model.getVars() if var.x > 0.5]
    total_weight = sum(edge['Weight'] for edge in edges for vertex in edge['vertices'] if vertex in selected_vertices)
    return selected_vertices

def evaluate_ILP(individual, edges):
    sum_of_weight = sum(edge['Weight'] for edge in edges)
    vertices = set()
    for edge in edges:
        vertices.update(edge['vertices'])
    vertices_index = {vertex: index for index, vertex in enumerate(vertices)}
    selected_vertices = [vertex for vertex, selected in zip(vertices, individual) if selected]
    total_weight = sum(edge['Weight'] for edge in edges for vertex in edge['vertices'] if vertex in selected_vertices)
    num_selected_vertices = sum(individual)
    penalty = 0
    for edge in edges:
        selected_edge_vertices = [vertex for vertex in edge['vertices'] if individual[vertices_index[vertex]] == 1]
        if len(selected_edge_vertices) > 1:
            total_weight += - penalty_c
            penalty = 1
    return total_weight, num_selected_vertices, penalty,

def process_floorplan(row):
    polygon = row['geometry']
    letter = row['distinct']

    floorplan_intersecting = floorplan_skeleton[floorplan_skeleton.intersects(polygon)]
    floorplan_border_intersecting = floorplan_border[floorplan_border.intersects(polygon)]
    floorplan_Corners = floorplan_corner[floorplan_corner.intersects(polygon)]

    G = momepy.gdf_to_nx(floorplan_intersecting, approach="primal", length="length")
    H = G

    i = 0

    while i < 10:  # Merge nodes that are closer than some meter
        for u, v, key, attr in list(G.edges(keys=True, data=True)):
            if attr['length'] < 100:
                # check if the edge still exists in the graph
                if G.has_edge(u, v):
                    if G.degree(u) >= G.degree(v):
                        # merge the nodes
                        G = nx.contracted_nodes(G, u, v, self_loops=False)
                    else:
                        # merge the nodes
                        G = nx.contracted_nodes(G, v, u, self_loops=False)

        i = i + 1

    if nx.utils.graphs_equal(G, H):
        same += 1
    else:
        different += 1

    for node in G.nodes():
        edgecount[G.degree(node)] += 1
        all_edges_in_dateset += G.degree(node)

    All_uncertain_edges = {}
    for Nodes in G.nodes():
        bearings = all_bearings_sorted(G, Nodes)
        if G.degree(Nodes) > 0:
            nearest_quarters = []
            for ang in bearings:
                nearest_quarters.append(nearest_quarter(ang))

            if check_if_the_node_is_uncertain(bearings):
                uncertain_intersection[G.degree(Nodes)] += 1
                u_edge, u_coord = check_edge_uncertainty(G, Nodes)
                all_uncertain_edges_in_dateset += len(u_edge.keys())
                All_uncertain_edges[Nodes] = u_coord

    edges = uncertainty_weights_by_ID(G.edges(data=True), All_uncertain_edges)

    ILP = []

    for edge in edges:
        if 'Weight' in edge[2]:
            if edge[2]['Weight'] > 1:
                vertices, vertices_ID = find_vertices_within_distance_from_points(floorplan_Corners, gpd.GeoDataFrame(geometry=[LineString([edge[0], edge[1]])], crs=32639), 150)
                ILP.append({'ID': edge[2]['ID'], 'Weight': edge[2]['Weight'], 'vertices': vertices_ID})

    pareto_fronts = []

    for index, row in floorplan_Bbox.iterrows():
        counter += 1
        if counter > 20:
            break

        result = process_floorplan(index, row)
        if result is not None:
            pareto_fronts.append(result)


if __name__ == '__main__':
    # Load your shapefiles and prepare your data here
    floorplan_Bbox = gpd.read_file("Curated7/ICD_CalculatedV7.shp")
    floorplan_border = gpd.read_file("Curated7/OneLetterV7.shp")
    floorplan_skeleton = gpd.read_file("Curated7/SkelV7.shp")
    floorplan_corner = gpd.read_file("Curated7/CornersV7.shp")

    row = [row for index, row in floorplan_Bbox.iterrows()]

    # Clear all files in Geojson folder
    folder = 'Geojson/'
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))

    counter = 0
    same = 0
    different = 0
    due = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    all_edges_in_dateset = 0
    all_uncertain_edges_in_dateset = 0
    uncertain_intersection = {}
    for i in range(10):
        uncertain_intersection[i] = 0
    edgecount = {i: 0 for i in range(101)}

    # Define the number of processes
    num_processes = os.cpu_count() - 1  # You can adjust this number

    # Pool of workers for parallel execution
    with Pool(processes=num_processes) as pool:
        # Use partial to create a function with some arguments fixed
        partial_process_floorplan = partial(process_floorplan,)

        # Iterate over the shapefile rows and map them to the pool
        # Collect any results you need from the map function
        results = pool.map(partial_process_floorplan, enumerate(floorplan_Bbox.iterrows()))
