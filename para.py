import json
import math
import multiprocessing

import geopandas as gpd
import numpy as np
from shapely.geometry import LineString, Point, Polygon
import momepy
import networkx as nx
import numpy
import os
import osmnx as ox
from collections import Counter


import random
from deap import base, creator, tools, algorithms
from computational_geometry_functions import *


def evaluate_ILP(individual, edges,vertices,penalty_c=50):
    try:
        # Extract the selected vertices from the individual
        Xc = [vertex for vertex, is_selected in zip(vertices, individual) if is_selected]

        # see if a vertex is among all_vertices of an edge (if vertex is in the buffer of an edge)
        Yodc = {}
        for edgesss in edges:
            for vertex in edgesss['all_vertices']:
                if vertex in Xc and vertex in edgesss['vertices']:
                    Yodc[(edgesss['ID'], vertex)] = 1
                else:
                    Yodc[(edgesss['ID'], vertex)] = 0

        # Calculate the total weight and the number of selected vertices
        total_weight = sum(edge['Weight'] * Yodc[(edge['ID'], vertex)] for edge in edges for vertex in edge['all_vertices'])
        num_selected_vertices = sum(individual)

        # Penalty for violating the constraint: For every edge, only one vertex must be selected
        penalty = 0
        for edge in edges:
            if sum(Yodc[(edge['ID'], vertex)] for vertex in edge['vertices']) > 1:
                penalty = (sum(Yodc[(edge['ID'], vertex)] for vertex in edge['vertices'])) - 1
                total_weight += - penalty_c * penalty
        # penalty for violating the constraint: the Yodc can be 1 if only corresponding vertex is both selected and edge[
        # 'vertices'], otherwise it should be 0
        for edge in edges:
            for vertex in edge['vertices']:
                if Yodc[(edge['ID'], vertex)] == 1 and vertex not in Xc:
                    penalty += - 1
                    total_weight += - penalty_c * penalty

        # return the total_weight, num_selected_vertices, penalty, and the ID of selected vertices
        return total_weight, num_selected_vertices, penalty,
    except Exception as e:
        import traceback

        traceback.print_exc()
        raise

def process_floorplan(index, row):
    polygon = row['geometry']
    letter = row['distinct']
    print(letter)
    floorplan_skeleton = gpd.read_file("data_realworld/shpa/chadstone_graph_utm.shp")
    floorplan_corner = gpd.read_file("data_realworld/shpa/chadstone_Corners_utm.shp")

    edgecount = {}
    for i in range(101):
        edgecount[i] = 0

    uncertain_intersection = {}
    for i in range(10):
        uncertain_intersection[i] = 0


    floorplan_intersecting = floorplan_skeleton[floorplan_skeleton.intersects(polygon)]
    # floorplan_border_intersecting = floorplan_border[floorplan_border.intersects(polygon)]
    floorplan_Corners = floorplan_corner[floorplan_corner.intersects(polygon)]

    G = momepy.gdf_to_nx(floorplan_intersecting, approach="primal", length="length")
    print(G)
    H = G

    i = 0

    while i < 10:  # Merge nodes that are closer than some meter
        for u, v, key, attr in list(G.edges(keys=True, data=True)):
            if attr['length'] < 10:
                # check if the edge still exists in the graph
                if G.has_edge(u, v):
                    if G.degree(u) >= G.degree(v):
                        # merge the nodes
                        G = nx.contracted_nodes(G, u, v, self_loops=False)
                    else:
                        # merge the nodes
                        G = nx.contracted_nodes(G, v, u, self_loops=False)
        i = i + 1
    print(G)

    for node in G.nodes():
        edgecount[G.degree(node)] += 1

    All_uncertain_edges = {}
    for Nodes in G.nodes():
        bearings = all_bearings_sorted(G, Nodes)
        if G.degree(Nodes) > 0:
            nearest_quarters = []
            for ang in bearings:
                nearest_quarters.append(nearest_quarter(ang))

            if check_if_the_node_is_uncertain(bearings):
                uncertain_intersection[G.degree(Nodes)] += 1
                u_edge, u_coord = check_edge_uncertainty(G, Nodes, Grammar='4sectors')
                All_uncertain_edges[Nodes] = u_coord

    edges = uncertainty_weights_by_ID(G.edges(data=True), All_uncertain_edges)

    ILP = []
    for edge in edges:
        if 'Weight' in edge[2]:
            if edge[2]['Weight'] > 1:
                vertices, vertices_ID = find_vertices_within_distance_from_points(floorplan_Corners, gpd.GeoDataFrame(
                    geometry=[LineString([edge[0], edge[1]])], crs=32639), 10)
                all_vertices_ID = list(floorplan_Corners['ID'])
                ILP.append({'ID': edge[2]['ID'], 'Weight': edge[2]['Weight'], 'vertices': vertices_ID,
                            'all_vertices': all_vertices_ID})

    penalty_c = 50
    edges = ILP

    sum_of_weight = sum(edge['Weight'] for edge in edges)
    vertices = set()
    for edge in edges:
        vertices.update(edge['vertices'])

    vertices_index = {vertex: index for index, vertex in enumerate(vertices)}

    creator.create("FitnessMulti", base.Fitness, weights=(1.0, -1.0, 500))
    creator.create("Individual", list, fitness=creator.FitnessMulti)

    toolbox = base.Toolbox()
    toolbox.register("attr_bool", random.randint, 0, 1)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, len(vertices))
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    toolbox.register("evaluate", evaluate_ILP, edges=edges, vertices=list(vertices), penalty_c=penalty_c)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
    toolbox.register("select", tools.selNSGA2)

    NGEN = 5000
    MU = 500

    pop = toolbox.population(n=MU)

    try:
        algorithms.eaMuPlusLambda(pop, toolbox, mu=MU, lambda_=MU, cxpb=0.7, mutpb=0.2, ngen=NGEN, verbose=False)
    except:
        print("An exception occurred")
        return []

    pareto_front = tools.sortNondominated(pop, len(pop), first_front_only=True)[0]
    pareto_fronts = []
    for ind in pareto_front:
        selected_vertex_ids = [list(vertices)[i] for i, is_selected in enumerate(ind) if is_selected]
        total_weight, num_selected_vertices, penalties = evaluate_ILP(ind, edges, list(vertices), penalty_c=penalty_c)
        pareto_fronts.append(
            {'letter': letter, 'total_weight': total_weight, 'num_selected_vertices': num_selected_vertices,
             'penalties': penalties, 'uncertain_edges': len(edges), 'sum_of_weight': sum_of_weight,
             'selected_vertex_ids': selected_vertex_ids})
    print(pareto_fronts)
    with open('pareto_fronts_Chadstone.json', 'w') as f:
        for element in pareto_fronts:
            json.dump(element, f)
            f.write('\n')
    unique_pareto_fronts = []
    for pareto_front in pareto_fronts:
        if pareto_front['total_weight'] > 0:
            if pareto_front not in unique_pareto_fronts:
                unique_pareto_fronts.append(pareto_front)

    return unique_pareto_fronts


def parallel_process():
    pool = multiprocessing.Pool()
    results = []

    for index, row in floorplan_Bbox.iterrows():
        result = pool.apply_async(process_floorplan, args=(index, row))
        results.append(result)

    pool.close()
    pool.join()

    pareto_fronts = []
    for result in results:
        pareto_fronts.extend(result.get())

    with open('pareto_fronts_Falcon.json', 'w') as f:
        for element in pareto_fronts:
            json.dump(element, f)
            f.write('\n')


if __name__ == '__main__':
    same = 0
    different = 0
    due = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    penalty_c = 50

    pareto_fronts = []

    edgecount = {}
    for i in range(101):
        edgecount[i] = 0

    all_edges_in_dateset = 0
    all_uncertain_edges_in_dateset = 0

    uncertain_intersection = {}
    for i in range(10):
        uncertain_intersection[i] = 0

    # floorplan_Bbox = gpd.read_file("data_realworld/shpa/Falcon_BBox.shp")
    # floorplan_border = gpd.read_file("data_realworld/shpa/Falcon_OneArea.shp")
    # floorplan_skeleton = gpd.read_file("data_realworld/shpa/falcon_graph.shp")
    # floorplan_corner = gpd.read_file("data_realworld/shpa/Falcon_Corners.shp")

    floorplan_Bbox = gpd.read_file("data_realworld/shpa/chadstone_BBox_utm.shp")
    floorplan_border = gpd.read_file("data_realworld/shpa/chadstone_area_utm.shp")
    floorplan_skeleton = gpd.read_file("data_realworld/shpa/chadstone_graph_utm.shp")
    floorplan_corner = gpd.read_file("data_realworld/shpa/chadstone_Corners_utm.shp")
    parallel_process()