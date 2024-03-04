# read shapefiles in Curated 7, convert then to Networkx graphs, and merge nodes that are closer than a measure
import json
import math

import geopandas as gpd
import numpy as np
import pulp
from shapely.geometry import LineString, Point, Polygon
import momepy
import networkx as nx
import numpy
import os
import osmnx as ox
from collections import Counter
import gurobipy as gp
from gurobipy import GRB

import random
from deap import base, creator, tools, algorithms
from computational_geometry_functions import *



# Define the ILP problem evaluation function
def evaluate_ILP(individual, edges):
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


    # Evaluate fitness
    return total_weight, num_selected_vertices, penalty,



same = 0
different = 0
due = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

# create a dict to store pareto fronts for each floorplan
pareto_fronts = []

# create a dictionary named edgecount for number of edges for each node from 1 to 100
# for each node, count the number of edges, and update the dictionary
edgecount = {}
for i in range(101):
    edgecount[i] = 0

all_edges_in_dateset = 0
all_uncertain_edges_in_dateset = 0

uncertain_intersection = {}
for i in range(10):
    uncertain_intersection[i] = 0

floorplan_Bbox = gpd.read_file("Curated7/ICD_CalculatedV7.shp")
floorplan_border = gpd.read_file("Curated7/OneLetterV7.shp")
floorplan_skeleton = gpd.read_file("Curated7/SkelV7.shp")
floorplan_corner = gpd.read_file("Curated7/CornersV7.shp")

# clear all files in Geojson folder
folder = 'Geojson/'
for filename in os.listdir(folder):
    file_path = os.path.join(folder, filename)
    try:
        if os.path.isfile(file_path) or os.path.islink(file_path):
            os.unlink(file_path)

    except Exception as e:
        print('Failed to delete %s. Reason: %s' % (file_path, e))

counter = 0
Grammar = '6sectors'

for index, row in floorplan_Bbox.iterrows():
    counter += 1
    if counter > 4:
        break

    polygon = row['geometry']
    letter = row['distinct']
    print(letter)

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

                    # merge the nodes
                    # G = nx.contracted_edge(G, (u, v), self_loops=False)

                else:
                    due[i] = due[i] + 1
        i = i + 1
    # # PRINT THE NUMBER OF REMAINING NODES AND EDGES
    # print("number of nodes: " + str(G.number_of_nodes()))
    # print("number of edges: " + str(G.number_of_edges()))
    # check if the merged graph is the same as the original graph
    if nx.utils.graphs_equal(G, H):
        same += 1
    else:
        different += 1
    # points, lines = momepy.nx_to_gdf(G, lines=True, points=True)
    # points.make_valid()
    # points = points[['geometry', 'nodeID']]
    # points.to_file("Geojson/"+letter+"po.geojson", driver="GeoJSON", crs=32639)
    # count the number of edges for each node and update the edgecount dictionary
    for node in G.nodes():
        edgecount[G.degree(node)] += 1
        all_edges_in_dateset += G.degree(node)

    # calculate the bearing of all edges connecting to a node
    # create a dictionary that appendes all the uncertain edges for all nodes
    All_uncertain_edges = {}
    for Nodes in G.nodes():
        bearings = all_bearings_sorted(G, Nodes)
        if G.degree(Nodes) > 0:

            # calculate the nearest quarter for each bearing
            nearest_quarters = []
            for ang in bearings:
                nearest_quarters.append(nearest_quarter(ang))
            # print("the bearing : " + str(bearings))
            # print("the nearest_quarters : " + str(nearest_quarters))

            # check if the node is uncertain
            if check_if_the_node_is_uncertain(bearings):
                uncertain_intersection[G.degree(Nodes)] += 1
                u_edge, u_coord = check_edge_uncertainty(G, Nodes,Grammar=Grammar)
                all_uncertain_edges_in_dateset += len(u_edge.keys())
                # print(str(G.degree(Nodes))+"from which"+str(len(check_edge_uncertainty(G, Nodes).keys())))
                #print(u_edge)
                #print(u_coord)
                All_uncertain_edges[Nodes] = u_coord
    #print(All_uncertain_edges)

    # for keys, values in
    edges = uncertainty_weights_by_ID(G.edges(data=True), All_uncertain_edges)

    #create a dictionary that has this structure, ID of uncertain edge, weight of uncertain edge, and the vertices that are within distance from the edge
    # we will solve a linear integer programming problem to find the best vertices to allocate
    ILP = []
    for edge in edges:
        #check if the aedge has attriuibte weight
        if 'Weight' in edge[2]:
            if edge[2]['Weight'] > 1:
                # find the vertices within distance from the edge
                vertices, vertices_ID = find_vertices_within_distance_from_points(floorplan_Corners, gpd.GeoDataFrame(geometry=[LineString([edge[0], edge[1]])], crs=32639), 150)
                # extract the list of all IDs of vertices in the floorplan
                all_vertices_ID = list(floorplan_Corners['ID'])
                ILP.append({'ID': edge[2]['ID'], 'Weight': edge[2]['Weight'], 'vertices': vertices_ID, 'all_vertices': all_vertices_ID })
    # Define the ILP problem evaluation function and solve it using DEAP
    penalty_c = 50
    edges = ILP

    sum_of_weight = sum(edge['Weight'] for edge in edges)
    # Extract vertices from edges
    vertices = set()
    for edge in edges:
        vertices.update(edge['vertices'])

    # Create a dictionary to quickly find the index of vertices
    vertices_index = {vertex: index for index, vertex in enumerate(vertices)}

    # Create fitness and individual classes for the problem
    creator.create("FitnessMulti", base.Fitness, weights=(1.0, -1.0, 500))  # Minimize both objectives and the penalty
    creator.create("Individual", list, fitness=creator.FitnessMulti)

    # Create toolbox for the problem
    toolbox = base.Toolbox()
    toolbox.register("attr_bool", random.randint, 0, 1)
    toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, len(vertices))
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    toolbox.register("evaluate", evaluate_ILP, edges=edges)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
    toolbox.register("select", tools.selNSGA2)

    # Define the number of generations and population size
    NGEN = 5000
    MU = 100

    # Initialize the population
    pop = toolbox.population(n=MU)

    # Run the evolution
    try:
        algorithms.eaMuPlusLambda(pop, toolbox, mu=MU, lambda_=MU, cxpb=0.7, mutpb=0.2, ngen=NGEN, verbose=False)
    except:
        continue

    # Extract the Pareto front
    pareto_front = tools.sortNondominated(pop, len(pop), first_front_only=True)[0]
    for ind in pareto_front:
        total_weight, num_selected_vertices, penalties = evaluate_ILP(ind, edges)

        # update the pareto_fronts dictionary
        pareto_fronts.append({'letter': letter, 'total_weight': total_weight, 'num_selected_vertices': num_selected_vertices,
                              'penalties': penalties, 'uncertain_edges': len(edges), 'sum_of_weight': sum_of_weight})

    # write down a single json containg letter, number of require vertices, number of selected vertices and the total weight
    # for each pareto fronts that has nonnegative total weight
    unique_pareto_fronts = []
    for pareto_front in pareto_fronts:
        if pareto_front['total_weight'] > 0:
            if pareto_front not in unique_pareto_fronts:
                unique_pareto_fronts.append(pareto_front)
    # print(unique_pareto_fronts)
    with open('pareto_fronts.json', 'w') as f:
        for element in unique_pareto_fronts:
            json.dump(element, f)
            f.write('\n')






