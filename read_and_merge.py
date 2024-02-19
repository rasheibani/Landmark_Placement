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

from multiprocessing import Pool, cpu_count




def edge_bearing(edge):
    # calculate the bearing of an edge
    tuple1 = ox.projection.project_geometry(Point(edge[0]), crs=32639, to_latlong=True)[0]
    tuple2 = ox.projection.project_geometry(Point(edge[1]), crs=32639, to_latlong=True)[0]
    bearing = ox.bearing.calculate_bearing(tuple1.y, tuple1.x, tuple2.y, tuple2.x)
    return bearing


def all_bearings_sorted(G, Node):
    # calculate the bearing of all edges connecting to a node
    bearings = []
    if G.degree(Node) > 0:
        # print(Nodes)
        Edges = G.edges(Node, data=True)
        # for edge in Edges:
        # print(edge)
        for u, v, attr in Edges:
            if u == Node:
                bearing = edge_bearing((u, v))
            elif v == Node:
                bearing = edge_bearing((v, u))
            bearings.append(bearing)
        # Sort the bearings
        for i in range(len(bearings)):
            bearings[i] = bearings[i] % 360

        bearings.sort()
    return bearings


def nearest_quarter(ang):
    # calculate the nearest quarter for an angle
    ang = ang % 360
    return (ang // 90) * 90


def calculate_edge_angles(bearings_sorted):
    num_edges = len(bearings_sorted)
    angles = []

    # Calculate angles between consecutive edges
    for i in range(num_edges):
        # Calculate the index of the next edge, considering circularity
        next_index = (i + 1) % num_edges

        # Calculate the angular difference between consecutive edges
        angle = (bearings_sorted[next_index] - bearings_sorted[i] + 360) % 360
        angles.append(angle)

    return angles


def check_if_the_node_is_uncertain(bearings_sorted):
    # check if the node is uncertain
    # print("the bearing : " + str(bearings_sorted))
    angles_between_edges = calculate_edge_angles(bearings_sorted)
    # print("the angles_between_edges : " + str(angles_between_edges))
    for bearingss in bearings_sorted:
        # subtract all elements of the list from the bearingss
        angles_relative = [- bearingss + angle for angle in bearings_sorted]
        # mod 360 all elements of the list
        angles_relative = [angle % 360 for angle in angles_relative]
        # claculate the quarter of angles_relative
        nearest_quarters = [nearest_quarter(angle) for angle in angles_relative]
        # check if there is a quarter that is repeated 2 times or more
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

    # remove duplicate and 0 from the list
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
    #print("Node: " + str(Node) + " has " + str(G.degree(Node)) + " edges")
    # return a dictionary of edges that are uncertain in this way: {incoimg_edge: [(outgoing_edge1,utgoing_edge2),...],...}
    bearings = all_bearings_sorted(G, Node)
    #print("the bearings : " + str(bearings))
    edge_uncertainty = {}
    edge_uncertainty_by_coordinates = {}

    # iterate over all edges as incoming edges
    for u, v, attr in G.edges(Node, data=True):
        # calculate the bearing of the incoming edge
        incoming_bearing = edge_bearing((v, u))

        # iterate over all edges as outcoming edges 1
        for u1, v1, attr1 in G.edges(Node, data=True):
            outgoing_bearing1 = edge_bearing((u1, v1))
            outgoing_bearing1_corrected = (outgoing_bearing1 - incoming_bearing) % 360


            # iterate over all edges as outcoming edges 2
            for u2, v2, attr2 in G.edges(Node, data=True):
                outgoing_bearing2 = edge_bearing((u2, v2))
                outgoing_bearing2_corrected = (outgoing_bearing2 - incoming_bearing) % 360

                incoming_bearing_rounded = round(incoming_bearing, 2)
                incoming_bearing_rounded = (incoming_bearing_rounded + 180) % 360
                outgoing_bearing1_rounded = round(outgoing_bearing1, 2)
                outgoing_bearing2_rounded = round(outgoing_bearing2, 2)

                #print the first element of tuple u

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


                # check if the two outcoming edges are uncertain
                if check_if_two_outgoing_are_uncertain(outgoing_bearing1_corrected, outgoing_bearing2_corrected) and not np.allclose(outgoing_bearing1, outgoing_bearing2, atol=0.1):
                    if incoming_bearing_rounded in edge_uncertainty:
                        if (outgoing_bearing1_rounded, outgoing_bearing2_rounded) not in edge_uncertainty[incoming_bearing_rounded] and (outgoing_bearing2_rounded, outgoing_bearing1_rounded) not in edge_uncertainty[incoming_bearing_rounded]:
                            edge_uncertainty[incoming_bearing_rounded].append((outgoing_bearing1_rounded, outgoing_bearing2_rounded))

                            edge_uncertainty_by_coordinates[((v_x, v_y), (u_x, u_y))].append(
                                (((u1_x, u1_y), (v1_x, v1_y)), ((u2_x, u2_y), (v2_x, v2_y))))
                    else:
                        edge_uncertainty[incoming_bearing_rounded] = [(outgoing_bearing1_rounded, outgoing_bearing2_rounded)]

                        edge_uncertainty_by_coordinates[((v_x, v_y), (u_x, u_y))] = [
                            (((u1_x, u1_y), (v1_x, v1_y)), ((u2_x, u2_y), (v2_x, v2_y)))]



    # look for similar elements with two decimal digits and remove them

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
                    # check if attribute Weight exists
                    if 'Weight' in edge[2]:
                        edge[2]['Weight'] += 1
                    else:
                        edge[2]['Weight'] = 1

    return edges


def find_vertices_within_distance(polygon_A, point_B, distance):

    # Ensure both datasets have the same CRS (Coordinate Reference System)
    polygon_A = polygon_A.to_crs(point_B.crs)

    # Get the geometry of point B
    point_B_geometry = point_B.geometry.iloc[0]

    # Create a buffer around point B with the specified radius
    buffered_point_B = point_B_geometry.buffer(distance)

    # Filter the vertices of polygon A that fall within the buffered area around point B
    vertices_within_distance = []
    vertices_within_distance_ID = []
    for index, row in polygon_A.iterrows():
        vertices = list(row['geometry'].exterior.coords)
        for vertex in vertices:
            if Point(vertex).within(buffered_point_B):
                vertices_within_distance.append(vertex)

    return vertices_within_distance

def find_vertices_within_distance_from_points(Corners, geomB, distance):

    # Ensure both datasets have the same CRS (Coordinate Reference System)
    Corners = Corners.to_crs(geomB.crs)

    # Get the geometry of point B
    point_B_geometry = geomB.geometry.iloc[0]

    # Create a buffer around point B with the specified radius
    buffered_point_B = point_B_geometry.buffer(distance, cap_style=3).buffer(1)

    # Filter the vertices of polygon A that fall within the buffered area around point B
    vertices_within_distance = []
    vertices_within_distance_ID = []
    for index, row in Corners.iterrows():
        if Point(row['geometry']).within(buffered_point_B):
            vertices_within_distance.append(row['geometry'])
            vertices_within_distance_ID.append(row['ID'])

    return vertices_within_distance, vertices_within_distance_ID
def ILP_solver(edges):


    # Create a binary variable for each vertex
    vertex_vars = {vertex: pulp.LpVariable(f'X_{vertex}', cat=pulp.LpBinary) for edge in edges for vertex in
                   edge['vertices']}

    # Create the problem instance
    prob = pulp.LpProblem("Vertex_Selection", pulp.LpMaximize)

    # Objective function 1
    # maximise the total weight of selected vertices
    prob += pulp.lpSum(
        [edge['Weight'] * vertex_vars[vertex] for edge in edges for vertex in edge['vertices']]), "Total_Weight"
    # Objective function 2
     # minimise the number of selected vertices
    prob += pulp.lpSum([vertex_vars[vertex] for edge in edges for vertex in edge['vertices']]), "Total_Vertices"


    # Constraints: Each edgeID should have one selected vertex
    for edge in edges:
        prob += pulp.lpSum([vertex_vars[vertex] for vertex in edge['vertices']]) == 1, f"Edge_{edge['ID']}_Selection"

    # Solve the ILP problem
    prob.solve()

    selected_vertices = [v.name.split('_')[1] for v in prob.variables() if v.varValue == 1]
    # print("Selected Vertices:", selected_vertices)

    # Print the total weight of selected vertices
    total_weight = sum(
        edges[0]['Weight'] for edge in edges for vertex in edge['vertices'] if vertex in selected_vertices)
    return selected_vertices

def ILP_solver_Gorubi(edges):
    # Create a new model
    model = gp.Model("Vertex_Selection")

    # Create binary variables for each vertex
    vertex_vars = {}
    for edge in edges:
        for vertex in edge['vertices']:
            vertex_vars[vertex] = model.addVar(vtype=GRB.BINARY, name=f'X_{vertex}')

    # Update model to integrate new variables
    model.update()

    # Set up objective function 1: Maximize total weight
    total_weight_expr = gp.LinExpr()
    for edge in edges:
        for vertex in edge['vertices']:
            total_weight_expr += edge['Weight'] * vertex_vars[vertex]
    model.setObjective(total_weight_expr, GRB.MAXIMIZE)

    # # Set up objective function 2: Minimize overall selected vertices
    # model.setObjective(gp.quicksum(vertex_vars[vertex] for vertex in vertex_vars.values()), GRB.MINIMIZE)

    # Add constraints: Each edgeID should have one selected vertex
    for edge in edges:
        model.addConstr(gp.quicksum(vertex_vars[vertex] for vertex in edge['vertices']) == 1, f"Edge_{edge['ID']}_Selection")

    # Optimize the model
    model.optimize()

    # Retrieve selected vertices
    selected_vertices = [var.varName.split('_')[1] for var in model.getVars() if var.x > 0.5]

    # Print selected vertices
    print("Selected Vertices:", selected_vertices)

    # Print total weight of selected vertices
    total_weight = sum(edge['Weight'] for edge in edges for vertex in edge['vertices'] if vertex in selected_vertices)
    print("Total Weight of Selected Vertices:", total_weight)

    return selected_vertices


# Define the ILP problem evaluation function
def evaluate_ILP(individual, edges):
    sum_of_weight = sum(edge['Weight'] for edge in edges)
    # Extract vertices from edges
    vertices = set()
    for edge in edges:
        vertices.update(edge['vertices'])
    # Create a dictionary to quickly find the index of vertices
    vertices_index = {vertex: index for index, vertex in enumerate(vertices)}

    selected_vertices = [vertex for vertex, selected in zip(vertices, individual) if selected]
    total_weight = sum(edge['Weight'] for edge in edges for vertex in edge['vertices'] if vertex in selected_vertices)
    num_selected_vertices = sum(individual)

    # Penalty for violating the constraint: For every edge, only one vertex must be selected
    penalty = 0
    for edge in edges:
        selected_edge_vertices = [vertex for vertex in edge['vertices'] if individual[vertices_index[vertex]] == 1]

        if len(selected_edge_vertices) > 1:
            total_weight += - penalty_c
            penalty = 1

    # Evaluate fitness
    return total_weight, num_selected_vertices, penalty

def process (row):


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



    # clear all files in Geojson folder
    folder = 'Geojson/'
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)

        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))


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


                # check if the node is uncertain
                if check_if_the_node_is_uncertain(bearings):
                    uncertain_intersection[G.degree(Nodes)] += 1
                    u_edge, u_coord = check_edge_uncertainty(G, Nodes)
                    all_uncertain_edges_in_dateset += len(u_edge.keys())

                    All_uncertain_edges[Nodes] = u_coord

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
                    ILP.append({'ID': edge[2]['ID'], 'Weight': edge[2]['Weight'], 'vertices': vertices_ID})


        edges = ILP

        sum_of_weight = sum(edge['Weight'] for edge in edges)
        # Extract vertices from edges
        vertices = set()
        for edge in edges:
            vertices.update(edge['vertices'])



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
        NGEN = 1000
        MU = 200

        # Initialize the population
        pop = toolbox.population(n=MU)

        # Run the evolution
        algorithms.eaMuPlusLambda(pop, toolbox, mu=MU, lambda_=MU, cxpb=0.7, mutpb=0.2, ngen=NGEN, verbose=False)

        # Extract the Pareto front
        pareto_front = tools.sortNondominated(pop, len(pop), first_front_only=True)[0]

        for ind in pareto_front:
            total_weight, num_selected_vertices, penalties = evaluate_ILP(ind, edges)

            # update the pareto_fronts dictionary
            pareto_fronts.append(
                {'letter': letter, 'total_weight': total_weight, 'num_selected_vertices': num_selected_vertices,
                 'penalties': penalties, 'uncertain_edges': len(edges), 'sum_of_weight': sum_of_weight})

    print(pareto_fronts)
    return pareto_fronts

if __name__ == '__main__':
    penalty_c = 100

    floorplan_Bbox = gpd.read_file("Curated7/ICD_CalculatedV7.shp")
    floorplan_border = gpd.read_file("Curated7/OneLetterV7.shp")
    floorplan_skeleton = gpd.read_file("Curated7/SkelV7.shp")
    floorplan_corner = gpd.read_file("Curated7/CornersV7.shp")

    rows = [row for index, row in floorplan_Bbox.iterrows()]

    print(rows[0])

    num_cores = cpu_count()
    #start parallel processing
    with Pool(processes=num_cores) as pool:
        results = pool.map(process, rows)

    # remove None rows from the results
    results = [x for x in results if x is not None]
    print(results)



    # write down a single json containg letter, number of require vertices, number of selected vertices and the total weight
    # for each pareto fronts that has nonnegative total weight
    unique_pareto_fronts = []
    for pareto_front in results:
        if pareto_front['total_weight'] > 0:
            if pareto_front not in unique_pareto_fronts:
                unique_pareto_fronts.append(pareto_front)

    with open('pareto_fronts.json', 'w') as f:
        json.dump(unique_pareto_fronts, f)






