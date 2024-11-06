import json
import geopandas as gpd
from shapely.geometry import LineString
import momepy
import networkx as nx
import os
import random
from deap import base, creator, tools, algorithms
from computational_geometry_functions import edge_bearing, all_bearings_sorted, nearest_quarter, calculate_edge_angles, find_vertices_within_distance_from_points, check_if_the_node_is_uncertain, check_edge_uncertainty, uncertainty_weights_by_ID
from multiprocessing import Pool

def collect_results(result):
    pareto_fronts.append(result)
# Define the function to evaluate ILP
def evaluate_ILP(individual_edges):
    individual, edges = individual_edges
    Xc = [vertex for vertex, is_selected in zip(vertices, individual) if is_selected]
    Yodc = {}
    for edge in edges:
        for vertex in edge['all_vertices']:
            if vertex in Xc and vertex in edge['vertices']:
                Yodc[(edge['ID'], vertex)] = 1
            else:
                Yodc[(edge['ID'], vertex)] = 0
    total_weight = sum(edge['Weight'] * Yodc[(edge['ID'], vertex)] for edge in edges for vertex in edge['all_vertices'])
    num_selected_vertices = sum(individual)
    penalty = 0
    for edge in edges:
        if sum(Yodc[(edge['ID'], vertex)] for vertex in edge['vertices']) > 1:
            penalty = (sum(Yodc[(edge['ID'], vertex)] for vertex in edge['vertices'])) - 1
            total_weight += - penalty_c * penalty
    for edge in edges:
        for vertex in edge['vertices']:
            if Yodc[(edge['ID'], vertex)] == 1 and vertex not in Xc:
                penalty += - 1
                total_weight += - penalty_c * penalty
    return total_weight, num_selected_vertices, penalty

def process_floorplan(index, row):
    global same, different, all_edges_in_dateset, all_uncertain_edges_in_dateset, pareto_fronts
    global edgecount, uncertain_intersection
    polygon = row['geometry']
    letter = row['distinct']
    floorplan_intersecting = floorplan_skeleton[floorplan_skeleton.intersects(polygon)]
    floorplan_border_intersecting = floorplan_border[floorplan_border.intersects(polygon)]
    floorplan_Corners = floorplan_corner[floorplan_corner.intersects(polygon)]
    G = momepy.gdf_to_nx(floorplan_intersecting, approach="primal", length="length")
    H = G
    i = 0
    while i < 10:
        for u, v, key, attr in list(G.edges(keys=True, data=True)):
            if attr['length'] < 100:
                if G.has_edge(u, v):
                    if G.degree(u) >= G.degree(v):
                        G = nx.contracted_nodes(G, u, v, self_loops=False)
                    else:
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
                all_vertices_ID = list(floorplan_Corners['ID'])
                ILP.append({'ID': edge[2]['ID'], 'Weight': edge[2]['Weight'], 'vertices': vertices_ID, 'all_vertices': all_vertices_ID })
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
    toolbox.register("evaluate", evaluate_ILP, edges=edges)
    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutFlipBit, indpb=0.05)
    toolbox.register("select", tools.selNSGA2)
    NGEN = 5000
    MU = 100
    pop = toolbox.population(n=MU)
    try:
        algorithms.eaMuPlusLambda(pop, toolbox, mu=MU, lambda_=MU, cxpb=0.7, mutpb=0.2, ngen=NGEN, verbose=False)
    except:
        pass
    pareto_front = tools.sortNondominated(pop, len(pop), first_front_only=True)[0]
    for ind in pareto_front:
        total_weight, num_selected_vertices, penalties = evaluate_ILP((ind, edges))
        pareto_fronts.append({'letter': letter, 'total_weight': total_weight, 'num_selected_vertices': num_selected_vertices, 'penalties': penalties, 'uncertain_edges': len(edges), 'sum_of_weight': sum_of_weight})
    print(pareto_fronts)
    return pareto_fronts

if __name__ == "__main__":
    penalty_c = 50
    same = 0
    different = 0
    due = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    pareto_fronts = []
    edgecount = {}
    for i in range(101):
        edgecount[i] = 0
    all_edges_in_dateset = 0
    all_uncertain_edges_in_dateset = 0
    uncertain_intersection = {}
    for i in range(10):
        uncertain_intersection[i] = 0
    floorplan_Bbox = gpd.read_file("../data_synthetic/ICD_CalculatedV7.shp")
    floorplan_border = gpd.read_file("../data_synthetic/OneLetterV7.shp")
    floorplan_skeleton = gpd.read_file("../data_synthetic/SkelV7.shp")
    floorplan_corner = gpd.read_file("../data_synthetic/CornersV7.shp")
    folder = 'Geojson/'
    for filename in os.listdir(folder):
        file_path = os.path.join(folder, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                os.unlink(file_path)
        except Exception as e:
            print('Failed to delete %s. Reason: %s' % (file_path, e))
    counter = 0
    pool = Pool()  # Initialize multiprocessing Pool
    for index, row in floorplan_Bbox.iterrows():
        counter += 1
        if counter > 100:
            break
        pool.apply_async(process_floorplan, args=(index, row), callback=collect_results)

    pool.close()
    pool.join()
    unique_pareto_fronts = []
    for pareto_front in pareto_fronts:
        if pareto_front['total_weight'] > 0:
            if pareto_front not in unique_pareto_fronts:
                unique_pareto_fronts.append(pareto_front)
    with open('../pareto_fronts.json', 'w') as f:
        for element in unique_pareto_fronts:
            json.dump(element, f)
