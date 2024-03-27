
from computational_geometry_functions import *


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
    if counter > 10000:
        break

    polygon = row['geometry']
    letter = row['distinct']
    # print(letter)

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
            if check_if_the_node_is_uncertain(bearings, Grammar=Grammar):
                uncertain_intersection[G.degree(Nodes)] += 1
                u_edge, u_coord = check_edge_uncertainty(G, Nodes, Grammar=Grammar)
                all_uncertain_edges_in_dateset += len(u_edge.keys())
                # print(str(G.degree(Nodes))+"from which"+str(len(check_edge_uncertainty(G, Nodes).keys())))
                # print(u_edge)
                # print(u_coord)
                All_uncertain_edges[Nodes] = u_coord
    # print('all:'+str (all_uncertain_edges_in_dateset))
    # print('uncertain:'+ str(all_edges_in_dateset) )

print("uncertain_intersection: " + str(uncertain_intersection))
print('all:'+str (all_uncertain_edges_in_dateset))
print('uncertain:'+ str(all_edges_in_dateset) )