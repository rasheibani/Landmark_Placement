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

from shapely.geometry import LineString, Point, Polygon
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


def nearest_quarter(ang, Grammar='4sectors'):
    # calculate the nearest quarter for an angle given a grammar
    if Grammar=='4sectors':
        ang = ang % 360
        return (ang // 90) * 90
    elif Grammar=='8sectors':
        ang = ang % 360
        return (ang // 45) * 45
    elif Grammar == '6sectors':
        ang = ang % 360

        if ang < 22.5:
            return 0
        elif ang < 67.5:
            return 22
        elif ang < 90 + 45:
            return 67.5
        elif ang < 225:
            return 135
        elif ang < 292.5:
            return 225
        else:
            return 0
    elif Grammar == 'Klippel':
        ang = ang % 360
        if ang < 5:
            return 0
        elif ang < 5+73:
            return 5
        elif ang < 5+73+20:
            return 5+73
        elif ang < 5+73+20+70:
            return 5+73+20
        elif ang < 5+73+20+70+24:
            return 5+73+20+70
        elif ang < 5+73+20+70+24+60:
            return 5+73+20+70+24
        elif ang < 5+73+20+70+24+60+39:
            return 5+73+20+70+24+60
        else:
            return 0

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


def check_if_the_node_is_uncertain(bearings_sorted, Grammar='4sectors'):
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
        nearest_quarters = [nearest_quarter(angle,Grammar=Grammar) for angle in angles_relative]
        # check if there is a quarter that is repeated 2 times or more
        if max(Counter(nearest_quarters).values()) >= 2 and len(bearings_sorted) > 1:
            return True
    return False

def check_if_two_outgoing_are_uncertain(out1, out2, Grammar='4sectors'):
    if nearest_quarter(out1 % 360,Grammar=Grammar) == nearest_quarter(out2 % 360,Grammar=Grammar):
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
            if check_if_two_outgoing_are_uncertain(bearings[i], bearings[i+1], Grammar='4sectors'):
                return True
        elif i == len(bearings)-1:
            if check_if_two_outgoing_are_uncertain(bearings[i], bearings[0], Grammar='4sectors'):
                return True
    return False

def check_edge_uncertainty(G, Node, Grammar='4sectors'):
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
                if check_if_two_outgoing_are_uncertain(outgoing_bearing1_corrected, outgoing_bearing2_corrected,Grammar=Grammar) and not np.allclose(outgoing_bearing1, outgoing_bearing2, atol=0.1):
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
