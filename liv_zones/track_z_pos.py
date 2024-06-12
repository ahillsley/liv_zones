# This script runs on motile
# https://github.com/funkelab/motile
import numpy as np
import glob
import pandas as pd
import networkx as nx
from typing import Any
import logging
from motile_toolbox.candidate_graph.utils import add_cand_edges
from motile_toolbox.candidate_graph.graph_attributes import NodeAttr
from motile_toolbox.candidate_graph import graph_to_nx
import motile
from motile.constraints import MaxParents, MaxChildren
from motile.costs import EdgeDistance, Appear
from motile import TrackGraph
import matplotlib.pyplot as plt

logger = logging.getLogger(__name__)

# Define functions to build a candidate graph from a list of points.
# in this case the points are the cell centers, and the time dimension represents 
# the stacks (z)
# the following are functions taken from motile_toolbox
# https://github.com/funkelab/motile_toolbox

def nodes_from_points_list(
    points_list: np.ndarray,
) -> tuple[nx.DiGraph, dict[int, list[Any]]]:
    """Extract candidate nodes from a list of points. Uses the index of the
    point in the list as its unique id.
    Returns a networkx graph with only nodes, and also a dictionary from frames to
    node_ids for efficient edge adding.

    Args:
        points_list (np.ndarray): An NxD numpy array with N points and D
            (3 or 4) dimensions. Dimensions should be in order (t, [z], y, x).

    Returns:
        tuple[nx.DiGraph, dict[int, list[Any]]]: A candidate graph with only nodes,
            and a mapping from time frames to node ids.
    """
    cand_graph = nx.DiGraph()
    # also construct a dictionary from time frame to node_id for efficiency
    node_frame_dict: dict[int, list[Any]] = {}
    print("Extracting nodes from points list")
    for i, point in enumerate(points_list):
        # assume t, [z], y, x
        t = point[0]
        pos = point[1:]
        node_id = i
        attrs = {
            NodeAttr.TIME.value: t,
            NodeAttr.POS.value: pos,
        }
        cand_graph.add_node(node_id, **attrs)
        if t not in node_frame_dict:
            node_frame_dict[t] = []
        node_frame_dict[t].append(node_id)
    return cand_graph, node_frame_dict

def get_candidate_graph_from_points_list(
    points_list: np.ndarray,
    max_edge_distance: float,
) -> nx.DiGraph:
    """Construct a candidate graph from a points list.

    Args:
        points_list (np.ndarray): An NxD numpy array with N points and D
            (3 or 4) dimensions. Dimensions should be in order  (t, [z], y, x).
        max_edge_distance (float): Maximum distance that objects can travel between
            frames. All nodes with centroids within this distance in adjacent frames
            will by connected with a candidate edge.

    Returns:
        tuple[nx.DiGraph, list[set[Any]] | None]: A candidate graph that can be passed
        to the motile solver, and a list of conflicting node ids.
    """
    # add nodes
    cand_graph, node_frame_dict = nodes_from_points_list(points_list)
    logger.info(f"Candidate nodes: {cand_graph.number_of_nodes()}")
    # add edges
    add_cand_edges(
        cand_graph,
        max_edge_distance=max_edge_distance,
        node_frame_dict=node_frame_dict,
    )
    return cand_graph


if __name__ == '__main__':
    # for a given folder, read in all the avg_properties_per_cell_positions
    acinus_path = '/groups/feliciano/felicianolab/For_Alex_and_Mark/FOR_ALEX/Male/Liv1/Lobule1/acinus0'
    stack_paths = glob.glob(f'{acinus_path}/*/average_properties_per_cell.csv')

    # Extract only the x (centroid-1) and y (centroid-0) positions from the 
    # individual stack csvs

    vol_positions = pd.DataFrame()
    acinus_props = pd.DataFrame()
    for i, stack in enumerate(stack_paths):
        cell_props = pd.read_csv(stack)
        stack_positions = cell_props[['centroid-0', 'centroid-1']]
        stack_positions = stack_positions.assign(stack=i)
        vol_positions = pd.concat((vol_positions, stack_positions))
        acinus_props = pd.concat((acinus_props, cell_props))

    node_list = np.asarray(vol_positions)
    node_list[:] = node_list[:,[2,0,1]]

    # create a candidate graph from the list of cell centers
    # where a potential edge is drawn
    graph = get_candidate_graph_from_points_list(node_list, 200)
    solver = motile.Solver(
        TrackGraph(graph, frame_attribute=NodeAttr.TIME.value)
        )

    # add a cost to the distance i.e. weight it so that
    # closer points are more likely to be part of the same track
    solver.add_costs(
        EdgeDistance(
            weight=1.0,
            constant=-200,
            position_attribute='pos'))

    # add a small penalty to start a new track
    solver.add_costs(Appear(constant=1.0))

    # add constraints on the solution (no splits, no merges)
    solver.add_constraints(MaxParents(1))
    solver.add_constraints(MaxChildren(1))

    # solve and convert output into a nx graph 
    solution = solver.solve()
    solution_graph = solver.get_selected_subgraph(solution=solution)
    output_graph = graph_to_nx(solution_graph)

    # find the tracks
    node_connections = np.zeros((len(node_list),1))
    for i, tracklet in enumerate(nx.weakly_connected_components(output_graph)):
        # label nodes by the track they are a part of
        node_connections[list(tracklet)] = i

    connected_nodes = np.hstack((node_list, node_connections))

    # visualize the results
    plt.scatter(connected_nodes[:,2], connected_nodes[:,1], c=connected_nodes[:,3], alpha=0.5)

    # save_results as a new csv
    # Note that cell_id_linked = 0 means that the cell is not part of a track
    # i.e. that cell was only found in a single slice
    acinus_props = acinus_props.assign(cell_id_linked=node_connections)
    acinus_props.to_csv(f'{acinus_path}/avg_cell_props_cell_ids_linked.csv')
