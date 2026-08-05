#%%
#  This script runs on motile and the motile toolbox
# https://github.com/funkelab/motile
# https://github.com/funkelab/motile_toolbox

# Motile z-stack tracking with 2D and interactive 3D visualization

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
from mpl_toolkits.mplot3d import Axes3D

logger = logging.getLogger(__name__)

def nodes_from_points_list(points_list: np.ndarray):
    cand_graph = nx.DiGraph()
    node_frame_dict = {}
    print("Extracting nodes from points list")

    for i, point in enumerate(points_list):
        t = point[0]
        pos = point[1:]
        cand_graph.add_node(i,
                            **{NodeAttr.TIME.value: t,
                               NodeAttr.POS.value: pos})
        node_frame_dict.setdefault(t, []).append(i)
    return cand_graph, node_frame_dict


def get_candidate_graph_from_points_list(points_list, max_edge_distance):
    cand_graph, node_frame_dict = nodes_from_points_list(points_list)

    add_cand_edges(
        cand_graph,
        max_edge_distance=max_edge_distance,
        node_frame_dict=node_frame_dict,
    )
    return cand_graph


if __name__ == "__main__":

    acinus_path = r"...add path.../Liv1/Lobule1/acinus0"

    stack_paths = glob.glob(f"{acinus_path}/*/average_properties_per_cell.csv")

    print(f"Found {len(stack_paths)} CSV files")

    vol_positions = pd.DataFrame()
    acinus_props = pd.DataFrame()

    for i, stack in enumerate(stack_paths):
        cell_props = pd.read_csv(stack)

        stack_positions = cell_props[["centroid-0", "centroid-1"]]
        stack_positions = stack_positions.assign(stack=i)

        vol_positions = pd.concat((vol_positions, stack_positions), ignore_index=True)
        acinus_props = pd.concat((acinus_props, cell_props), ignore_index=True)

    node_list = np.asarray(vol_positions)
    node_list[:] = node_list[:, [2, 0, 1]]

    graph = get_candidate_graph_from_points_list(node_list, 200)

    solver = motile.Solver(
        TrackGraph(graph, frame_attribute=NodeAttr.TIME.value)
    )

    solver.add_costs(
        EdgeDistance(weight=1.0,
                     constant=-200,
                     position_attribute="pos")
    )

    solver.add_costs(Appear(constant=1.0))
    solver.add_constraints(MaxParents(1))
    solver.add_constraints(MaxChildren(1))

    solution = solver.solve()

    solution_graph = solver.get_selected_subgraph(solution=solution)
    output_graph = graph_to_nx(solution_graph)

    node_connections = np.zeros((len(node_list), 1), dtype=int)

    for i, tracklet in enumerate(nx.weakly_connected_components(output_graph)):
        node_connections[list(tracklet)] = i

    connected_nodes = np.hstack((node_list, node_connections))

    ntracks = int(node_connections.max()) + 1
    cmap = plt.get_cmap("gist_ncar", max(ntracks, 2))

    # ===========================
    # 2D plot
    # ===========================
    fig, ax = plt.subplots(figsize=(9, 9))

    ax.scatter(
        connected_nodes[:, 2],
        connected_nodes[:, 1],
        c=connected_nodes[:, 3],
        cmap=cmap,
        s=20,
        alpha=0.8,
    )

    for u, v in output_graph.edges():
        color = cmap(int(connected_nodes[u, 3]))
        ax.plot(
            [connected_nodes[u, 2], connected_nodes[v, 2]],
            [connected_nodes[u, 1], connected_nodes[v, 1]],
            color=color,
            linewidth=1,
            alpha=0.9,
        )

    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title("Linked Cells with Selected Edges (2D)")
    ax.set_aspect("equal")

    # ===========================
    # Interactive 3D plot
    # ===========================
    fig3 = plt.figure(figsize=(10, 9))
    ax3 = fig3.add_subplot(111, projection="3d")

    ax3.scatter(
        connected_nodes[:, 2],
        connected_nodes[:, 1],
        connected_nodes[:, 0],
        c=connected_nodes[:, 3],
        cmap=cmap,
        s=18,
        depthshade=True,
    )

    for u, v in output_graph.edges():
        color = cmap(int(connected_nodes[u, 3]))
        ax3.plot(
            [connected_nodes[u, 2], connected_nodes[v, 2]],
            [connected_nodes[u, 1], connected_nodes[v, 1]],
            [connected_nodes[u, 0], connected_nodes[v, 0]],
            color=color,
            linewidth=1,
            alpha=0.9,
        )

    # Set limits from the data
    ax3.set_xlim(connected_nodes[:,2].min(), connected_nodes[:,2].max())
    ax3.set_ylim(connected_nodes[:,1].min(), connected_nodes[:,1].max())
    ax3.set_zlim(connected_nodes[:,0].min(), connected_nodes[:,0].max())

    # Make the 3D box proportional to the data

    ax3.set_box_aspect((5, 2, 3))
   

    ax3.set_xlabel("X")
    ax3.set_ylabel("Y")
    ax3.set_zlabel("Stack")
    ax3.set_title("Linked Cells with Selected Edges (3D)")
    ax3.view_init(elev=25, azim=-60)

    plt.show()

   
# ==========================================================
# Save the original table with linked cell IDs
# ==========================================================

acinus_props = acinus_props.assign(cell_id_linked=node_connections)
acinus_props = acinus_props.assign(stack_id=connected_nodes[:,0])

acinus_props.to_csv(
    f"{acinus_path}/cell_props_for_cell_id_linked.csv",
    index=False,
)

# ==========================================================
# Generate one averaged row per linked cell
# ==========================================================

print("Generating averaged cell table...")

# Average every numeric column
mean_cells = (
    acinus_props
    .groupby("cell_id_linked")
    .mean(numeric_only=True)
)

# Remove numeric columns that should not be averaged
mean_cells = mean_cells.drop(
    columns=["stack_id", "label", "Unnamed: 0"],
    errors="ignore"
)

# Keep the first value of all non-numeric columns
metadata = (
    acinus_props
    .groupby("cell_id_linked")
    .first()
    .drop(columns=mean_cells.columns, errors="ignore")
)

# Remove columns that should not appear in the averaged table
metadata = metadata.drop(
    columns=["stack_id", "label", "Unnamed: 0"],
    errors="ignore"
)

# Keep all stack IDs contributing to each reconstructed cell
stack_ids = (
    acinus_props
    .groupby("cell_id_linked")["stack_id"]
    .apply(list)
)

# Combine everything
mean_cells = pd.concat(
    [
        metadata,
        mean_cells,
        stack_ids.rename("stack_id"),
    ],
    axis=1,
).reset_index()

# Save averaged table
mean_cells.to_csv(
    f"{acinus_path}/mean_data_for_linked_cell_masks.csv",
    index=False,
)

print(f"Generated {len(mean_cells)} averaged cells.")
