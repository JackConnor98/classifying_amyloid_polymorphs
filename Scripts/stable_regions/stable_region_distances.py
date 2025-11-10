import pandas as pd
import numpy as np
import os
import sys
from plotnine import *
import matplotlib.colors as mcolors
from Bio.PDB import PDBParser
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Arc

################################
### Reading in sys arguments ###
################################

if len(sys.argv) > 1: 
    try: 
        distance_threshold = int(sys.argv[1]) 
        print(f"Distance Threshold: {distance_threshold}") 
    except ValueError: 
        print("Warning: invalid distance threshold provided — defaulting to 10.8.") 
        distance_threshold = 10.8 
else: 
    distance_threshold = 10.8 
    print("No distance threshold provided, defaulting to 10.8")

########################
### Creating folders ###
########################

# Creating output_dir
output_dir = os.path.join("Output", "stable_regions")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
region_dir = os.path.join(output_dir, "region_contacts")
if not os.path.exists(region_dir):
    os.makedirs(region_dir)
    
residue_dir = os.path.join(output_dir, "residue_contacts")
if not os.path.exists(residue_dir):
    os.makedirs(residue_dir)

#########################    
### Getting PDB files ###
#########################

# Finding all asym_units that need to be looped through
pdb_path = os.path.join("Output", "PDBs", "asymetric_unit")
pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]

# Getting the names of the PDB files in the asymmetric units folder
pdb_filenames = [os.path.splitext(os.path.basename(f))[0] for f in pdb_files]
pdb_names = [name.split("_")[0] for name in pdb_filenames]

######################
### Importing data ###
######################

# Loading High Resolution PDBs
high_resolution_PDBs = pd.read_csv(os.path.join("Output", "Validation", "high_resolution_pdb_ids.csv"))

# Loading RMSD Cluster Data
cluster_groups = pd.read_csv(os.path.join("Output", "RMSD", "data", "RMSD_cluster_groups.csv"))

# Loading Fibril Info
fibril_info = pd.read_csv(os.path.join("Output", "PDBs", "COM_and_fibril.csv"))
fibril_info = (
    fibril_info
    .drop(columns=["chain", "x_com", "y_com", "z_com", "Rg"])
    .drop_duplicates()
    )
fibril_df = pd.merge(fibril_info, cluster_groups, on = "pdb_id")

# Importing Stable Region Data
stable_regions = pd.read_csv(os.path.join("Output", "stable_regions", "stable_regions.csv"))
# Separating residues into start and end
stable_regions[["start", "end"]] = stable_regions["residues"].str.split("-", expand=True).astype(float)

########################################################################################
### Removing PDBs with a poor mean resolution - identified in the validation scripts ###
########################################################################################

# Creating PDB column
high_resolution_PDBs["pdb"] = high_resolution_PDBs["pdb_id"].str.split("_").str[0]

# Identifying the PDBs not found in high_resolution_PDBs
low_res = set(pdb_names) - set(high_resolution_PDBs["pdb"])

# Finding indices to remove
indices_to_remove = [i for i, name in enumerate(pdb_names) if name in low_res]

# Removing the corresponding elements from pdb_files
pdb_files_filtered = [f for i, f in enumerate(pdb_files) if i not in indices_to_remove]

##########################
### Defining Functions ###
##########################

def stable_region_distance(CA_filtered, stable_regions, selected_pdb):
    
    #######################################################
    ##### Calculating Distance Between Stable Regions #####
    #######################################################

    # Initialise empty region columns
    CA_filtered = CA_filtered.copy()
    CA_filtered["region1"] = np.nan
    CA_filtered["region2"] = np.nan

    # Assign region labels based on residue ranges
    for _, region_row in stable_regions.iterrows():
        start = region_row["start"]
        end = region_row["end"]
        region_name = region_row["region"]

        mask1 = (CA_filtered["resno1"] >= start) & (CA_filtered["resno1"] <= end)
        mask2 = (CA_filtered["resno2"] >= start) & (CA_filtered["resno2"] <= end)

        CA_filtered.loc[mask1, "region1"] = region_name
        CA_filtered.loc[mask2, "region2"] = region_name

    # Remove rows without stable region assignments
    stable_region_CA_dist = CA_filtered.dropna(subset=["region1", "region2"])

    # Remove comparisons between same region on same fibril
    stable_region_CA_dist = stable_region_CA_dist[
        ~((stable_region_CA_dist["region1"] == stable_region_CA_dist["region2"]) &
          (stable_region_CA_dist["which_fibril"] == "same"))
    ]

    # Find the minimum distance for each pair of regions per fibril
    min_distances = (
        stable_region_CA_dist
        .groupby(["region1", "region2", "which_fibril", "pdb_id"], as_index=False)
        .agg(min_distance=("distance", "min"))
    )

    # Add PDB name column
    min_distances["pdb"] = selected_pdb

    return min_distances

################################
### Looping through all PDBs ###
################################

parser = PDBParser(QUIET=True)

# Initialising dataframe
combined_residue_distance = pd.DataFrame()
combined_region_distance = pd.DataFrame()

for file in pdb_files_filtered:
    # Extracting PDB code from file name
    tmp = os.path.basename(file)
    tmp = os.path.splitext(tmp)[0]
    selected_pdb = tmp.split("_")[0]

    print(f"\nAnalysing: {selected_pdb}\n")

    # --- Read PDB file ---
    structure = parser.get_structure(selected_pdb, os.path.join(pdb_path, file))

    # --- Extract atom info ---
    atom_data = []
    for atom in structure.get_atoms():
        parent = atom.get_parent()
        if atom.get_name() == "CA":  # Only CA atoms
            atom_data.append({
                "PDB": selected_pdb,
                "fibril": atom.bfactor,   # fibril number is stored as the B-factor
                "resno": parent.get_id()[1],
                "x": atom.coord[0],
                "y": atom.coord[1],
                "z": atom.coord[2]
            })

    df = pd.DataFrame(atom_data)

    # Merge with fibril_df
    df = df.merge(fibril_df, on=["PDB", "fibril"])
    polymorph_count = df["polymorph"].max()

    # --- Calculate pairwise distances between CA atoms ---
    coords = df[["x", "y", "z"]].to_numpy()
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    CA_distances = np.linalg.norm(diff, axis=2)

    # Convert to long format
    n = len(df)
    CA_distance_df = pd.DataFrame({
        "resno1": np.repeat(df["resno"].values, n),
        "fibril1": np.repeat(df["fibril"].values, n),
        "polymorph1": np.repeat(df["polymorph"].values, n),
        "resno2": np.tile(df["resno"].values, n),
        "fibril2": np.tile(df["fibril"].values, n),
        "polymorph2": np.tile(df["polymorph"].values, n),
        "distance": CA_distances.flatten()
    })

    # Remove self-comparisons and near neighbours on same fibril
    CA_distance_df = CA_distance_df[
        ~((CA_distance_df["resno1"] == CA_distance_df["resno2"]) &
          (CA_distance_df["fibril1"] == CA_distance_df["fibril2"]))
    ]
    CA_distance_df = CA_distance_df[
        ~((abs(CA_distance_df["resno1"] - CA_distance_df["resno2"]) <= 3) &
          (CA_distance_df["fibril1"] == CA_distance_df["fibril2"]))
    ]

    # Add which_fibril column
    CA_distance_df["which_fibril"] = np.where(
        CA_distance_df["fibril1"] == CA_distance_df["fibril2"],
        "same", "different"
    )

    # --- If only one polymorph ---
    if polymorph_count == 1:
        min_distances = (
            CA_distance_df.groupby(["resno1", "resno2"], as_index=False)
            .agg(min_distance=("distance", "min"))
        )

        CA_filtered = (
            CA_distance_df.merge(min_distances, on=["resno1", "resno2"])
            .query("distance == min_distance")
            .drop(columns=["min_distance", "fibril1", "fibril2"])
        )

        if len(CA_filtered) > 0:
            current_pdb_info = fibril_df.query("PDB == @selected_pdb")
            CA_filtered["pdb_id"] = current_pdb_info["pdb_id"].iloc[0]

        combined_residue_distance = pd.concat([combined_residue_distance, CA_filtered], ignore_index=True)

        # Calculate stable region distances
        min_region_distances = stable_region_distance(CA_filtered, stable_regions, selected_pdb)
        combined_region_distance = pd.concat([combined_region_distance, min_region_distances], ignore_index=True)

    else:
        # Multiple polymorphs
        remaining_polymorphs = CA_distance_df["polymorph1"].unique()

        for i in remaining_polymorphs:
            CA_distance_polymorph = CA_distance_df.query("polymorph1 == @i")

            min_distances = (
                CA_distance_polymorph.groupby(["resno1", "resno2"], as_index=False)
                .agg(min_distance=("distance", "min"))
            )

            CA_filtered = (
                CA_distance_polymorph.merge(min_distances, on=["resno1", "resno2"])
                .query("distance == min_distance")
                .drop(columns=["min_distance", "fibril1", "fibril2"])
            )

            if len(CA_filtered) > 0:
                current_pdb_info = fibril_df.query("PDB == @selected_pdb") \
                                            .drop_duplicates(subset=["pdb_id", "polymorph"]) \
                                            .drop(columns=["fibril"])
                CA_filtered["pdb_id"] = current_pdb_info.loc[
                    current_pdb_info["polymorph"] == i, "pdb_id"
                ].values[0]

            combined_residue_distance = pd.concat([combined_residue_distance, CA_filtered], ignore_index=True)
            min_region_distances = stable_region_distance(CA_filtered, stable_regions, selected_pdb)
            combined_region_distance = pd.concat([combined_region_distance, min_region_distances], ignore_index=True)

# Saving combined_residue_distances
combined_residue_distance.to_csv(os.path.join(residue_dir, "residue_distances.csv"), index = False)

# Saving combined_region_distances
combined_region_distance.to_csv(os.path.join(region_dir, "region_distances.csv"), index = False)


#######################################################################
##### Counting Stable Region Contacts for Each RMSD Cluster Group #####
#######################################################################

# Creating a single comparison column
combined_region_distance["comparison"] = (
    combined_region_distance["region1"].fillna(0).astype(int).astype(str) + "--" +
    combined_region_distance["region2"].fillna(0).astype(int).astype(str)
)
# Merging stable region distances and cluster group dataframe
stable_dist_by_cluster = combined_region_distance.merge(cluster_groups, on="pdb_id")

# Saving data
stable_dist_by_cluster.to_csv(os.path.join(region_dir, "region_distances_by_cluster_group.csv"), index=False)


# Keep contacts with min_distance < distance_threshold 
contact_df = stable_dist_by_cluster[stable_dist_by_cluster["min_distance"] < distance_threshold]

# Count number of each contact per group
contact_counts = (
    contact_df.groupby(["group", "comparison"])
    .size()
    .reset_index(name="n")
    .sort_values("group")
)

# Count number of unique PDBs in each group
pdbs_per_group = (
    stable_dist_by_cluster.groupby("group")["pdb_id"]
    .nunique()
    .reset_index(name="unique_pdb_count")
)

# Merge counts with PDB totals
percent_df = contact_counts.merge(pdbs_per_group, on="group")

# Calculate frequency percentage
percent_df["freq"] = (percent_df["n"] / percent_df["unique_pdb_count"]) * 100

# Saving to .csv
percent_df.to_csv(os.path.join(region_dir, "region_contacts_percentage_per_group.csv"), index=False)

#####################
### Network Plots ###
#####################

### Plotting contacts between stable regions

stable_dist_by_cluster["group"] = stable_dist_by_cluster["group"].astype(int)
stable_dist_by_cluster["min_distance"] = stable_dist_by_cluster["min_distance"].astype(int)
stable_dist_by_cluster["region1"] = stable_dist_by_cluster["region1"].astype(int)
stable_dist_by_cluster["region2"] = stable_dist_by_cluster["region2"].astype(int)

# Load the qualitative palettes
set1 = plt.get_cmap("Set1")
dark2 = plt.get_cmap("Dark2")

# Convert to hex
set1_colours = [mcolors.to_hex(set1(i)) for i in range(set1.N)]
dark2_colours = [mcolors.to_hex(dark2(i)) for i in range(dark2.N)]

# Combine them
base_colours = set1_colours + dark2_colours

# Expand to match number of stable regions
region_colours = np.tile(base_colours, int(np.ceil(len(stable_regions) / len(base_colours))))[:len(stable_regions)]

for group in np.unique(stable_dist_by_cluster["group"]):

    # Filter network data by distance threshold
    network_data = stable_dist_by_cluster[
        (stable_dist_by_cluster["min_distance"] <= distance_threshold) &
        (stable_dist_by_cluster["group"] == group)
    ]

    if not network_data.empty:
        # --- Extract unique vertex names from your data ---
        vertex_names = sorted(np.unique(stable_dist_by_cluster["region1"]))

        # --- Create complete graph with string nodes ---
        G_complete = nx.complete_graph(len(vertex_names))
        mapping = dict(zip(range(len(vertex_names)), vertex_names))
        G_complete = nx.relabel_nodes(G_complete, mapping)

        # --- Remove all edges initially ---
        G_complete.remove_edges_from(list(G_complete.edges))

        # --- Node colours ---
        node_colors = {v: region_colours[i % len(region_colours)] for i, v in enumerate(vertex_names)}

        # --- Custom circular layout (1 at top, increasing clockwise) ---
        n_nodes = len(vertex_names)
        angle_step = 2 * np.pi / n_nodes
        angles = [np.pi / 2 - i * angle_step for i in range(n_nodes)]  # start from top
        pos = {node: (np.cos(angle), np.sin(angle)) for node, angle in zip(vertex_names, angles)}

          # --- Add weighted edges ---
        for _, row in network_data.iterrows():
            u, v = row["region1"], row["region2"]
            if G_complete.has_edge(u, v):
                G_complete[u][v]["weight"] += 1
            else:
                G_complete.add_edge(u, v, weight = 1, which_fibril = row["which_fibril"])
                
        # --- Draw network ---
        plt.figure(figsize = (5, 5))

        # Draw normal edges (excluding self-loops)
        edges_to_draw = [(u, v, d) for u, v, d in G_complete.edges(data = True) if u != v]
        nx.draw_networkx_edges(
            G_complete, pos,
            edgelist = [(u, v) for u, v, _ in edges_to_draw],
            width = [1 + d["weight"] * 0.3 for _, _, d in edges_to_draw],
            edge_color = ["grey" if d["which_fibril"] == "same" else "red" for _, _, d in edges_to_draw],
            alpha = 0.6
        )

        # Draw simple outward-facing self-loops
        ax = plt.gca()
        for node, (x, y) in pos.items():
            if G_complete.has_edge(node, node):
                data = G_complete[node][node]
                weight = data.get("weight", 1)
                fibril_type = data.get("which_fibril", "same")

                # Get outward angle (from origin)
                angle_deg = np.degrees(np.arctan2(y, x))

                # Slightly push the loop away from the node
                offset = 0.05
                loop_radius = 0.1
                loop_center_x = x + offset * np.cos(np.radians(angle_deg))
                loop_center_y = y + offset * np.sin(np.radians(angle_deg))

                # Adjust arc so it connects visually to node edge
                loop = Arc(
                    (loop_center_x, loop_center_y),
                    width = loop_radius, height = loop_radius * 5,
                    angle = angle_deg - 270,
                    theta1 = 200, theta2 = 340,
                    color = "grey" if fibril_type == "same" else "red",
                    lw = 1 + weight * 0.3,  # <--- scale linewidth like normal edges
                    alpha = 0.9
                )
                ax.add_patch(loop)

        # Draw nodes and labels
        nx.draw_networkx_nodes(G_complete, pos, node_color = [node_colors[n] for n in G_complete.nodes()], node_size = 500)
        nx.draw_networkx_labels(G_complete, pos, font_size = 12, font_color = "black")

        plt.title(f"Group {group}", fontsize = 20)
        plt.axis("off")
        plt.xlim(-1.6, 1.6)
        plt.ylim(-1.6, 1.6)

        # --- Save network plot ---
        plt.savefig(os.path.join(region_dir, f"group_{group}_network.png"), dpi = 300)
        plt.close()
        
        
        
        
        
        
        
        
        
#################################################################
##### Counting Residue Contacts for Each RMSD Cluster Group #####
#################################################################

# Merging stable residue distances and cluster group dataframe
residue_dist_by_cluster = combined_residue_distance.merge(cluster_groups, on="pdb_id")

# Saving data
residue_dist_by_cluster.to_csv(os.path.join(residue_dir, "residue_distances_by_cluster_group.csv"), index=False)

#####################
### Network Plots ###
#####################

### Plotting contacts between residues

residue_dist_by_cluster["group"] = residue_dist_by_cluster["group"].astype(int)
residue_dist_by_cluster["distance"] = residue_dist_by_cluster["distance"].astype(int)
residue_dist_by_cluster["resno1"] = residue_dist_by_cluster["resno1"].astype(int)
residue_dist_by_cluster["resno2"] = residue_dist_by_cluster["resno2"].astype(int)


for group in np.unique(residue_dist_by_cluster["group"]):

    # Filter network data by distance threshold
    network_data = residue_dist_by_cluster[
        (residue_dist_by_cluster["distance"] <= distance_threshold) &
        (residue_dist_by_cluster["group"] == group)
    ]

    if not network_data.empty:
        # --- Extract unique vertex names from your data ---
        vertex_names = sorted(np.unique(residue_dist_by_cluster["resno1"]))

        # --- Create complete graph with string nodes ---
        G_complete = nx.complete_graph(len(vertex_names))
        mapping = dict(zip(range(len(vertex_names)), vertex_names))
        G_complete = nx.relabel_nodes(G_complete, mapping)

        # --- Remove all edges initially ---
        G_complete.remove_edges_from(list(G_complete.edges))

        # --- Node colours ---
        node_colors = {v: region_colours[i % len(region_colours)] for i, v in enumerate(vertex_names)}

        # --- Custom circular layout (1 at top, increasing clockwise) ---
        n_nodes = len(vertex_names)
        angle_step = 2 * np.pi / n_nodes
        angles = [np.pi / 2 - i * angle_step for i in range(n_nodes)]  # start from top
        pos = {node: (np.cos(angle), np.sin(angle)) for node, angle in zip(vertex_names, angles)}

          # --- Add weighted edges ---
        for _, row in network_data.iterrows():
            u, v = row["resno1"], row["resno2"]
            if G_complete.has_edge(u, v):
                G_complete[u][v]["weight"] += 1
            else:
                G_complete.add_edge(u, v, weight = 1, which_fibril = row["which_fibril"])
                
        # --- Draw network ---
        plt.figure(figsize = (10, 10))

        # Draw normal edges (excluding self-loops)
        edges_to_draw = [(u, v, d) for u, v, d in G_complete.edges(data = True) if u != v]
        nx.draw_networkx_edges(
            G_complete, pos,
            edgelist = [(u, v) for u, v, _ in edges_to_draw],
            width = [1 + d["weight"] * 0.15 for _, _, d in edges_to_draw],
            edge_color = ["grey" if d["which_fibril"] == "same" else "red" for _, _, d in edges_to_draw],
            alpha = 0.6
        )

        # Draw simple outward-facing self-loops
        ax = plt.gca()
        for node, (x, y) in pos.items():
            if G_complete.has_edge(node, node):
                data = G_complete[node][node]
                weight = data.get("weight", 1)
                fibril_type = data.get("which_fibril", "same")

                # Get outward angle (from origin)
                angle_deg = np.degrees(np.arctan2(y, x))

                # Slightly push the loop away from the node
                offset = 0.01
                loop_radius = 0.03
                loop_center_x = x + offset * np.cos(np.radians(angle_deg))
                loop_center_y = y + offset * np.sin(np.radians(angle_deg))

                # Adjust arc so it connects visually to node edge
                loop = Arc(
                    (loop_center_x, loop_center_y),
                    width = loop_radius, height = loop_radius * 8,
                    angle = angle_deg - 270,
                    theta1 = 200, theta2 = 340,
                    color = "grey" if fibril_type == "same" else "red",
                    lw = 1 + weight * 0.15,  # <--- scale linewidth like normal edges
                    alpha = 0.9
                )
                ax.add_patch(loop)

        # Draw nodes and labels
        nx.draw_networkx_nodes(G_complete, pos, node_color = [node_colors[n] for n in G_complete.nodes()], node_size = 250)
        nx.draw_networkx_labels(G_complete, pos, font_size = 9, font_color = "black")

        plt.title(f"Group {group}", fontsize = 20)
        plt.axis("off")
        plt.xlim(-1.2, 1.2)
        plt.ylim(-1.2, 1.2)

        # --- Save network plot ---
        plt.savefig(os.path.join(residue_dir, f"group_{group}_network.png"), dpi = 300)
        plt.close()