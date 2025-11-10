import os
import sys
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, fcluster
import matplotlib.pyplot as plt
from plotnine import *

# Setting up paths
output_dir = os.path.join("Output", "RMSD")
data_path = os.path.join(output_dir, "data")
single_pdb_dir = os.path.join(output_dir, "single_reference")

# Create output directories if they don't exist
if not os.path.exists(single_pdb_dir):
    os.mkdir(single_pdb_dir)

### Reading in command line arguments ###

if len(sys.argv) > 1: 
    try: 
        custom_cut_height = int(sys.argv[1]) 
        print(f"Cut Height: {custom_cut_height}") 
    except ValueError: 
        print("Warning: invalid cut height provided — defaulting to mean Euclidean distance") 
        custom_cut_height = 0 
else: 
    custom_cut_height = 0 
    print("No cut height provided, defaulting to mean Euclidean distance")


# Importing data
df = pd.read_csv(os.path.join(data_path, "pairwise_rmsd.csv"), sep=",")
high_resolution_residues = pd.read_csv(os.path.join("Output", "Validation", "high_resolution_residues.csv"), sep=",")


### Removing low resolution PDBs from df ###

amyloid_names = high_resolution_residues['pdb_id'].astype(str).unique()

filtered_df = df[
    df['ref_name'].isin(amyloid_names) &
    df['mob_name'].isin(amyloid_names)
].copy()

##############################################################################

### Plotting the RMSD for each reference structure ###

for i in amyloid_names:
    
    # Selecting data for the current reference structure
    x = filtered_df[filtered_df['ref_name'] == i]

    # Sorting RMSD in ascending order
    x = x.sort_values(by='rmsd', ascending=True)

    # Making pdb_id an ordered categorical variable for plotting
    x["mob_name"] = pd.Categorical(
        x["mob_name"],
        categories=x["mob_name"].tolist(),
        ordered=True
    )

    p = (
        ggplot(x, aes(x='mob_name', y='rmsd')) 
        + geom_point(size = 4)
        + labs(title=f'RMSD values for reference structure {i}', 
               x='PDB ID', 
               y='RMSD (Å)')
        + scale_y_continuous(limits=(0, np.max(x['rmsd']) * 1.05), expand = (0, 0, 0.1, 0), 
                             breaks=range(0, 999999, 5)) 
        + theme_bw()
        + theme(
            panel_grid_major_y = element_line(color = "grey", size = 0.25, linetype = "dashed"),
            panel_grid_major_x = element_blank(),
            panel_border = element_rect(color = "black", size = 1.5),
            plot_title = element_text(size = 20, face = "bold", ha = "center"),
            axis_title = element_text(size = 18, face = "bold"),
            axis_text_y = element_text(size = 14, colour = "black"),
            axis_text_x = element_text(size = 13, angle = 90, vjust = 0.5, hjust=1, colour = "black")
        )
    )

    # Saving plot
    p.save(os.path.join(single_pdb_dir, f'{i}_RMSD.png'), width=20, height=8, dpi=300)

##############################################################################

# Pivot the DataFrame to wide format
df_wide = filtered_df.pivot(index='ref_name', columns='mob_name', values='rmsd').reset_index()

# Rename 'ref_name' to 'pdb_id'
df_wide = df_wide.rename(columns={'ref_name': 'pdb_id'})

# Set 'pdb_id' as the index
cluster_data = df_wide.set_index('pdb_id')

### Handling PDBs with no overlapping residues by setting the RMSD to mean + 3*SD ###

# Compute mean + 3*SD for all numeric entries, ignoring NaN
NA_distance = np.nanmean(cluster_data.values) + 3 * np.nanstd(cluster_data.values)
# Replace NaNs with NA_distance
cluster_data = cluster_data.fillna(NA_distance)

# Scale the data (mean = 0, std = 1)
scaled_cluster_data = (cluster_data - cluster_data.mean()) / cluster_data.std()

# Save the scaled data to a CSV file
scaled_cluster_data.to_csv(os.path.join(data_path, "scaled_cluster_data.csv"))

############################################################################################

################################
### Clustering Based On RMSD ###
################################

# Calculate Euclidean distance matrix
cluster_dist = pdist(scaled_cluster_data, metric='euclidean')

# Perform hierarchical clustering using average linkage
hc_average = linkage(cluster_dist, method='average')


### Comparing number of clusters vs cut height ###

# Extract linkage heights
heights = hc_average[:, 2]  # third column of linkage matrix

# Create scree_data equivalent
scree_data = pd.DataFrame({
    "height": heights,
    "groups": np.arange(len(heights), 0, -1)
})

# Save to CSV
scree_data.to_csv(os.path.join(data_path, "scree_data.csv"), index=False)

# Plotting the scree plot
p = (
    ggplot(scree_data, aes(x='groups', y='height'))
    + geom_point(size=1.5) 
    + geom_line(size=0.5)
    + annotate("segment", x=1, xend=len(heights), y=float(custom_cut_height), yend=float(custom_cut_height), color="red", linetype="dashed", size=0.75)
    + scale_x_continuous(breaks=range(0, len(heights)+1, 10), expand = (0, 0.05))
    + scale_y_continuous(expand = (0, 0, 0.05, 0))
    + labs(title="", 
           x="Groups", 
           y="Euclidean Distance (Å)")
    + theme_bw()
    + theme(
        panel_border = element_rect(color = "black", size = 1.5),
        panel_grid_major = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
        panel_grid_minor = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
        axis_title = element_text(size = 16, face = "bold"),
        axis_text = element_text(size = 14, colour = "black")
    )
)

# Saving plot
p.save(os.path.join(output_dir, "scree_plot.png"), height=4, width=8, dpi=300)

##################
### Dendrogram ###
##################

# Get the dendrogram order
leaf_order = leaves_list(hc_average)

# Getting a list of PDB IDs in the order of the dendrogram leaves
pdb_id_labels = df_wide["pdb_id"].iloc[leaf_order].reset_index(drop=True)

# Create the dendrogram order DataFrame
dendrogram_order = pd.DataFrame({
    "pdb_id": pdb_id_labels,
    "order": np.arange(1, len(pdb_id_labels) + 1)
})

#### Cutting the dendrogram ###

# Setting the cut height
if float(custom_cut_height) == 0:
    cut_height = np.mean(cluster_dist)
    print(f"Using mean Euclidean distance as the cut height: {cut_height}")
else:
    cut_height = float(custom_cut_height)
    print(f"Using custom cut height: {cut_height}")

# Assign clusters based on the cut height
cluster_assignments = fcluster(hc_average, t=cut_height, criterion='distance')

# Create a DataFrame for cluster assignments
cluster_assignments = pd.DataFrame({
    "pdb_id": df_wide["pdb_id"],
    "cluster": cluster_assignments
})

# Merge dendrogram order with cluster assignments
cluster_groups = pd.merge(dendrogram_order, cluster_assignments, on="pdb_id")

# Identify the order in which clusters appear from left to right
unique_cluster_order = cluster_groups["cluster"].drop_duplicates().tolist()

# Create a mapping: first seen cluster -> 1, next -> 2, etc.
cluster_mapping = {old: new for new, old in enumerate(unique_cluster_order, start=1)}

# Apply the renumbering
cluster_groups["group"] = cluster_groups["cluster"].map(cluster_mapping)

# Remove the old cluster column
cluster_groups = cluster_groups.drop(columns=["cluster"])

# Save the final cluster groups to a CSV file
cluster_groups.to_csv(os.path.join(data_path, "RMSD_cluster_groups.csv"), index=False)

### Plotting the dendrogram ###

# plot dendrogram
plt.figure(figsize=(12,6))
d = dendrogram(hc_average, 
               labels=df_wide["pdb_id"].tolist(), 
               leaf_rotation=90, 
               color_threshold=cut_height,
               above_threshold_color="black")
plt.axhline(y=cut_height, color='grey', linestyle='--', linewidth=1.5)
plt.title("Hierarchical Clustering Dendrogram")
plt.xlabel("PDB ID")
plt.ylabel("Distance")
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "RMSD_cluster_dendrogram.png"), dpi=300, bbox_inches="tight")

##############################################################################

#############################
### RMSD Heatmap Plotting ###
#############################

# Convert df_wide to long format
mean_distance_heatmap = df_wide.melt(
    id_vars='pdb_id',
    var_name='comparison',
    value_name='mean_distance'
)

# Specifying plot order based on RMSD cluster order stored in pdb_id_labels
mean_distance_heatmap['pdb_id'] = pd.Categorical(mean_distance_heatmap['pdb_id'], categories=pdb_id_labels, ordered=True)
mean_distance_heatmap['comparison'] = pd.Categorical(mean_distance_heatmap['comparison'], categories=pdb_id_labels, ordered=True)

# Applying the order to the DataFrame (so it is saved in the correct order in the CSV file)
mean_distance_heatmap = mean_distance_heatmap.sort_values(['pdb_id', 'comparison']).reset_index(drop=True)

# Saving the mean distance heatmap data to a CSV file
mean_distance_heatmap.to_csv(os.path.join(data_path, "RMSD_heatmap.csv"), index=False)

# Heatmap plot
p = (
    ggplot(mean_distance_heatmap, aes(x='comparison', y='pdb_id', fill='mean_distance'))
    + geom_tile(color='white')
    + scale_fill_gradient(low='white', high='red', name='Mean RMSD (Å)', limits=(0, np.max(mean_distance_heatmap['mean_distance']) * 1.05))
    + labs(title='Mean Pairwise RMSD Heatmap', x='PDB ID', y='PDB ID')
    + theme_minimal()
    + theme(
        panel_border = element_rect(color="black", size=1.5),
        panel_grid_major = element_blank(),
        panel_grid_minor = element_blank(),
        plot_title = element_text(size=20, face="bold", ha="center"),
        axis_title = element_text(size=25, face="bold"),
        axis_text_x = element_text(size=12, colour="black", angle=90, vjust=0.5, hjust=1),
        axis_text_y = element_text(size=12, colour="black"),
        legend_title = element_text(size=20, face="bold"),
        legend_text = element_text(size=16)
    )
)

# Saving plot
p.save(os.path.join(output_dir, "RMSD_heatmap.png"), width=18, height=15, dpi=300)

