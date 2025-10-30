import pandas as pd
import numpy as np
import os
import sys
from scipy.signal import find_peaks
from plotnine import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

#############
### Setup ###
#############

# Creating output_dir
output_dir = os.path.join("Output", "stable_regions")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)
    
# Importing data
filtered_df = pd.read_csv(os.path.join("Output", "thermodynamics", "foldx_stability_filtered.csv"))
mean_dG_per_res_by_PDB = pd.read_csv(os.path.join("Output", "thermodynamics", "mean_deltaG_per_residue.csv"))

# Read window_size argument 
if len(sys.argv) > 1: 
    try: 
        window_size = int(sys.argv[1]) 
        print(f"Window size: {window_size}") 
    except ValueError: 
        print("Warning: invalid window size provided — defaulting to 3.") 
        window_size = 3 
else: 
    window_size = 3 
    print("No window size provided, defaulting to 3")

################
### Analysis ###
################

# Calculating the rolling average 
mean_dG_per_res_by_PDB["sliding_window"] = (
    mean_dG_per_res_by_PDB
    .groupby("pdb_id")["mean_energy"]
    .transform(lambda x: x.rolling(window = window_size, center = True, min_periods = window_size).mean())
)

# Calculating the mean dG per residue across all PDBs
mean_dG_by_pos = (
    mean_dG_per_res_by_PDB
    .groupby("Pos", as_index = False)["sliding_window"]
    .mean()
    .dropna()
)

# Calculating the mean stability
mean_stability = mean_dG_per_res_by_PDB["sliding_window"].mean()
#mean_stability = mean_dG_per_res_by_PDB["mean_energy"].mean()

#######################################
### Finding Local Minima and Maxima ###
#######################################

x = mean_dG_by_pos["Pos"]
y = mean_dG_by_pos["sliding_window"]

# Find local maxima (peaks)
peaks, peak_props = find_peaks(y, prominence = 0.01)  
# Find local minima (by inverting the signal)
pits, pit_props = find_peaks(-y, prominence = 0.01)

# Combine into a single dataframe of turning points
tp_df = pd.concat([
    pd.DataFrame({
        "Pos": x.iloc[peaks].values,
        "type": "peak",
        "prominence": peak_props["prominences"]
    }),
    pd.DataFrame({
        "Pos": x.iloc[pits].values,
        "type": "pit",
        "prominence": pit_props["prominences"]
    })
])

# Sort by position
tp_df = tp_df.sort_values("Pos").reset_index(drop = True)

# filtering turning points by prominence
threshold = 0.15  
tp_df = tp_df[tp_df["prominence"] > threshold]

# Adding in sliding_window to tp_df
tp_df = pd.merge(tp_df, mean_dG_by_pos, on = "Pos", how = "right").drop(columns = ["prominence"])

# Converting NaN to strings
tp_df["type"] = tp_df["type"].fillna("NA")

###############################
### Defining stable regions ###
###############################

# Selecting peaks above the mean line
significant_peaks = tp_df.loc[(tp_df["type"] == "peak") & (tp_df["sliding_window"] > mean_stability), "Pos"].tolist()

# Initialise variables
stabilisers = []
region_list = []
region = 1
stabiliser_count = 0

# Read min_stable_region_size argument
if len(sys.argv) > 2: 
    try: 
        min_region_size = int(sys.argv[2]) 
        print(f"Window size: {min_region_size}") 
    except ValueError: 
        print("Warning: invalid min region size provided — defaulting to 0") 
        min_region_size = 0 
else: 
    min_region_size = 0 
    print("No min region size provided, defaulting to 0")


# Convert to list for faster iteration
positions = mean_dG_by_pos["Pos"].tolist()
energies = mean_dG_by_pos["sliding_window"].tolist()

# Looping through each residue position
for i in range(len(positions)):
    pos = positions[i]
    energy = energies[i]

    # Check if the current position is stabilising
    if (energy < mean_stability):
        stabilisers.append(pos)
        region_list.append(region)
        stabiliser_count += 1

    # If a significant peak is hit
    if pos in significant_peaks:
        # Remove regions smaller than min_region_size
        if stabiliser_count < min_region_size:
            stabilisers = [s for s, r in zip(stabilisers, region_list) if r != region]
            region_list = [r for r in region_list if r != region]

        # Increment region for next potential stable region
        region += 1
        stabiliser_count = 0

# Handle the last region in case the sequence ends without a peak
if stabiliser_count < min_region_size:
    stabilisers = [s for s, r in zip(stabilisers, region_list) if r != region]
    region_list = [r for r in region_list if r != region]

# Create a DataFrame of stabilising positions and their assigned region
stabilising_positions = pd.DataFrame({
    "Pos": stabilisers,
    "region": region_list
})

# Get the unique region numbers, sorted in ascending order
unique_regions = sorted(stabilising_positions["region"].unique())

# Make sure regions start at 1 and increase by 1 integer between each region
region_map = {old: new for new, old in enumerate(unique_regions, start = 1)}
stabilising_positions["region"] = stabilising_positions["region"].map(region_map)

# Group stabilising positions by region and compute min and max Pos in each region
stable_regions = (
    stabilising_positions
    .groupby("region", as_index = False)
    .agg(min = ("Pos", "min"), max = ("Pos", "max"))
)

# Create a string describing the residue range in each region
stable_regions["residues"] = stable_regions["min"].astype(str) + "-" + stable_regions["max"].astype(str)

# Convert 'region' to a categorical type with sorted order
stable_regions["region"] = stable_regions["region"].astype(str)
# Sort the regions numerically, but keep as strings
numeric_sorted_regions = [str(i) for i in sorted(stable_regions["region"].astype(int))]
# Assign categorical with numeric order
stable_regions["region"] = pd.Categorical(
    stable_regions["region"],
    categories=numeric_sorted_regions,
    ordered=True
)

# Saving stable_regions
stable_regions_minimal = stable_regions[["region", "residues"]]
stable_regions_minimal.to_csv(os.path.join(output_dir, "stable_regions.csv"), index = False)

# Using matplotlib colour palette tab20
n_regions = 20
cmap = plt.get_cmap("tab20")  # qualitative palette
# get first n_regions colours
colours = [cmap(i) for i in range(n_regions)]
# convert to hex
colours = [mcolors.to_hex(c) for c in colours]

# Extend the color vector to match the length of `region` if needed
region_colours = np.tile(colours, int(np.ceil(len(stable_regions) / len(colours))))[:len(stable_regions)]

# Creating a new column to adjust the width of the rectangles in the plot
stable_regions["min_adj"] = stable_regions["min"] - 0.5
stable_regions["max_adj"] = stable_regions["max"] + 0.5

#######################################
### Visualising Stabilising Regions ###
#######################################

# Plotting turning points
p = (
    ggplot(tp_df, aes(x = "Pos", y = "sliding_window")) 
    + geom_rect(stable_regions, aes(xmin="min_adj", xmax="max_adj", ymin=-np.inf, ymax=np.inf, fill="region"),
        alpha=0.75, inherit_aes=False)
    + geom_point(aes(colour = "type"), size = 2) 
    + geom_line(size = 1) 
    + scale_x_continuous(breaks = range(0, 9999, 10))
    + scale_y_continuous(expand = (0.05, 0.05))
    + geom_hline(yintercept = mean_stability, linetype = "dashed", size = 0.75, colour = "black")
    + scale_fill_manual(values = region_colours)
    + scale_color_manual(values = {"pit": "red", "peak": "blue", "NA": "black"},
                         breaks = ["pit", "peak", "NA"],
                         labels = ["Pit", "Peak", "NA"]) 
    + labs(
        title = f"Turning Points: Sliding Window of {window_size}",
        y = r'$\mathbf{\Delta G^{\circ}\ per\ residue\ (kcal\cdot mol^{-1})}$',
        x = "Residue Number",
        colour = "Type",
        fill = "Region"
    )
    + theme_minimal()
    + theme(
        panel_border = element_rect(size = 1, colour = "black"),
        panel_grid_major = element_line(colour = "grey", size = 0.25, linetype = "dashed"),
        panel_grid_minor = element_line(colour = "grey", size = 0.25, linetype = "dashed"),
        plot_title = element_text(size = 14, colour = "black", face = "bold", ha = "center"),
        axis_title = element_text(size = 12, colour = "black", face = "bold"),
        axis_text = element_text(size = 10, colour = "black", face = "bold"),
        legend_title = element_text(size = 12, colour = "black", face = "bold"),
        legend_text = element_text(size = 10, colour = "black")
    )
)
# Saving turning points plot
p.save(os.path.join(output_dir, "turning_points_and_regions.png"), height = 5, width = 10, dpi = 300)

# Creating a group variable based on pdb_id to include line breaks if residue numbers are not continuous
mean_dG_per_res_by_PDB = mean_dG_per_res_by_PDB.sort_values(['pdb_id', 'Pos'])
mean_dG_per_res_by_PDB['segment'] = mean_dG_per_res_by_PDB.groupby('pdb_id')['Pos'].diff().ne(1).cumsum()
mean_dG_per_res_by_PDB['interaction'] = mean_dG_per_res_by_PDB['pdb_id'].astype(str) + "_" + mean_dG_per_res_by_PDB['segment'].astype(str)

# Plotting stable regions and sliding window data
p = (
    ggplot(mean_dG_per_res_by_PDB, aes(x = "Pos")) 
    + geom_rect(stable_regions, aes(xmin="min_adj", xmax="max_adj", ymin=-np.inf, ymax=np.inf, fill="region"),
        alpha=0.75, inherit_aes=False)
    + geom_hline(yintercept = 0, colour = "black", size = 1, linetype = "solid", alpha = 0.9) 
    + geom_line(aes(y = "sliding_window", group = "interaction"),
                size = 0.5, alpha = 0.25) 
    + scale_x_continuous(breaks = range(0, 9999, 10))
    + scale_y_continuous(expand = (0.05, 0.05))
    + scale_fill_manual(values = region_colours)
    + scale_color_manual(values = {"pit": "red", "peak": "blue", "NA": "black"},
                         breaks = ["pit", "peak", "NA"],
                         labels = ["Pit", "Peak", "NA"]) 
    + labs(
        title = f"Rolling Average Using Sliding Window of Size {window_size}",
        y = r'$\mathbf{\Delta G^{\circ}\ per\ residue\ (kcal\cdot mol^{-1})}$',
        x = "Residue Number",
        colour = "Type",
        fill = "Region"
    )
    + theme_minimal()
    + theme(
        panel_border = element_rect(size = 1, colour = "black"),
        panel_grid_major = element_line(colour = "grey", size = 0.25, linetype = "dashed"),
        panel_grid_minor = element_line(colour = "grey", size = 0.25, linetype = "dashed"),
        plot_title = element_text(size = 14, colour = "black", face = "bold", ha = "center"),
        axis_title = element_text(size = 12, colour = "black", face = "bold"),
        axis_text = element_text(size = 10, colour = "black", face = "bold"),
        legend_title = element_text(size = 12, colour = "black", face = "bold"),
        legend_text = element_text(size = 10, colour = "black")
    )
)
p.save(os.path.join(output_dir, "sliding_window_and_regions.png"), height = 5, width = 10, dpi = 300)

# Plotting stable regions and mean deltaG per residue for each PDB
p = (
    ggplot(mean_dG_per_res_by_PDB, aes(x = "Pos")) 
    + geom_rect(stable_regions, aes(xmin="min_adj", xmax="max_adj", ymin=-np.inf, ymax=np.inf, fill="region"),
        alpha=0.75, inherit_aes=False)
    + geom_hline(yintercept = 0, colour = "black", size = 1, linetype = "solid", alpha = 0.9) 
    + geom_line(aes(y = "mean_energy", group = "interaction"),
                size = 0.5, alpha = 0.25) 
    + scale_x_continuous(breaks = range(0, 9999, 10))
    + scale_y_continuous(expand = (0.05, 0.05))
    + scale_fill_manual(values = region_colours)
    + scale_color_manual(values = {"pit": "red", "peak": "blue", "NA": "black"},
                         breaks = ["pit", "peak", "NA"],
                         labels = ["Pit", "Peak", "NA"]) 
    + labs(
        title = "",
        y = r'$\mathbf{\Delta G^{\circ}\ per\ residue\ (kcal\cdot mol^{-1})}$',
        x = "Residue Number",
        colour = "Type",
        fill = "Region"
    )
    + theme_minimal()
    + theme(
        panel_border = element_rect(size = 1, colour = "black"),
        panel_grid_major = element_line(colour = "grey", size = 0.25, linetype = "dashed"),
        panel_grid_minor = element_line(colour = "grey", size = 0.25, linetype = "dashed"),
        plot_title = element_text(size = 14, colour = "black", face = "bold", ha = "center"),
        axis_title = element_text(size = 12, colour = "black", face = "bold"),
        axis_text = element_text(size = 10, colour = "black", face = "bold"),
        legend_title = element_text(size = 12, colour = "black", face = "bold"),
        legend_text = element_text(size = 10, colour = "black")
    )
)
p.save(os.path.join(output_dir, "mean_delta_G_and_regions.png"), height = 5, width = 10, dpi = 300)
