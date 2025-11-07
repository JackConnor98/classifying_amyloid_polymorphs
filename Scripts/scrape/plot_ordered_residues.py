import os
import pandas as pd
import numpy as np
import math
from plotnine import *
from Bio.PDB import PDBParser

print("Running: plot_ordered_residues.py")

# =============================================================================
# Load and clean data
# =============================================================================
data = pd.read_csv(os.path.join("Output", "selected_pdbs_metadata.csv"))

# Remove NMR structures
data = data[data["Method"] != "ssNMR"]

# Keep relevant columns
residues = data[["PDB ID", "Residues Ordered"]].copy()
residues.columns = ["pdb", "ordered_residues"]

# Split by commas (indicating discontinuities)
residues["ordered_residues"] = residues["ordered_residues"].str.split(",")

# Flatten the lists into long form
residues = residues.explode("ordered_residues").reset_index(drop=True)

# Clean messy characters
residues["ordered_residues"] = (
    residues["ordered_residues"]
    .astype(str)
    .str.replace(r'[c()"()]', "", regex=True)
    .str.strip()
)

# Split comma-separated regions
residues["ordered_residues"] = residues["ordered_residues"].str.split(",")

# Expand into separate rows
residues = residues.explode("ordered_residues").reset_index(drop=True)

# Clean up stray characters
residues["ordered_residues"] = (
    residues["ordered_residues"]
    .astype(str)
    .str.replace(r'[c()"()]', "", regex=True)
    .str.strip()
)

# Split each stretch (e.g. 20-96) into start and end
split_cols = residues["ordered_residues"].str.split("-", expand=True)
split_cols.columns = ["start", "end"]

# Convert to numeric (ignore errors from NaN or empty values)
split_cols = split_cols.apply(pd.to_numeric, errors="coerce")

# Add a group ID for each continuous region per PDB
residues = pd.concat([residues[["pdb"]], split_cols], axis=1)
residues["group"] = residues.groupby("pdb").cumcount() + 1

# Expand start/end into separate rows
residues_long = pd.melt(
    residues,
    id_vars=["pdb", "group"],
    value_vars=["start", "end"],
    var_name="boundary",
    value_name="ordered_residues"
).drop(columns=["boundary"])

# Convert to numeric
residues_long["ordered_residues"] = pd.to_numeric(residues_long["ordered_residues"], errors="coerce")
residues_long = residues_long.dropna(subset=["ordered_residues"]).reset_index(drop=True)


residues_long["pdb_group"] = residues_long["pdb"].astype(str) + "_" + residues_long["group"].astype(str)
  
# =============================================================================
# Plotting
# =============================================================================
num_pdbs = len(data)
protein_names = data["PDB ID"].unique()

# Set base width per label (tweak this depending on how dense your labels are)
width_per_label = 0.2  # inches per PDB label
min_width = 8  # minimum plot width

# Calculate final width
plot_width = max(min_width, num_pdbs * width_per_label)
plot_width = min(plot_width, 20)  # cap maximum width to 25 inches

# Set base font size and minimum font size
base_font_size = 14
min_font_size = 8

# Scale font size down as number of labels increases
# The more labels, the smaller the text
x_axis_font_size = max(min_font_size, base_font_size - 0.05*num_pdbs)

# Create an interaction column
residues_long["pdb_group_interaction"] = residues_long["pdb"].astype(str) + "_" + residues_long["group"].astype(str)

### Creating protein specific annotations ###
protein_type = data["Protein"].str.lower()  # lowercase for easier matching

# Initialize annotation and label DataFrames
annotations = pd.DataFrame(columns=["xmin","xmax","ymin","ymax","fill"])
labels = pd.DataFrame(columns=["x","y","label"])

# a-Synuclein case
if protein_type.str.contains("synuclein").any():
    annotations = pd.DataFrame({
        "xmin": [-math.inf, num_pdbs + 1, -math.inf, num_pdbs + 1, -math.inf, num_pdbs + 1],
        "xmax": [math.inf, num_pdbs + 7, math.inf, num_pdbs + 7, math.inf, num_pdbs + 7],
        "ymin": [1, 13, 61, 73, 96, 113],
        "ymax": [60, 23, 95, 83, 140, 123],
        "fill": ["darkblue","darkblue","red","red","green","green"]
    })
    labels = pd.DataFrame({
        "x": [num_pdbs + 4, num_pdbs + 4, num_pdbs + 4],
        "y": [18, 78, 118],
        "label": ["N-Term","NAC","C-Term"]
    })

# Amyloid-beta case
elif protein_type.str.contains("amyloid-").any():
    annotations = pd.DataFrame({
        "xmin": [-math.inf, num_pdbs + 0.75, -math.inf, num_pdbs + 0.75, -math.inf, num_pdbs + 0.75, -math.inf, num_pdbs + 0.75, -math.inf, num_pdbs + 0.75],
        "xmax": [math.inf, num_pdbs + 2.75, math.inf, num_pdbs + 2.75, math.inf, num_pdbs + 2.75, math.inf, num_pdbs + 2.75, math.inf, num_pdbs + 2.75],
        "ymin": [1, 6, 17, 17, 24, 23.5, 28, 29.5, 36, 37],
        "ymax": [16, 10, 21, 21, 27, 27.5, 35, 33.5, 42, 41],
        "fill": ["darkblue","darkblue","black","black","red","red","orange","orange","green","green"]
    })
    labels = pd.DataFrame({
        "x": [num_pdbs + 1.75]*5,
        "y": [8, 19, 25.5, 31.5, 39],
        "label": ["N-Term","Hydrophobic\nCore","Turn\nRegion","Hydrophobic\nRegion","C-Term"]
    })

# Tau case
elif protein_type.str.contains("tau").any():
    annotations = pd.DataFrame({
        "xmin": [-math.inf, num_pdbs + 1, -math.inf, num_pdbs + 1, -math.inf, num_pdbs + 1, -math.inf, num_pdbs + 1],
        "xmax": [math.inf, num_pdbs + 5, math.inf, num_pdbs + 5, math.inf, num_pdbs + 5, math.inf, num_pdbs + 5],
        "ymin": [243, 253, 274, 284, 305, 315, 336, 346.5],
        "ymax": [273, 263, 304, 294, 335, 325, 367, 356.5],
        "fill": ["red","red","green","green","orange","orange","darkblue","darkblue"]
    })
    labels = pd.DataFrame({
        "x": [num_pdbs + 3]*4,
        "y": [258, 289, 320, 351.5],
        "label": ["R1","R2","R3","R4"]
    })

# Create the plot
p = (
    ggplot() 

    + geom_rect(
        aes(xmin="xmin", xmax="xmax", ymin="ymin", ymax="ymax", fill="fill"),
        alpha=0.5,
        colour="black",
        data=annotations
    )
    + geom_text(
        aes(x="x", y="y", label="label"),
        data=labels,
        size=8,
        fontweight="bold"
    )

    + geom_point(residues_long, aes(x="pdb", y="ordered_residues"), size=3)
    + geom_line(residues_long, aes(x="pdb", y="ordered_residues", group="pdb"), size=1, linetype="dotted", colour="black")
    + geom_line(residues_long, aes(x="pdb", y="ordered_residues", group="pdb_group_interaction"), size=2) 
    + geom_line(residues_long, aes(x="pdb", y="ordered_residues", group = "pdb"), size = 1, linetype = "dotted", colour = "black") 
    + geom_line(residues_long, aes(x="pdb", y="ordered_residues", group = "pdb_group"), size = 2)
    + labs(y="Residue Position", x="PDB") 
    + scale_x_discrete() 
    + scale_y_continuous(breaks=range(0, 10000, 10)) 
    + theme_bw() 
    + theme(
        panel_border=element_rect(colour="black", fill=None, size=1.5),
        axis_text_y=element_text(size=14, colour="black"),
        axis_text_x=element_text(size=x_axis_font_size, colour="black", angle=90, va="center"),
        axis_title=element_text(size=20, weight="bold"),
        axis_line=element_line(colour="black", size=1.5),
        legend_position="none"
    )
)

# Save plot
p.save("Output/Ordered_Residues.png", width=plot_width, height=8, dpi=300)
print("Plot saved to Output/Ordered_Residues.png")









