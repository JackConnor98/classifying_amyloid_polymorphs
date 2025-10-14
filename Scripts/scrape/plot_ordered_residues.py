import os
import pandas as pd
import numpy as np
from plotnine import *
from Bio.PDB import PDBParser

print("Running: plot_ordered_residues.py")

# =============================================================================
# Load and clean data
# =============================================================================
data = pd.read_csv(os.path.join("Output", "selected_pdbs_metadata.txt"), sep="\t")

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

p = (
    ggplot(residues_long, aes(x="pdb", y="ordered_residues")) 
    + geom_point(size=3)
    + geom_line(aes(group="pdb"), size=1, linetype="dotted", colour="black")
    + geom_line(aes(y="ordered_residues", group="pdb"), size=2) 
    + geom_line(aes(group = "pdb"), size = 1, linetype = "dotted", colour = "black") 
    + geom_line(aes(group = "pdb_group"), size = 2)
    + labs(y="Residue Position", x="PDB") 
    + scale_x_discrete() 
    + scale_y_continuous(breaks=range(0, 10000, 10)) 
    + theme_bw() 
    + theme(
        panel_border=element_rect(colour="black", fill=None, size=1.5),
        plot_title=element_text(size=20, weight="bold", colour="black", ha="center"),
        axis_text_y=element_text(size=14, colour="black", weight="bold"),
        axis_text_x=element_text(size=14, colour="black", weight="bold", angle=90, ha="right", va="center"),
        axis_title=element_text(size=20, weight="bold"),
        axis_line=element_line(colour="black", size=1.5)
    ) 
)

# Save
p.save("Output/Ordered_Residues.png", width=20, height=8, dpi=300)
print("Plot saved to Output/Ordered_Residues.png")









