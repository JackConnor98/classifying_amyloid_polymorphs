# This anaylsis was inspired by Errico et al 2025 https://doi.org/10.1042/bcj20240602
import os
import sys
import numpy as np
import pandas as pd
from Bio.PDB import *
from plotnine import *
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

################################
### Reading in sys arguments ###
################################

if len(sys.argv) > 1: 
    try: 
        min_length = int(sys.argv[1]) 
        print(f"Min B-strand length: {min_length}") 
    except ValueError: 
        print("Warning: invalid minimum length provided — defaulting to 4") 
        min_length = 4 
else: 
    min_length = 4 
    print("No minimum length provided, defaulting to 4")

############################
### Creating Directories ###
############################

output_dir = os.path.join("Output", "b_sheets")
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

ramachandran_dir = os.path.join(output_dir, "ramachandran_plots")
if not os.path.exists(ramachandran_dir):
    os.makedirs(ramachandran_dir)
    
plot_dir = os.path.join(output_dir, "plots")
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)

# Specifying range for B-sheet conformation
# Reference: 
# B-Turns and their distortions: a proposed new nomenclature (Wilmot and Thornton 1990) - Figure 2
# https://doi.org/10.1093/protein/3.6.479

psi_min = 60
psi_max = 180
psi_min_extended = -180
psi_max_extended = -150
phi_min = -180
phi_max = -90


# Locate all pdbs in folder that need to be looped through
pdb_path = os.path.join("Output", "PDBs", "unique_chains")

# Get a list of all .pdb files in the folder
pdb_files = [
    os.path.join(pdb_path, file)
    for file in os.listdir(pdb_path)
    if file.endswith(".pdb")
]
# Getting the base names of the PDB files
pdb_names = [
    os.path.splitext(os.path.basename(f))[0]
    for f in pdb_files
]
# Getting the list of high resolution PDBs
high_resolution_PDBs = pd.read_csv(os.path.join("Output", "Validation", "high_resolution_pdb_ids.csv"))

# Removing the PDBs that are not found in high_resolution_PDBs
pdb_files_filtered = [
    f for f, name in zip(pdb_files, pdb_names)
    if name in high_resolution_PDBs["pdb_id"].values
]

# DataFrames for storing all torsion data and proportions
all_torsion_data = pd.DataFrame()
b_strand_data = pd.DataFrame()

parser = PDBParser(QUIET=True)

for i, pdb_path in enumerate(pdb_files_filtered):
    pdb_name = pdb_names[i]
    print(f"Analysing - {pdb_name}")

    # Parse the PDB file into a hierarchical structure (Model → Chain → Residue → Atom)
    structure = parser.get_structure(pdb_name, pdb_path)

    # List to store torsion angles
    torsion_records = []

    for model in structure:
        for chain in model:
            try:
                # Build all continuous polypeptide segments within this chain
                polypeptides = Polypeptide.PPBuilder().build_peptides(chain)
            except Exception:
                # If something goes wrong (e.g. invalid chain), skip this one
                continue
            
            for poly in polypeptides:
                phi_psi = poly.get_phi_psi_list()
                residues = poly

                for residue, (phi, psi) in zip(residues, phi_psi):
                    # Extract residue information
                    resname = residue.get_resname()
                    resnum = residue.get_id()[1]
                    chain_id = chain.get_id()

                    # Skip unknown residues or missing torsions
                    if resname == "UNK" or (phi is None and psi is None):
                        continue

                    torsion_records.append({
                        "pos": resnum,
                        "chain": chain_id,
                        "residue": resname,
                        "phi": np.degrees(phi) if phi is not None else np.nan,
                        "psi": np.degrees(psi) if psi is not None else np.nan,
                        "pdb_id": pdb_name
                    })

    # Convert to DataFrame
    torsion_df = pd.DataFrame(torsion_records)
    if torsion_df.empty:
        print(f"Skipping {pdb_name}: PDB file is empty or contains no atomic coordinates.")
        continue

    # Append to global torsion dataset
    all_torsion_data = pd.concat([all_torsion_data, torsion_df], ignore_index=True)

    # Plotting Ramachandran
    # p = (
    #     ggplot(torsion_df, aes(x = "phi", y = "psi"))
    #     + geom_rect(xmin = phi_min, xmax = phi_max,
    #                 ymin = psi_min, ymax = psi_max,
    #                 fill = "skyblue", alpha = 0.01)
    #     + geom_rect(xmin = phi_min, xmax = phi_max,
    #                 ymin = psi_min_extended, ymax = psi_max_extended,
    #                 fill = "skyblue", alpha = 0.01)
    #     + geom_point(colour = "black", alpha = 0.5, size = 2)
    #     + scale_x_continuous(breaks = range(-180, 181, 30), limits = [-180, 180], expand = [0, 0])
    #     + scale_y_continuous(breaks = range(-180, 181, 30), limits = [-180, 180], expand = [0, 0])
    #     + labs(
    #         title = f"Ramachandran Plot - {pdb_name}",
    #         x = "Phi (φ) Angles", 
    #         y = "Psi (ψ) Angles")
    #     + theme_bw()
    #     + theme(
    #         panel_grid_major = element_line(colour = "grey", size = 0.5),
    #         plot_title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
    #         axis_title = element_text(size = 16, colour = "black", face = "bold"),
    #         axis_text = element_text(size = 12, colour = "black")
    #     )
    # )

    # p.save(os.path.join(ramachandran_dir, f"{pdb_name}.png"), height = 8, width = 8)

    # --- Filter residues by torsion angles ---
    b_strands = torsion_df[
        ((torsion_df["psi"] > psi_min) & (torsion_df["psi"] < psi_max) &
        (torsion_df["phi"] > phi_min) & (torsion_df["phi"] < phi_max)) |
        ((torsion_df["psi"] > psi_min_extended) & (torsion_df["psi"] < psi_max_extended) &
        (torsion_df["phi"] > phi_min) & (torsion_df["phi"] < phi_max))
    ]

    # Append to global proportion DataFrame
    b_strand_data = pd.concat([b_strand_data, b_strands], ignore_index=True)

# Saving data
all_torsion_data.to_csv(os.path.join(output_dir, "all_torsion_data.csv"), index = False)
b_strand_data.to_csv(os.path.join(output_dir, "b_strands.csv"), index = False)

###############################################
### Ramachandran Plot for all PDBs combined ###
###############################################

# p = (
#     ggplot(all_torsion_data, aes(x = "phi", y = "psi"))
#     + geom_rect(xmin = phi_min, xmax = phi_max,
#                 ymin = psi_min, ymax = psi_max,
#                 fill = "skyblue", alpha = 0.01)
#     + geom_rect(xmin = phi_min, xmax = phi_max,
#                 ymin = psi_min_extended, ymax = psi_max_extended,
#                 fill = "skyblue", alpha = 0.01)
#     + geom_point(colour = "black", alpha = 0.5, size = 2)
#     + scale_x_continuous(breaks = range(-180, 181, 30), limits = [-180, 180], expand = [0, 0])
#     + scale_y_continuous(breaks = range(-180, 181, 30), limits = [-180, 180], expand = [0, 0])
#     + labs(
#         x = "Phi (φ) Angles", 
#         y = "Psi (ψ) Angles")
#     + theme_bw()
#     + theme(
#         panel_grid_major = element_line(colour = "grey", size = 0.5),
#         plot_title = element_text(size = 20, colour = "black", face = "bold", hjust = 0.5),
#         axis_title = element_text(size = 16, colour = "black", face = "bold"),
#         axis_text = element_text(size = 12, colour = "black")
#     )
# )

# p.save(os.path.join(output_dir, f"all_pdbs_ramachandran.png"), height = 8, width = 8)

##########################################################
### Identifying B-sheets longer than a given threshold ###
##########################################################

# Sort by PDB and residue position
b_strand_data = b_strand_data.sort_values(["pdb_id", "pos"]).reset_index(drop=True)

# Identify consecutive sequences using a group id
def mark_consecutive(group):
    # compute the difference between consecutive positions
    diff = group["pos"].diff().fillna(1)
    # increment group id when the difference is not 1
    group_id = (diff != 1).cumsum()
    group["seq_group"] = group_id
    # keep only sequences of length >= min_length
    keep = group.groupby("seq_group")["pos"].transform("count") >= min_length
    return group[keep]

if min_length <= 1:
    sig_b_strands = b_strand_data
    
else:
    # Apply per pdb_id
    sig_b_strands = b_strand_data.groupby("pdb_id", group_keys = False).apply(mark_consecutive, include_groups = True)

    # Drop the helper column
    sig_b_strands = sig_b_strands.drop(columns = "seq_group")

# Saving significant data
sig_b_strands.to_csv(os.path.join(output_dir, "significant_beta_strands.csv"), index = False)


##########################################################################
### Calculating the total b-strand proportion across all unique chains ###
##########################################################################

total_PDBs = len(pdb_files_filtered)

# Counting the total occurences of each residue
pos_counts_total = all_torsion_data["pos"].value_counts().reset_index()
pos_counts_total.columns = ["pos", "total_count"]
pos_counts_total = pos_counts_total.sort_values(by = "pos", ascending = True).reset_index(drop = True)

# Counting the num. times a residue is found in a b-strand
pos_counts_b_strand = b_strand_data["pos"].value_counts().reset_index()
pos_counts_b_strand.columns = ["pos", "strand_count"]
pos_counts_b_strand = pos_counts_b_strand.sort_values(by = "pos", ascending = True).reset_index(drop = True)

# Counting the num. times a residue is found in a significant b-strand
pos_counts_signif = sig_b_strands["pos"].value_counts().reset_index()
pos_counts_signif.columns = ["pos", "signif_count"]
pos_counts_signif = pos_counts_signif.sort_values(by = "pos", ascending = True).reset_index(drop = True)

count_df = (
    pos_counts_total
    .merge(pos_counts_b_strand, on = "pos", how = "outer")
    .merge(pos_counts_signif, on = "pos", how = "outer")
    .fillna(0)
)

count_df["norm_strand_count"] = count_df["strand_count"] / count_df["total_count"]
count_df["norm_signif_count"] = count_df["signif_count"] / count_df["total_count"]

count_df.to_csv(os.path.join(output_dir, "residue_counts.csv"), index = False)


### Plotting ###

# Reshape for plotting
count_long = count_df.melt(id_vars = ["pos"],
                           value_vars = ["total_count", "strand_count", "signif_count"],
                           var_name = "type",
                           value_name = "count")

# Setting the legend order
count_long["type"] = pd.Categorical(
    count_long["type"],
    categories = ["total_count", "strand_count", "signif_count"],
    ordered = True
)

# Make it pretty
p = (
    ggplot(count_long, aes(x = "pos", y = "count", colour = "type"))
    + geom_line(size = 0.75)
    + scale_colour_manual(
        values = {
            "total_count": "black",
            "strand_count": "blue",
            "signif_count": "red"
        },
        name = "Count Type",
        labels = {
            "total_count": "Total",
            "strand_count": "β-strand",
            "signif_count": "Significant β-strand"
        }
    )
    + scale_x_continuous(breaks = range(0, 9999, 10))
    + scale_y_continuous(breaks = range(0, 9999, 10))
    + labs(
        y = "Count",
        x = "Residue Number"
    )
    + theme_classic()
    + theme(
        panel_border = element_rect(size = 1),
        panel_grid = element_line(colour = "grey", size = 0.25, linetype = "dashed"),
        axis_title = element_text(size = 12, weight = "bold"),
        axis_text = element_text(size = 10),
        legend_title = element_text(size = 12, weight = "bold"),
        legend_text = element_text(size = 10)
    )
)
p.save(os.path.join(plot_dir, "residue_counts.png"), height = 5, width = 10, dpi = 300)


### Plotting B-strands and stable regions to see if they overlap ###

# Importing stabilising regions
stable_regions = pd.read_csv(os.path.join("Output", "stable_regions", "stable_regions.csv"))
# Separating residues into start and end
stable_regions[["min", "max"]] = stable_regions["residues"].str.split("-", expand=True).astype(float)
# Creating a new column to adjust the width of the rectangles in the plot
stable_regions["min_adj"] = stable_regions["min"] - 0.5
stable_regions["max_adj"] = stable_regions["max"] + 0.5
# Making region categorical
stable_regions["region"] = pd.Categorical(stable_regions["region"])


# Using matplotlib colour palette tab20
n_regions = 20
cmap = plt.get_cmap("tab20")  # qualitative palette
# get first n_regions colours
colours = [cmap(i) for i in range(n_regions)]
# convert to hex
colours = [mcolors.to_hex(c) for c in colours]
# Extend the color vector to match the length of `region` if needed
region_colours = np.tile(colours, int(np.ceil(len(stable_regions) / len(colours))))[:len(stable_regions)]

# Removing total count
b_strand_count_long = count_long[count_long["type"] != "total_count"]
b_strand_count_long["type"] = pd.Categorical(
    b_strand_count_long["type"],
    categories = ["strand_count", "signif_count"],
    ordered = True
)

# Plotting
p = (
    ggplot()
    + geom_line(b_strand_count_long, aes(x = "pos", y = "count", linetype = "type"), size = 0.75, colour = "black")
    + geom_rect(stable_regions, aes(xmin="min_adj", xmax="max_adj", ymin=-np.inf, ymax=np.inf, fill="region"),
        alpha=0.75, inherit_aes=False)
    + scale_linetype_manual(
        values = {
            "strand_count": "solid",
            "signif_count": "dashed"
        },
        labels = {
            "strand_count": "β-strand",
            "signif_count": "Significant β-strand"
        },
        name = "Count Type"
    )
    + scale_fill_manual(values = region_colours)
    + scale_x_continuous(breaks = range(0, 9999, 10))
    + scale_y_continuous(breaks = range(0, 9999, 10))
    + labs(
        y = "Count",
        x = "Residue Number"
    )
    + theme_classic()
    + theme(
        panel_border = element_rect(size = 1),
        panel_grid = element_line(colour = "grey", size = 0.25, linetype = "dashed"),
        axis_title = element_text(size = 12, weight = "bold"),
        axis_text = element_text(size = 10),
        legend_title = element_text(size = 12, weight = "bold"),
        legend_text = element_text(size = 10)
    )
)
p.save(os.path.join(plot_dir, "strands_vs_regions.png"), height = 5, width = 10, dpi = 300)