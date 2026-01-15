from png_functions import process_pdb_directory, set_coloring, create_grouped_png_grid, create_png_grid
import pymol
from pymol import cmd
import os 
import pandas as pd
import sys
from PIL import Image, ImageDraw, ImageFont
import math
from itertools import cycle
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

#############################################################################################################################
 # Reading in system arguments

if len(sys.argv) > 1: 
    try: 
        png_colouring = int(sys.argv[1]) 
    except ValueError: 
        print("Warning: invalid value provided for png_colouring — defaulting to 0") 
        png_colouring = 0 
else: 
    png_colouring = 0 
    print("No value for png_colouring provided, defaulting to 0")
    
if len(sys.argv) > 2:
    png_palette = sys.argv[2]
    # Parse the comma-separated list into a Python list
    user_palette = [c.strip() for c in png_palette.split(",") if c.strip()]
else:
    user_palette = []

#############################################################################################################################

# Loading in metadata
metadata = pd.read_csv(os.path.join("Output", "selected_pdbs_metadata.csv"))

# Loading in fibril_info
fibril_info = (
    pd.read_csv(os.path.join("Output", "PDBs", "COM_and_fibril.csv"))
    .drop(columns = ["x_com", "y_com", "z_com"])
)
fibril_info["PDB"] = fibril_info["pdb_id"].str.split("_").str[0]

# Getting path to stable regions
stable_regions_path = os.path.join(
    "Output", "stable_regions", "stable_regions.csv"
)

# Getting RMSD cluster groups
rmsd_cluster_groups_path = os.path.join("Output", "RMSD", "data", "RMSD_cluster_groups.csv")

#############################################################################################################################

pymol.pymol_argv = ["pymol", "-qc"]
pymol.finish_launching()

unique_pdb_path = os.path.join("Output", "PDBs", "unique_chains")
unique_png_path = os.path.join("Output", "PNG", "unique_chains")
colour_residues = set_coloring(png_colouring, user_palette, 
                               metadata, fibril_info, 
                               stable_regions = stable_regions_path)
process_pdb_directory(unique_pdb_path, unique_png_path, 
                      colour_residues, apply_transparency=False, 
                      fibril_info=fibril_info, reference_pdb=0)
create_png_grid(
    pdb_path = unique_pdb_path,
    png_dir = unique_png_path,
    output_name = "combined_grid.png",
    transparent = False,
    font_size = 140
    )
create_grouped_png_grid(
    pdb_path = unique_pdb_path,
    png_dir = unique_png_path,
    cluster_path = rmsd_cluster_groups_path,
    output_name = "RMSD_cluster_groups.png",
    transparent = False,
    font_size = 140,
    title_size = 140,
    line_thickness = 10,
    top_padding = 50
)


asym_pdb_path = os.path.join("Output", "PDBs", "asymetric_unit")
asym_png_path = os.path.join("Output", "PNG", "asymetric_unit")
colour_residues = set_coloring(png_colouring, user_palette, 
                               metadata, fibril_info, 
                               stable_regions = stable_regions_path)
process_pdb_directory(asym_pdb_path, asym_png_path, 
                      colour_residues, apply_transparency=True, 
                      fibril_info=fibril_info, reference_pdb=0)

create_png_grid(
    pdb_path = asym_pdb_path,
    png_dir = asym_png_path,
    output_name = "combined_grid.png",
    transparent = True,
    font_size = 140
    )

create_grouped_png_grid(
    pdb_path= asym_pdb_path,
    png_dir = asym_png_path,
    cluster_path = rmsd_cluster_groups_path,
    output_name = "RMSD_cluster_groups.png",
    transparent = True,
    font_size = 140,
    title_size = 140,
    line_thickness = 10,
    top_padding = 50
)

cmd.quit()
