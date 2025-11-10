import pymol
from pymol import cmd
import os 
import pandas as pd
from itertools import cycle
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

#############################################################################################################################

save_dir = os.path.join("Output", "PNG", "stable_regions")
# Making output directory
if not os.path.exists(save_dir):
    os.mkdir(save_dir)


save_path = os.path.join("Output", "PNG", "stable_regions", "unique_chains")

# Making output directory
if not os.path.exists(save_path):
    os.mkdir(save_path)

#############################################################################################################################

# Reading in all the .pdb files in the unique chains directory
pdb_path = os.path.join("Output", "PDBs", "unique_chains")

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]

# Reading in stable regions data
df = pd.read_csv(os.path.join("Output", "stable_regions", "stable_regions.csv"))

# Extract the first column
stable_regions = df.iloc[:, 0].tolist()

print(f"\n Stable Regions: {stable_regions}")

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

# Register colours in PyMOL
for idx, hex_colour in enumerate(region_colours):
    rgb = mcolors.to_rgb(hex_colour)
    cmd.set_color(f"region_colour_{idx}", rgb)


# Set PyMOL to run without GUI
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

# Load a reference structure (just take the first one, for example)
reference_path = os.path.join(pdb_path, pdb_files[0])

for i in pdb_files:

    # Joining path and filename
    current_path = os.path.join(pdb_path, i)
    
    # Setting name of loaded structures to be the pdb code
    # removing the .pdb extension
    object_name = i.split(".")[0]

    # Loading in the structure and setting the name 
    cmd.load(current_path, object_name)

    cmd.load(reference_path, "reference")

    # Align the structure to the reference using alpha carbons
    cmd.align(f"{object_name} and name CA", "reference and name CA")

    # Making the structure green
    #cmd.color("green", object_name)
    cmd.color("white", object_name)

    # Set the cartoon representation
    cmd.show("cartoon")

    # Set the cartoon loop radius
    cmd.set("cartoon_loop_radius", 1)

    # Removing dashed lines between disconnected residues
    cmd.set("cartoon_gap_cutoff", 0)

    # Iterate over each region and color it
    for j, row in df.iterrows():
        region_range = row["residues"] # Getting the residues in the current stable region
        cmd.color(f"region_colour_{j}", f"{object_name} and resi {region_range}")

    # Set the background color to transparent
    cmd.bg_color("white")  # Set to white to make sure transparency works
    cmd.set("ray_opaque_background", 0)

    # Orienting camera to focus on the reference so it is consistent
    cmd.orient(object_name)

    # Deleting the reference structure
    cmd.delete("reference")

    # Refresh makes sure everything is updated before saving
    cmd.refresh()

    # Save the image as a PNG with a transparent background
    png_path = os.path.join(save_path, f"{object_name}.png")
    cmd.png(png_path, dpi=300, ray=1, quiet=1)

    # Clear the current object for the next iteration
    cmd.delete(object_name)

############################################################################
############################################################################
############################################################################

save_path = os.path.join("Output", "PNG", "stable_regions", "asym_unit")

# Making output directory
if not os.path.exists(save_path):
    os.mkdir(save_path)


# Reading in all the .pdb files in the unique chains directory
pdb_path = os.path.join("Output", "PDBs", "asymetric_unit")

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]

# Reading in stable regions data
df = pd.read_csv(os.path.join("Output", "stable_regions", "stable_regions.csv"))

# Reading in fibril_info
fibril_info = pd.read_csv(os.path.join("Output", "PDBs", "COM_and_fibril.csv"))
current_fibril_info = fibril_info # Assiging for use in the loop

# Extract the first column
stable_regions = df.iloc[:, 0].tolist()

print(f"\n Stable Regions: {stable_regions}")

# Defining colours
colours_list = [
    "red", "blue", "gold", "br4", 
    "orange", "cyan", "green", "pink", "brown",
    "forest", "magenta", "deepteal", "br3", "salmon", 
    "grey"]

colours = [color for _, color in zip(stable_regions, cycle(colours_list))]

# Load a reference structure (just take the first one, for example)
reference_path = os.path.join(pdb_path, pdb_files[0])


for i in pdb_files:

    # Joining path and filename
    current_path = os.path.join(pdb_path, i)
    
    # Setting name of loaded structures to be the pdb code
    # removing the .pdb extension
    object_name = i.split(".")[0]

    # Keep everything before the first
    object_name = object_name.split("_")[0]

    # Getting the fibril information for the current structure
    current_fibril_info = fibril_info[fibril_info['PDB'].str.contains(object_name, na=False)]

    # Identifying if intra-PDB is present
    unique_pdb_ids = current_fibril_info['pdb_id'].unique().tolist()

    # Get the number of monomers for image size scaling
    num_monomers = current_fibril_info["fibril"].max()

    # Define a fixed monomer size (in pixels per monomer)
    monomer_size = 150  # Adjust this as needed

    # Set the PyMOL viewport size
    if len(unique_pdb_ids) > 1:
        for j in unique_pdb_ids:
            
            object_name = j

            # Loading in the structure and setting the name
            cmd.load(current_path, object_name)

            # Get list of chains in the loaded structure
            chains = cmd.get_chains(object_name)
            #print("Loaded chains:", chains)

            # Filter fibril_info to get the chains for the current pdb_id
            expected_chains = current_fibril_info[current_fibril_info['pdb_id'] == j]['chain'].unique().tolist()
            #print("Expected chains:", expected_chains)

            # Identify matching and non-matching chains
            matching_chains = [chain for chain in chains if chain in expected_chains]
            non_matching_chains = [chain for chain in chains if chain not in expected_chains]

            # Set non-matching chains to 50% transparency
            for chain in non_matching_chains:
                cmd.set("cartoon_transparency", 0.75, f"{object_name} and chain {chain}")

            # Loading in the reference structure and setting the name
            cmd.load(reference_path, "reference")

            # Align the structure to the reference using alpha carbons
            cmd.align(f"{object_name} and name CA", "reference and name CA")

            # Making the structure green
            cmd.color("white", object_name)

            # Set the cartoon representation
            cmd.show("cartoon")

            # Set the cartoon loop radius
            cmd.set("cartoon_loop_radius", 1)

            # Removing dashed lines between disconnected residues
            cmd.set("cartoon_gap_cutoff", 0)

            # Iterate over each region and color it
            for j, row in df.iterrows():
                region_range = row["residues"] # Getting the residues in the current stable region
                cmd.color(f"region_colour_{j}", f"{object_name} and resi {region_range}")

            # Set the background color to transparent
            cmd.bg_color("white")  # Set to white to make sure transparency works
            cmd.set("ray_opaque_background", 0)

            # Orienting camera to focus on the reference so it is consistent
            cmd.orient(object_name)

            # Deleting the reference structure
            cmd.delete("reference")

            # Refresh makes sure everything is updated before saving
            cmd.refresh()

            # Save the image as a PNG with a transparent background
            png_path = os.path.join(save_path, f"{object_name}.png")
            cmd.png(png_path, dpi=300, ray=1, quiet=1)

            # Clear the current object for the next iteration
            cmd.delete(object_name)
    
    else:
        # Loading in the structure and setting the name 
        cmd.load(current_path, object_name)

        # Loading in the reference structure and setting the name
        cmd.load(reference_path, "reference")

        # Align the structure to the reference using alpha carbons
        cmd.align(f"{object_name} and name CA", "reference and name CA")

        # Making the structure green
        cmd.color("white", object_name)

        # Set the cartoon representation
        cmd.show("cartoon")

        # Set the cartoon loop radius
        cmd.set("cartoon_loop_radius", 1)

        # Removing dashed lines between disconnected residues
        cmd.set("cartoon_gap_cutoff", 0)

        # Iterate over each region and color it
        for j, row in df.iterrows():
            region_range = row["residues"] # Getting the residues in the current stable region
            cmd.color(f"region_colour_{j}", f"{object_name} and resi {region_range}")

        # Set the background color to transparent
        cmd.bg_color("white")  # Set to white to make sure transparency works
        cmd.set("ray_opaque_background", 0)

        # Orienting camera to focus on the reference so it is consistent
        cmd.orient(object_name)

        # Deleting the reference structure
        cmd.delete("reference")

        # Refresh makes sure everything is updated before saving
        cmd.refresh()

        # Save the image as a PNG with a transparent background
        png_path = os.path.join(save_path, f"{object_name}.png")
        cmd.png(png_path, dpi=300, ray=1, quiet=1)

        # Clear the current object for the next iteration
        cmd.delete(object_name)

# Close PyMol
cmd.quit()




