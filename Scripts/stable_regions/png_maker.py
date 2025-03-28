import pymol
from pymol import cmd
import os 
import pandas as pd

#############################################################################################################################

save_path = "Output/stable_regions/pngs"

# Making output directory
if not os.path.exists(save_path):
    os.mkdir(save_path)

#############################################################################################################################

# Reading in all the .pdb files in the unique chains directory
pdb_path = "Output/PDBs/unique_chains"
#pdb_path = "Output/PDBs/asymetric_unit"

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]

# Reading in stable regions data
df = pd.read_csv("Output/stable_regions/stable_regions.csv")

# Extract the first column
stable_regions = df.iloc[:, 0].tolist()

print(f"\n Stable Regions: {stable_regions}")

# Defining colours
colours_list = [
    "red", "blue", "gold", "br4", 
    "orange", "cyan", "green", "pink", "brown",
    "forest", "magenta", "deepteal", "br3", "salmon", 
    "grey"]

from itertools import cycle

colours = [color for _, color in zip(stable_regions, cycle(colours_list))]

# Set PyMOL to run without GUI
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

# Load a reference structure (just take the first one, for example)
reference_path = os.path.join(pdb_path, pdb_files[5])


for i in pdb_files:

    # Joining path and filename
    current_path = os.path.join(pdb_path, i)
    
    # Setting name of loaded structures to be the pdb code
    # removing the .pdb extension
    object_name = i.split(".")[0]
    #object_name = object_name.split("_")[0]

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

    # Iterate over each region and color it
    for j, region in enumerate(stable_regions):
        colour = colours[j % len(colours)]  # Cycle through colours
        cmd.color(colour, f"{object_name} and resi {region}")  # Apply the color

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
