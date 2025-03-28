import pymol
from pymol import cmd
import os 
import pandas as pd

#############################################################################################################################

save_path = "Output/PDBs/coloured_by_stable_regions"

# Making output directory
# Check if the directory already exists
if not os.path.exists(save_path):
    os.mkdir(save_path)

#############################################################################################################################

# Reading in all the .pdb files in the unique chains directory
pdb_path = "Output/PDBs/unique_chains"

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]

#############################################################################################################################

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

for i in pdb_files:

    # Refreshing the pymol session
    cmd.reinitialize()

    # Joining path and filename
    current_path = os.path.join(pdb_path, i)

    # Setting name of loaded structures to be the PDB code
    object_name = i.split(".")[0].split("_")[0]  # Remove .pdb and _asym_unit

    # Printing progress information
    print("\nLoading: " + object_name)

    # Loading in the structure and setting the name 
    cmd.load(current_path, object_name)

    # Making the structure green
    #cmd.color("green", object_name)
    cmd.color("white", object_name)

    cmd.hide("cartoon")
    cmd.show("sticks")

    # Iterate over each region and color it
    for j, region in enumerate(stable_regions):
        colour = colours[j % len(colours)]  # Cycle through colours
        cmd.color(colour, f"{object_name} and resi {region}")  # Apply the color

    # Saving Coloured Structure
    cmd.save(os.path.join(save_path, object_name + ".pse"))

    # Delete everything to clean the session for the next PDB
    cmd.delete('all')




# Close PyMol
cmd.quit()
