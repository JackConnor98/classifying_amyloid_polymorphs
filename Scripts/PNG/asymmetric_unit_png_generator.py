import pymol
from pymol import cmd
import os 

#############################################################################################################################
png_path = os.path.join("Output", "PNG")
save_path = os.path.join("Output", "PNG", "asymetric_unit")

if not os.path.exists(png_path):
    os.mkdir(png_path)
if not os.path.exists(save_path):
    os.mkdir(save_path)

#############################################################################################################################

# Reading in all the .pdb files in the asymmetric units directory
pdb_path = os.path.join("Output", "PDBs", "asymetric_unit")

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]

def color_white():
    cmd.color("white")

def color_residues():
    cmd.color("white")

# Asyn Colouring
if any("6cu7" in file for file in pdb_files):
    # Function to set color for specific residue ranges
    def color_residues():
        # N-terminus
        cmd.color("tv_blue", "resi 1-60")
        # NAC
        cmd.color("tv_red", "resi 61-95")
        # C-terminus
        cmd.color("tv_green", "resi 96-140")

# Tau Colouring
if any("5o3l" in file for file in pdb_files):
    # Function to set color for specific residue ranges
    def color_residues():
        # R1
        cmd.color("tv_red", "resi 243-273")
        # R2
        cmd.color("tv_green", "resi 274-304")
        # R3
        cmd.color("tv_orange", "resi 305-335")
        # R4
        cmd.color("tv_blue", "resi 336-367")

# Set PyMOL to run without GUI
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

for i in pdb_files:

    # Joining path and filename
    current_path = os.path.join(pdb_path, i)

    # Setting name of loaded structures to be the pdb code
    # removing the .pdb extension
    object_name = i.split(".")[0]

    # Loading in the structure and setting the name 
    cmd.load(current_path, object_name)

    # Set the cartoon representation
    cmd.show("cartoon")

    # Set the cartoon loop radius
    cmd.set("cartoon_loop_radius", 1)

    # Set all residues to white
    color_white()

    # Color specific residues
    color_residues()

    # Set the background color to transparent
    cmd.bg_color("white")  # Set to white to make sure transparency works
    cmd.set("ray_opaque_background", 0)

    # Refresh makes sure everything is updated before saving
    cmd.refresh()

    # Save the image as a PNG with a transparent background
    png_path = os.path.join(save_path, f"{object_name}.png")
    cmd.png(png_path, dpi=300, ray=1, quiet=1)

    # Clear the current object for the next iteration
    cmd.delete(object_name)
    
# Close PyMol
cmd.quit()

