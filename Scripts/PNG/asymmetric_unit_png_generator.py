import pymol
from pymol import cmd
import os 
import pandas as pd
import sys

#############################################################################################################################
png_path = os.path.join("Output", "PNG")
save_path = os.path.join("Output", "PNG", "asymetric_unit")

if not os.path.exists(png_path):
    os.mkdir(png_path)
if not os.path.exists(save_path):
    os.mkdir(save_path)

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

if png_colouring == 0: 
    selected_colour = sys.argv[2]
    print(f"Selected colour = {selected_colour}")
    
if len(sys.argv) > 3:
    png_palette = sys.argv[3]
    # Parse the comma-separated list into a Python list
    user_palette = [c.strip() for c in png_palette.split(",") if c.strip()]
else:
    user_palette = []

#############################################################################################################################

# Reading in all the .pdb files in the asymmetric units directory
pdb_path = os.path.join("Output", "PDBs", "asymetric_unit")

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]

# Loading in metadata
metadata = pd.read_csv(os.path.join("Output", "selected_pdbs_metadata.csv"))

# Loading in fibril_info
fibril_info = (
    pd.read_csv(os.path.join("Output", "PDBs", "COM_and_fibril.csv"))
    .drop(columns = ["x_com", "y_com", "z_com"])
)
fibril_info["PDB"] = fibril_info["pdb_id"].str.split("_").str[0]

#############################################################################################################################

def color_white():
    cmd.color("white")

def validate_colour(colour):
    valid_colours = {name for name, _ in cmd.get_color_indices()}
    if colour not in valid_colours:
        print(f"Warning: '{colour}' is not a recognised PyMOL colour. Reverting to 'white'.")
        return "white"
    return colour

if png_colouring == 0:
    # Checking if colour provided is valid
    selected_colour_valid = validate_colour(selected_colour)

    def color_residues():
        cmd.color(selected_colour_valid)
else:
    def color_residues():
        cmd.color("white")
        
if png_colouring == 1:              

    protein_column = metadata["Protein"].str.lower()

    asyn = protein_column.str.contains("-synuclein").any()
    tau = protein_column.str.contains("tau").any()
    abeta = protein_column.str.contains("amyloid-").any()
    
    
    # Asyn Colouring
    if asyn: 
        # Function to set color for specific residue ranges
        def color_residues():
            # N-terminus
            cmd.color("tv_blue", "resi 1-60")
            # NAC
            cmd.color("tv_red", "resi 61-95")
            # C-terminus
            cmd.color("tv_green", "resi 96-140")

    # Tau Colouring
    elif tau: 
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
            
    # A-beta Colouring
    elif abeta: 
        # Function to set color for specific residue ranges
        def color_residues():
            # N-term
            cmd.color("tv_blue", "resi 1-16")
            # Central Hydrophobic Core
            cmd.color("black", "resi 17-21")
            # Turn Region
            cmd.color("tv_red", "resi 24-27")
            # Second Hydrophobic Region
            cmd.color("tv_orange", "resi 28-35")
            # C-term
            cmd.color("tv_green", "resi 36-42")
            
    else:
        print("Domain colouring is not supported for this protein, feel free to add in your own using the examples in Scripts/PNG/asymmetric_unit_png_generator.py")

if png_colouring in [2, 3]:

    # Pick which metadata column drives the colouring
    target_col = "fibril" if png_colouring == 2 else "polymorph"

    # Pre-generate a colour palette so each group gets its own colour
    fibril_groups = fibril_info[target_col].unique()
    colour_map = {}

    # default colours
    default_colours = [
        "tv_red", "tv_blue", "tv_green", "tv_orange", "yellow",
        "purple", "cyan", "salmon", "lime", "deepblue"
    ]

    # ensure user_palette exists
    user_palette = user_palette or []

    # validate the user colours BEFORE trying to use them
    validated_palette = [validate_colour(c) for c in user_palette]

    # only keep default colours that aren't already in the validated palette
    filtered_defaults = [c for c in default_colours if c not in validated_palette]

    # combine user colours + filtered defaults
    pymol_colours = validated_palette + filtered_defaults

    for idx, f in enumerate(fibril_groups):
        colour_map[f] = pymol_colours[idx % len(pymol_colours)]

    def color_residues():

        curr = fibril_info[fibril_info["PDB"] == pdb_name]

        if curr.empty:
            print(f"No fibril_info entry found for {pdb_name}, nothing coloured.")
            return

        for _, row in curr.iterrows():
            chain = row["chain"]
            group_val = row[target_col]
            colour = colour_map[group_val]

            cmd.color(colour, f"{object_name} and chain {chain}")


# Set PyMOL to run without GUI
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

for i in pdb_files:

    # Joining path and filename
    current_path = os.path.join(pdb_path, i)

    # Setting name of loaded structures to be the pdb code
    # removing the .pdb extension
    object_name = i.split(".")[0]

    # Removing _asym_unit
    pdb_name = object_name.split("_")[0]

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

