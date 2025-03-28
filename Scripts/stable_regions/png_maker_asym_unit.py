import pymol
from pymol import cmd
import os 
import pandas as pd
from PIL import Image

#############################################################################################################################

save_path = "Output/stable_regions/pngs_asym_unit"

# Making output directory
if not os.path.exists(save_path):
    os.mkdir(save_path)

#############################################################################################################################

#### TODO ####
# Scale image size by the number of protofilaments

# Reading in all the .pdb files in the unique chains directory
pdb_path = "Output/PDBs/asymetric_unit"

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]

# Reading in stable regions data
df = pd.read_csv("Output/stable_regions/stable_regions.csv")

# Reading in fibril_info
fibril_info = pd.read_csv("Output/PDBs/COM_and_fibril.csv")
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

from itertools import cycle

colours = [color for _, color in zip(stable_regions, cycle(colours_list))]

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

    # Calculate the required viewport size
    #image_width = monomer_size * num_monomers
    #image_height = image_width  # Keep it square for consistency

    # Set the PyMOL viewport size
    #cmd.viewport(image_width, image_height)

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

            #print("Matching chains:", matching_chains)
            #print("Non-matching chains (to be set transparent):", non_matching_chains)

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

        # Load the saved image
        #img = Image.open(png_path)

        # Auto-crop white space (finds non-white bounding box)
        #img = img.crop(img.getbbox())

        # Compute the new image size
        #new_width = int(img.width * num_monomers)
        #new_height = int(img.height * num_monomers)

        # Resize the image while maintaining aspect ratio
        #resized_img = img.resize((new_width, new_height), Image.LANCZOS)

        # Save the resized image
        #resized_img.save(png_path)

# Close PyMol
cmd.quit()
