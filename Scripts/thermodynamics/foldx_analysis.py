# Working with PDB structures using PyFoldX

# Importing Packages
from pyfoldx.structure import Structure
from pyfoldx.foldx import foldxHandler
import sys
import pandas as pd
import os
import re
import shutil

#############################################################################################################################

# Setting FoldX Parameters

# Read command line arguments
ph = sys.argv[1] if len(sys.argv) > 1 else None
temp = sys.argv[2] if len(sys.argv) > 2 else None
ionstrength = sys.argv[3] if len(sys.argv) > 3 else None

# Update parameters only if valid numeric input is provided
def try_float(value, default):
    try:
        return float(value)
    except (ValueError, TypeError):
        return default

parameters = {
    "pH": try_float(ph, 7.4),
    "temperature": try_float(temp, 298),
    "ionStrength": try_float(ionstrength, 0.150)
}

print(parameters)

#############################################################################################################################

    
# Path to unique PDB chains
file_path = os.path.join("Output", "PDBs", "fibrils_extended")

# Path to save repaired extended fibrils
save_path = os.path.join("Output", "PDBs", "fibrils_extended", "repaired")

# Making output directory
# Check if the directory already exists
if not os.path.exists(save_path):
    os.mkdir(save_path)

thermodynamic_path = os.path.join("Output", "thermodynamics")

# Check if the directory already exists
if not os.path.exists(thermodynamic_path):
    os.mkdir(thermodynamic_path)

#############################################################################################################################

# Get all the .pdb files in the directory
pdb_files = [os.path.join(file_path, file) for file in os.listdir(file_path) if file.endswith(".pdb")]

# Initialize data frame to store the deltaG per residue 
combined_residue_energy = pd.DataFrame()

for pdb in pdb_files:

    # Extract the pdb_id
    pdb_id = os.path.splitext(os.path.basename(pdb))[0]

    # Split the string after the last "_"
    pdb_code = re.split(r'_(?!.*_)', pdb_id)[0]

    # Creating path for repaired structure to be saved to
    save_pdb_path = os.path.join(save_path, pdb_id + "_repaired.pdb")

    print("\n\nAnalysing: " + pdb_id + "\n\n")

    # Initializing seqRepaired
    seqRepaired = pd.DataFrame()

    # Repeating getResiduesEnergy() as for some reason it can error on different structures for each run
    # Check if seqRepaired is empty and repeat if necessary
    retry_count = 0
    max_retries = 5
    while seqRepaired.empty and retry_count < max_retries:
        retry_count += 1

        if not os.path.exists(save_pdb_path):

            # Load the structure using the Structure class
            st = Structure(code=pdb_id, path=pdb)
            
            # We create a new Structure object with the result of the repair() method
            stRepaired = st.repair()

            # Saving the repaired structure
            stRepaired.toPdbFile(save_pdb_path)
        
        else:
            
            # If repaired structure has already been created, load it directly
            stRepaired = Structure(code=pdb_id, path = save_pdb_path)

        # Extract all energy values
        seqRepaired = foldxHandler.getResiduesEnergy(stRepaired, consider_waters = False, other_parameters=parameters)

    print("\nEnergy values calculated for: " + pdb_id + "\n")
    
    # Reset the index to convert existing index columns into regular columns
    seqRepaired.reset_index(inplace=True)
    
    # Adding pdb_id to the start of the data frame
    seqRepaired.insert(0, "PDB", pdb_id)

    # Print the first 5 rows of seqRepaired
    print("\n", seqRepaired.head())

    # Appending current pdb results to combined results
    combined_residue_energy = pd.concat([combined_residue_energy, seqRepaired], ignore_index=True)

# Saving the combined results
combined_residue_energy.to_csv(os.path.join("Output", "thermodynamics", "foldx_stability.csv"), index=False)

# Deleting temporary .foldx folder
for folder in os.listdir("."):
    if folder.startswith(".foldx") and os.path.isdir(folder):
        try:
            shutil.rmtree(folder)
        except OSError:
            pass