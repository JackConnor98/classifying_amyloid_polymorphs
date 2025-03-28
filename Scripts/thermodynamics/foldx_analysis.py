# Working with PDB structures using PyFoldX

# Importing Packages
from pyfoldx.structure import Structure
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import re

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
    pdb_id = pdb.split('.')[0]
    pdb_id = pdb_id.split("/")[-1]

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
        seqRepaired = stRepaired.getResiduesEnergy()

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
combined_residue_energy.to_csv('Output/thermodynamics/foldx_stability.csv', index=False)
