import pymol
# Import PyMOL modules
from pymol import cmd
import os
import pandas as pd
import numpy as np

#############################################################################################################################

save_path = os.path.join("Output", "PDBs", "fibrils_extended")

# Making output directory
# Check if the directory already exists
if not os.path.exists(save_path):
    os.mkdir(save_path)

#############################################################################################################################

# importing EMD_data
file_path = os.path.join("Output", "emdb_data.txt")

# Read the text file into a DataFrame
df = pd.read_csv(file_path, sep="\t") 

# Setting the path to the PDB asymetric units
pdb_path = os.path.join("Output", "PDBs", "asymetric_unit")

# Initialize a DataFrame to store the extracted PDB, chain, and B-factor values
combined_chain_fibril = pd.DataFrame()
combined_exterior_chains = pd.DataFrame()

# Set PyMOL to run without GUI
pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

for index, row in df.iterrows():

    i = index

    # Setting PDB to import
    pdb = df["pdb"].iloc[i]
    print("PDB: ", pdb)
    twist = float(df["twist"].iloc[i])
    rise = float(df["rise"].iloc[i])
    
    # Applying transformations to a complete layer so converting twist and rise to per 4.8A
    ratio = rise / 4.8
    
    if ratio >= 0.2 and ratio < 0.3:
        rise = rise * 4
        twist = twist * 4
    
    if ratio >= 0.3 and ratio < 0.4:
        rise = rise * 3
        twist = twist * 3
        
    if ratio >= 0.4 and ratio < 0.6:
        rise = rise * 2
        twist = twist * 2    
    
    if ratio >= 1.5 and ratio <= 2.2:
        rise = rise / 2
        twist = twist / 2

    if ratio >=2.8 and ratio <= 3.2:
        rise = rise / 3
        twist = twist / 3
    
    # Converting twist to    
    if twist > 180:
        twist = twist - 180
    if twist < -180:
        twist = twist + 180
        
    print("Twist: ", twist)    
    print("Rise: ", rise, "\n")

    if np.isnan(twist) or np.isnan(rise):
        pass
    else:

        # Refreshing pymol session
        cmd.reinitialize()
        
        # Creating path to PDB asymetric unit
        current_pdb_path = os.path.join(pdb_path, pdb + "_asym_unit.pdb")

        # Load the PDB file
        cmd.load(current_pdb_path, pdb)

        print("Loaded PDB file: ", pdb, "\n")

        # Duplicating the asymetric unit 9 times to create a layer depth of 10 chains

        # Duplicating the monomer
        cmd.create("tmp", pdb)

        # Creating new layers
        for i in range(0, 9):

            # Moving the monomer down by 4.8A
            cmd.translate([0, 0, -rise], "tmp")
            
            # Roate around the y-axis
            cmd.rotate([0, 0, 1], -twist, "tmp")

            # Adding new chain to the original structure
            cmd.copy_to(pdb, "tmp")

        # Deleting  the temporary dublicated object
        cmd.delete("tmp")

        # Saving the extended fibril structure
        pymol.cmd.save(os.path.join(save_path, pdb + ".pdb"))

        # Suppressing the output from PyMOL's iterate command
        with open(os.devnull, 'w') as f:
            old_stdout = os.dup(1)
            os.dup2(f.fileno(), 1)

            # Extract unique chains and B-factor
            chain_bfactor = {}
            cmd.iterate(f"{pdb} and name CA", "chain_bfactor.setdefault(chain, b)", space={"chain_bfactor": chain_bfactor})

            os.dup2(old_stdout, 1)
            os.close(old_stdout)

        # Creating a pandas dataframe with the extracted data
        pdb_info = pd.DataFrame([{"PDB": pdb, "chain": chain, "fibril": bfactor} for chain, bfactor in chain_bfactor.items()])

        # Append the current pdb info to the combined DataFrame
        combined_chain_fibril = pd.concat([combined_chain_fibril, pdb_info], ignore_index=True)

        ########################################################
        ### Identifying head and tail chains for each fibril ###
        ########################################################

        # Calculate the center of mass for each chain
        chain_com = []
        for chain in chain_bfactor.keys():
            selection = f"{pdb} and chain {chain}"
            com = cmd.centerofmass(selection)
            chain_com.append({"PDB": pdb, "chain": chain, "fibril": chain_bfactor[chain], "com_z": com[2]})

        chain_com_df = pd.DataFrame(chain_com)

        # Determine chains with the highest and lowest center of mass z-coordinate for each fibril
        max_com = chain_com_df.loc[chain_com_df.groupby("fibril")["com_z"].idxmax()]
        min_com = chain_com_df.loc[chain_com_df.groupby("fibril")["com_z"].idxmin()]

        # Add a new column to indicate if it's max or min
        max_com['location'] = 'head'
        min_com['location'] = 'tail'

        # Combine max and min into a single DataFrame
        combined_exterior_chains = pd.concat([combined_exterior_chains, max_com, min_com], ignore_index=True)
        
# Drop the com_z column
combined_exterior_chains = combined_exterior_chains.drop(columns=['com_z'])

# Save the combined results
combined_chain_fibril.to_csv(os.path.join("Output", "PDBs", "fibrils_extended", "chain_fibril.csv"), index=False)
combined_exterior_chains.to_csv(os.path.join("Output", "PDBs", "fibrils_extended", "exterior_chains.csv"), index=False)
        
# Close PyMol
cmd.quit()

