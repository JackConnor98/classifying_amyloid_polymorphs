import re
import pandas as pd

# Open the CIF file and read all lines
with open('Output/PDBs/Alignments/all_unique_chains_aligned.cif', 'r') as f:
    lines = f.readlines()

# Initialize lists to store data
amyloid_names = []
atom_lines = []

current_name = ""

# Process each line
for line in lines:
    if line.startswith("data"):
        current_name = re.search(r'data_(\w+)', line).group(1)
        amyloid_names.append(current_name)
    elif line.startswith("ATOM"):
        if " CA " in line:
            line = line.strip()
            # Extract necessary fields and add current_name
            fields = line.split()
            fields.append(current_name)
            atom_lines.append(fields)

# Convert to DataFrame
df = pd.DataFrame(atom_lines, columns=['ATOM', 'atom_number', 'type_symbol', 'atom_id', 'label_alt_id',
                                       'residue_name', 'label_asym_id', 'label_entity_id',
                                       'residue_number', 'pdbx_PDB_ins_code', 'x', 'y',
                                       'z', 'occupancy', 'fibril_number', 'charge',
                                       'chain', 'pdbx_PDB_model_num', 'pdb_id'])

# Create new column 'pdb' by splitting 'pdb_id'
df['pdb'] = df['pdb_id'].str.split('_').str[0]

# Selecting columns of interest
df = df[['pdb_id', 'pdb', 'chain', 'atom_id', 'residue_name', 'residue_number', 'x', 'y', 'z', 'fibril_number']]

# Save DataFrame to CSV file
df.to_csv('Output/PDBs/Alignments/all_unique_chains_aligned.csv', index=False)