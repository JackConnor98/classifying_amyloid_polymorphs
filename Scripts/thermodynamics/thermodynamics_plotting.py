import sys
import os
import re
import pandas as pd
import numpy as np
import plotnine

# Read command line arguments
remove_poorly_resolved = int(sys.argv[1])

################################################################################

# Creating file paths and directories
thermo_path = os.path.join("Output", "thermodynamics")

deltaG_analysis_path = os.path.join(thermo_path, "deltaG_analysis")
if not os.path.exists(deltaG_analysis_path):
    os.makedirs(deltaG_analysis_path)

single_pdb_path = os.path.join(thermo_path, "single_PDB_plots")
if not os.path.exists(single_pdb_path):
    os.makedirs(single_pdb_path)

residue_analysis_path = os.path.join(thermo_path, "residue_analysis")
if not os.path.exists(residue_analysis_path):
    os.makedirs(residue_analysis_path)

pearson_path = os.path.join(thermo_path, "pearson_correlation")
if not os.path.exists(pearson_path):
    os.makedirs(pearson_path)

deltaG_comp_path = os.path.join(thermo_path, "deltaG_comparisons")
if not os.path.exists(deltaG_comp_path):
    os.makedirs(deltaG_comp_path)

################################################################################

### Importing data ###
df = pd.read_csv(os.path.join(thermo_path, "foldx_stability.csv"))
metadata = pd.read_csv(os.path.join("Output", "selected_pdbs_metadata.txt"), sep="\t")
chain_fibril = pd.read_csv(os.path.join("Output", "PDBs", "fibrils_extended", "chain_fibril.csv"))

####################################
### Tidying foldx stability data ###
####################################

### Removing non-cryoEM structures ###
# Extracting PDB IDs with Method "cryoEM" from metadata
cryoem_ids = metadata[metadata["Method"] == "cryoEM"]["PDB ID"].tolist()
# Removing non-cryoEM structures from df
df = df[df["PDB"].isin(cryoem_ids)].reset_index(drop=True)



### Removing residues outside the fibril core ###



# Select and rename columns
metadata = metadata[['PDB ID', 'Residues Ordered']].rename(columns={'PDB ID': 'PDB', 'Residues Ordered': 'ordered_residues'})

# Add new columns for core start and end
metadata['core_start'] = None
metadata['core_end'] = None

# Check which PDBs are in metadata but not in df
missing_pdbs = set(metadata['PDB']) - set(df['PDB'])
print(missing_pdbs)

# Parse residue ranges
def parse_residue_range(residue_str):
    # Split on '-' or ','
    parts = re.split(r'[-,]', residue_str)
    parts = [p.strip() for p in parts if p.strip()]  # clean up empty strings
    return float(parts[0]), float(parts[-1])

metadata[['core_start', 'core_end']] = metadata['ordered_residues'].apply(
    lambda s: pd.Series(parse_residue_range(s))
)


print(metadata)