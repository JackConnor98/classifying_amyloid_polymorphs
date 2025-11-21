# Working with PDB structures using PyFoldX

# Importing Packages
from pyfoldx.structure import Structure
from pyfoldx.foldx import foldxHandler
import sys
import pandas as pd
import os
import zipfile
import re
import subprocess
from concurrent.futures import ProcessPoolExecutor
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

# Setting the number of processes
if len(sys.argv) > 4: 
    try:
        processes = int(sys.argv[4]) 
        print(f"User has specified {processes} processes")
    except ValueError: 
        print("Warning: invalid argument — defaulting to 4 processes") 
        processes = 4 
else: 
    processes = 4 
    print("No input given, defaulting to 4 processes")
    
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

# Get all the .pdb files in the directory
pdb_files = [os.path.join(file_path, file) for file in os.listdir(file_path) if file.endswith(".pdb")]

#############################################################################################################################

##########################################
### Extracting pre-repaired structures ###
##########################################

def extract_matching_pdbs(zip_path, extract_to, pdb_files):
    # Build a set of expected base names (no extension)
    expected = {os.path.splitext(os.path.basename(p))[0] for p in pdb_files}

    # Ensure destination folder exists
    os.makedirs(extract_to, exist_ok = True)

    with zipfile.ZipFile(zip_path, "r") as zip_ref:
        for member in zip_ref.namelist():
            # Only consider .pdb files
            if not member.lower().endswith(".pdb"):
                continue

            # Extract the name without extension
            base = os.path.splitext(os.path.basename(member))[0]

            # Extract only if it matches one of the expected names
            if base in expected:
                zip_ref.extract(member, extract_to)

    print(f"Extraction complete. Files saved to: {extract_to}")

extract_matching_pdbs("Repaired_pdbs.zip", save_path, pdb_files)

#########################################
### Repairing PDBs for foldx analysis ###
#########################################

# Getting the path to the foldx.exe
foldx_path = os.environ.get("FOLDX_BINARY")

def repair_single(pdb):
    
    pdb_base = os.path.splitext(os.path.basename(pdb))[0]
    final_file = os.path.join(save_path, f"{pdb_base}_repaired.pdb")
    
    # Skip if already repaired
    if os.path.exists(final_file):
        print(f"✅  Skipping {pdb_base}, repaired file already exists.")
        return pdb
    
    retry_count = 0
    max_retries = 5

    while retry_count < max_retries:
        print(f"Attempt {retry_count}/{max_retries} to repair and extract energies for {pdb}")
        cmd = [
            foldx_path,
            "--command=RepairPDB",
            f"--pdb-dir={file_path}",
            f"--output-dir={save_path}",
            f"--pdb={os.path.basename(pdb)}"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            print(f"⚠️  Error repairing {pdb}")
            retry_count += 1
        else:
            print(f"✅  Repaired {pdb}")
            retry_count = max_retries  # Exit loop
        
    # Rename repaired file
    repaired_file = os.path.join(save_path, f"{pdb_base}_Repair.pdb")
    if os.path.exists(repaired_file):
        new_name = os.path.join(save_path, f"{pdb_base}_repaired.pdb")
        os.rename(repaired_file, new_name)
        print(f"✅  Repaired and renamed: {new_name}")
    else:
        print(f"⚠️  Could not find repaired file for {pdb_base}")

    # Manually delete FoldX clutter
    for filename in os.listdir(save_path):
        if (
            filename.startswith("FXout")
            or filename.endswith(".fxout")
        ):
            try:
                os.remove(os.path.join(save_path, filename))
            except OSError:
                pass

    # Remove .foldx folder manually
    if os.path.isdir(".foldx"):
        for root, dirs, files in os.walk(".foldx", topdown=False):
            for f in files:
                try:
                    os.remove(os.path.join(root, f))
                except OSError:
                    pass
            for d in dirs:
                try:
                    os.rmdir(os.path.join(root, d))
                except OSError:
                    pass
        try:
            os.rmdir(".foldx")
        except OSError:
            pass    
        
    return pdb

with ProcessPoolExecutor(max_workers=processes) as executor:
    executor.map(repair_single, pdb_files)

#############################################################################################################################

######################################
### Calculating deltaG per residue ###
######################################

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
        
# Deleting temporary rotabase.txt
if os.path.exists("rotabase.txt"):
    try:
        os.remove("rotabase.txt")
    except FileNotFoundError:
        pass

# Deleting temporary Unrecognized_molecules.txt
if os.path.exists("Unrecognized_molecules.txt"):
    try:
        os.remove("Unrecognized_molecules.txt")
    except FileNotFoundError:
        pass
