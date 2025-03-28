import pymol
from pymol import cmd
import os 
import sys

#############################################################################################################################

save_path = "Output/PDBs/Alignments"

# Making output directory
# Check if the directory already exists
if not os.path.exists(save_path):
    os.mkdir(save_path)

#############################################################################################################################

# Reading in all the .pdb files in the unique chains directory
pdb_path = "Output/PDBs/unique_chains"

# Get all the .pdb files in the directory
pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]

# Redirect stdout to a file
output_file = os.path.join(save_path, "rmsd_output.txt")
stdout_orig = sys.stdout
with open(output_file, 'w') as f:
    sys.stdout = f

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

    # Aligning all single chains from the same pdb file
    pymol.cmd.extra_fit(selection='(all)', reference=None, method='align')

    # Saving Alignment
    pymol.cmd.save(os.path.join(save_path, "all_unique_chains_aligned.cif"))
    pymol.cmd.save(os.path.join(save_path, "all_unique_chains_aligned.pse"))

    # Close PyMol
    cmd.quit()

# Restore stdout to its original state
sys.stdout = stdout_orig