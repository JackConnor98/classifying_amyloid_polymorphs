import os
import pandas as pd
import subprocess
from Bio import PDB
from Bio.SeqUtils import seq1
from Bio import AlignIO
import json

root_folder = os.path.join("Output", "PDBs")
unique_chains_dir = os.path.join(root_folder, "unique_chains")
asym_dir = os.path.join(root_folder, "asymetric_unit")
output_dir = os.path.join(root_folder, "alignment")

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

# =========================================
# 1. Parse PDBs and build sequence dataframe
# =========================================
unique_chains = [file for file in os.listdir(unique_chains_dir) if file.endswith(".pdb")]

data = {'PDB': [], 'pos': [], 'residue': []}

parser = PDB.PDBParser(QUIET=True)

for pdb_file in unique_chains:
    path = os.path.join(unique_chains_dir, pdb_file)
    pdb_name = pdb_file.split(".")[0]

    try:
        structure = parser.get_structure('protein', path)
        if not any(structure):
            print(f"No data found in structure for {pdb_name}")
            continue
    except Exception as e:
        print(f"Error parsing file {pdb_name}: {e}")
        continue

    for model in structure:
        for chain in model:
            for residue in chain:
                if PDB.is_aa(residue):
                    data['PDB'].append(pdb_name)
                    data['pos'].append(residue.id[1])
                    data['residue'].append(residue.get_resname())

df = pd.DataFrame(data)
outfile = os.path.join(output_dir, 'sequence_alignment.csv')
df.to_csv(outfile, index=False)

# =========================================
# 2. Build sequence strings for each PDB
# =========================================
seqs = {}
for pdb, group in df.groupby("PDB"):
    seq = ''.join(seq1(r) for r in group.sort_values("pos")["residue"])
    seqs[pdb] = seq

# Write sequences to FASTA
fasta_file = os.path.join(output_dir, "sequences.fasta")
with open(fasta_file, "w") as f:
    for name, seq in seqs.items():
        f.write(f">{name}\n{seq}\n")

# =========================================
# 3. Run MUSCLE alignment
# =========================================
aligned_file = os.path.join(output_dir, "aligned_muscle.fasta")
cmd = ["muscle", "-align", fasta_file, "-output", aligned_file]

try:
    subprocess.run(cmd, check=True)
    print(f"Alignment complete! Aligned FASTA saved to {aligned_file}\n")
except subprocess.CalledProcessError as e:
    print("Error running MUSCLE:\n", e)

# =========================================
# 4. Read alignment
# =========================================
alignment = AlignIO.read(aligned_file, "fasta")

# =========================================
# 5. Update PDB residue numbers based on alignment
# =========================================

# # Dictionaries to store data
# sequence_data = {}

# io = PDB.PDBIO()

# for record in alignment:
#     pdb_name = record.id
#     aligned_seq = str(record.seq)
#     pdb_file = os.path.join(unique_chains_dir, f"{pdb_name}.pdb")
#     structure = parser.get_structure(pdb_name, pdb_file)
    
#     # Build new residue numbering taking gaps into account
#     residue_numbers = []
#     current_number = 1  # adjust to starting residue number if needed
#     for aa in aligned_seq:
#         if aa != "-":
#             residue_numbers.append(current_number)
#         current_number += 1  # increment for gaps too
    
#     # Extracting PDB from pdb_id
#     current_pdb = pdb_name.split('_', 1)[0]
    
#     # # Storing sequence alignment data in a dict
#     # current_entry = {"PDB": current_pdb,
#     #                  "sequence": aligned_seq,
#     #                  "residue_numbers": residue_numbers}
    
#     # # Appending current dictionary to combined dictionary
#     # sequence_data[pdb_name] = current_entry
    
#     # Store sequence alignment data grouped by PDB
#     if current_pdb not in sequence_data:
#         sequence_data[current_pdb] = []

#     sequence_data[current_pdb].append({
#         "pdb_id": pdb_name,
#         "sequence": aligned_seq,
#         "residue_numbers": residue_numbers
#     })

#     # Assign new numbers to residues
#     for model in structure:
#         for chain in model:
#             seq_index = 0
#             for residue in chain:
#                 if PDB.is_aa(residue):
#                     residue.id = (residue.id[0], residue_numbers[seq_index], residue.id[2])
#                     seq_index += 1

#     # Save updated PDB to repaired folder
#     output_path = os.path.join(unique_chains_dir, f"{pdb_name}.pdb")
#     io.set_structure(structure)
#     io.save(output_path)
#     #print(f"Renumbered PDB saved: {output_path}")



# =========================================
# 5. Update PDB residue numbers based on alignment and record changes
# =========================================

# Dictionaries to store data
sequence_data = {}
mapping_records = []  # to store old/new numbering info

io = PDB.PDBIO()

for record in alignment:
    pdb_name = record.id
    aligned_seq = str(record.seq)
    pdb_file = os.path.join(unique_chains_dir, f"{pdb_name}.pdb")
    structure = parser.get_structure(pdb_name, pdb_file)
    
    # Build new residue numbering taking gaps into account
    residue_numbers = []
    current_number = 1  # adjust to starting residue number if needed
    for aa in aligned_seq:
        if aa != "-":
            residue_numbers.append(current_number)
        current_number += 1
    
    # Extract base PDB ID
    current_pdb = pdb_name.split('_', 1)[0]
    
    # Store alignment info
    if current_pdb not in sequence_data:
        sequence_data[current_pdb] = []

    sequence_data[current_pdb].append({
        "pdb_id": pdb_name,
        "sequence": aligned_seq,
        "residue_numbers": residue_numbers
    })

    # Assign new numbers AND record the mapping
    for model in structure:
        for chain in model:
            seq_index = 0
            for residue in chain:
                if not PDB.is_aa(residue):
                    continue

                old_num = residue.id[1]
                if seq_index < len(residue_numbers):
                    new_num = residue_numbers[seq_index]
                    mapping_records.append({
                        "PDB": current_pdb,
                        "pdb_id": pdb_name,
                        "chain": chain.id,
                        "old_residue_number": old_num,
                        "new_residue_number": new_num
                    })
                    residue.id = (residue.id[0], new_num, residue.id[2])
                    seq_index += 1
                else:
                    print(f"⚠️ Warning: Not enough new residue numbers for {pdb_name} chain {chain.id}")
                    break

    # Save updated PDB
    output_path = os.path.join(unique_chains_dir, f"{pdb_name}.pdb")
    io.set_structure(structure)
    io.save(output_path)


# Save residue renumbering mapping
mapping_df = pd.DataFrame(mapping_records)
mapping_file = os.path.join(output_dir, "residue_number_mapping.csv")
mapping_df.to_csv(mapping_file, index=False)

# =========================================
# 6. Add chain and polymorph info to sequence_data
# =========================================

# Loading polymorph data
polymorph_df = pd.read_csv(os.path.join(root_folder, "COM_and_fibril.csv"), sep=",") 

for pdb_name, entries in sequence_data.items():
    for entry in entries:
        pdb_id = entry["pdb_id"]
        match = polymorph_df.loc[polymorph_df["pdb_id"] == pdb_id]

        if not match.empty:
            # Get all chains and polymorphs for that pdb_id
            chains = match["chain"].tolist()
            polymorphs = match["polymorph"].unique().tolist()

            # If there's only one polymorph, just assign the value; otherwise, store the list
            entry["chain"] = chains
            entry["polymorph"] = polymorphs[0] if len(polymorphs) == 1 else polymorphs
        else:
            # Handle missing matches gracefully
            entry["chain"] = None
            entry["polymorph"] = None
            print(f"Warning: No polymorph match found for {pdb_id}")

# Save aligned sequences
sequence_dict = os.path.join(output_dir, "seq_align_data.json")
with open(sequence_dict, "w") as f:
    json.dump(sequence_data, f, indent=4)


# =========================================
# 7. Applying sequence order to asym_units
# =========================================

asym_units = [file for file in os.listdir(asym_dir) if file.endswith(".pdb")]

for asym in asym_units:
    asym_name = asym.split('_', 1)[0]
    asym_file = os.path.join(asym_dir, asym)

    structure = parser.get_structure(asym_name, asym_file)
    PDBmatch = sequence_data[asym_name]

    for model in structure:
        for chain in model:
            chain_id = chain.id
            matching_entry = None

            for entry in PDBmatch:
                chains = entry.get("chain")
                if chains is None:
                    continue
                if isinstance(chains, list) and chain_id in chains:
                    matching_entry = entry
                    break
                elif isinstance(chains, str) and chain_id == chains:
                    matching_entry = entry
                    break

            if not matching_entry:
                print(f"No residue mapping found for {asym_name} chain {chain_id}")
                continue

            residue_numbers = matching_entry["residue_numbers"]
            #print(f"{asym_name} chain {chain_id}: residue numbers {residue_numbers}")

            seq_index = 0
            for residue in chain:
                if not PDB.is_aa(residue):
                    continue

                # Safety check before indexing
                if seq_index >= len(residue_numbers):
                    print(f"⚠️ Warning: ran out of residue numbers for {asym_name} chain {chain_id}")
                    break

                residue.id = (residue.id[0], residue_numbers[seq_index], residue.id[2])
                seq_index += 1

    # Save updated PDB to repaired folder
    output_path = os.path.join(asym_dir, f"{asym_name}_asym_unit.pdb")
    io.set_structure(structure)
    io.save(output_path)




















