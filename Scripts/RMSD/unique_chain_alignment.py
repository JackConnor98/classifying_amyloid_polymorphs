import pymol
from pymol import cmd
import os 
import sys
import csv
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

output_dir = os.path.join("Output", "RMSD")

if not os.path.exists(output_dir):
    os.mkdir(output_dir)

save_path = os.path.join("Output", "RMSD", "data")

if not os.path.exists(save_path):
    os.mkdir(save_path)

output_file = os.path.join(save_path, "rmsd_output.txt")
csv_file = os.path.join(save_path, "pairwise_rmsd_by_reference.csv")

# Finding all PDB files in the specified directory
pdb_path = os.path.join("Output", "PDBs", "unique_chains")
pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]
object_names = [i.split(".")[0] for i in pdb_files]

# print(len(pdb_files))
# print(len(object_names))

### Removing PDBs with poor resolution ###
# Loading high resolution PDB IDs
high_res_df = pd.read_csv(os.path.join("Output", "Validation", "high_resolution_pdb_ids.csv"), sep=",")
high_res_ids = set(high_res_df['pdb_id'].astype(str))  # Just in case you're storing ints like a noob
# Filter pdb_files and object_names based on high_res_ids
filtered = [(f, name) for f, name in zip(pdb_files, object_names) if name in high_res_ids]
pdb_files, object_names = zip(*filtered) if filtered else ([], [])

# print(len(pdb_files))
# print(len(object_names))


# Initialise variable to store rmsd abd ca positions for each comparison
pairwise_rmsd_data = []
data = []

# Storing a list of names already used as the reference to prevent duplicate comparisons
already_used = []

pymol.pymol_argv = ['pymol', '-qc']
pymol.finish_launching()

for i in range(len(pdb_files)):
    ref_path = os.path.join(pdb_path, pdb_files[i])
    ref_name = pdb_files[i].split(".")[0]
    cmd.load(ref_path, ref_name)

    print("Analysing reference structure:", ref_name)

    #for j in range(i + 1, len(pdb_files)):
    for j in range(len(pdb_files)):    

        mob_path = os.path.join(pdb_path, pdb_files[j])
        mob_name = pdb_files[j].split(".")[0]
        mob_obj_name = f"mob_{mob_name}" # create a unique name for the mobile object to avoid error when doing a self alignment

        cmd.load(mob_path, mob_obj_name)
        
        # Align mobile to reference and get RMSD
        rmsd = cmd.align(ref_name, mob_obj_name)[0]

        pairwise_rmsd_data.append({"ref_name": ref_name,
                                   "mob_name": mob_name, 
                                   "rmsd": rmsd})

        model_ref = cmd.get_model(f"{ref_name} and name CA")
        model_mob = cmd.get_model(f"{mob_obj_name} and name CA")

        # Make dicts with keys = (resi, chain), values = atom object
        ref_atoms = { (a.resi): a for a in model_ref.atom }
        mob_atoms = { (a.resi): a for a in model_mob.atom }

        # Get the full union of all residues seen in either structure
        all_keys = sorted(set(ref_atoms.keys()) | set(mob_atoms.keys()))

        for key in all_keys:
            a_ref = ref_atoms.get(key)
            a_mob = mob_atoms.get(key)

            row = {
            "resi": key[0],
            "ref_name": ref_name,
            "ref_resn": a_ref.resn if a_ref else "NA",
            "ref_x": a_ref.coord[0] if a_ref else "NA",
            "ref_y": a_ref.coord[1] if a_ref else "NA",
            "ref_z": a_ref.coord[2] if a_ref else "NA",
            "mob_name": mob_name,
            "mob_resn": a_mob.resn if a_mob else "NA",
            "mob_x": a_mob.coord[0] if a_mob else "NA",
            "mob_y": a_mob.coord[1] if a_mob else "NA",
            "mob_z": a_mob.coord[2] if a_mob else "NA",
            }

            data.append(row)

        # Unload the mobile structure after alignment
        cmd.delete(mob_obj_name)

cmd.quit()

# Creating an RMSD dataframe
pairwise_rmsd_df = pd.DataFrame(pairwise_rmsd_data)
pairwise_rmsd_df.loc[pairwise_rmsd_df["ref_name"] == pairwise_rmsd_df["mob_name"], "rmsd"] = 0

# Generating a pandas data frame for ca positions
df = pd.DataFrame(data)

###########################################
### Calculating Distance and NA Penalty ###
###########################################

df["distance"] = "NA"  # default to "NA" for everything

# Calculating NA penalty
filtered_df = pairwise_rmsd_df[pairwise_rmsd_df["ref_name"] != pairwise_rmsd_df["mob_name"]] # Removing self comparisons as they ruin the mean
mean_rmsd = filtered_df["rmsd"].mean()
std_rmsd = filtered_df["rmsd"].std()

# Getting settings for NA penalty defined in run_analysis.sh
penalty = int(sys.argv[1])

if penalty == 0:
    na_penalty = "NA"  # No penalty
else:
    # Calculate NA penalty based on mean and standard deviation
    na_penalty = mean_rmsd + (penalty * std_rmsd)

print(f"NA Penalty = {na_penalty}")

for i in range(len(df)):

    row = df.iloc[i]

    if "NA" in [row["ref_x"], row["ref_y"], row["ref_z"],
                row["mob_x"], row["mob_y"], row["mob_z"]]:

        # Assigning penalty for non-overlapping regions
        distance = na_penalty

    else:
        # Calculating distance between C-alphas
        distance = np.sqrt(
            (float(row["ref_x"]) - float(row["mob_x"]))**2 +
            (float(row["ref_y"]) - float(row["mob_y"]))**2 +
            (float(row["ref_z"]) - float(row["mob_z"]))**2)
        
    df.at[i, "distance"] = distance

# Saving ca coordinate data
df.to_csv(os.path.join(save_path, "ca_distances.csv"), index=False)

###############################
### Calculating RMSD Scores ###
###############################

# Filter out rows with NA distances
df_valid = df[df["distance"] != "NA"].copy()
df_valid["distance"] = df_valid["distance"].astype(float)

# Group by ref and mob
grouped = df_valid.groupby(["ref_name", "mob_name"])

# Calculate RMSD per group
rmsd_df = grouped["distance"].apply(lambda x: np.sqrt(np.mean(x**2))).reset_index()
rmsd_df.columns = ["ref_name", "mob_name", "rmsd"]
rmsd_df.loc[rmsd_df["ref_name"] == rmsd_df["mob_name"], "rmsd"] = 0

rmsd_df.to_csv(os.path.join(save_path, "pairwise_rmsd.csv"), index=False)


# Plotting the impact of NA penalty
merged_df = pairwise_rmsd_df.merge(
    rmsd_df,
    left_on=["ref_name", "mob_name"],
    right_on=["ref_name", "mob_name"],
    suffixes=('', '_with_na')
)


plt.figure(figsize=(6,6))
plt.scatter(merged_df['rmsd'], merged_df['rmsd_with_na'], alpha=0.7)

# Plot y = x line
lims = [
    np.min([merged_df['rmsd'].min(), merged_df['rmsd_with_na'].min()]),
    np.max([merged_df['rmsd'].max(), merged_df['rmsd_with_na'].max()])
]
plt.plot(lims, lims, "r--")

plt.xlabel('rmsd')
plt.ylabel('rmsd_with_na')
plt.title('RMSD vs RMSD with NA penalty')
plt.grid(True)

plt.savefig(os.path.join(save_path, "rmsd_and_na_penalty.png"), dpi=300, bbox_inches='tight')