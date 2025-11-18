import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, Superimposer
import math
import os
from statistics import mode

def unit_vector(v):
    norm = np.linalg.norm(v)
    return v / norm if norm != 0 else v


def extract_twist_and_rise(chain1, chain2):
    """
    Compute twist, rise, and handedness between two chains
    using backbone CA atom superposition.
    """

    # Extract corresponding CA atoms
    atoms1 = [a for a in chain1.get_atoms() if a.get_id() == "CA"]
    atoms2 = [a for a in chain2.get_atoms() if a.get_id() == "CA"]

    if len(atoms1) != len(atoms2) or len(atoms1) == 0:
        raise ValueError("Mismatched or empty atom sets for alignment.")

    # Superimpose
    sup = Superimposer()
    sup.set_atoms(atoms1, atoms2)
    R, t = sup.rotran

    # Identify rotation axis (eigenvector where eigenvalue ≈ 1)
    w, v = np.linalg.eig(R)
    axis = v[:, np.isclose(w, 1.0)].flatten()

    axis = unit_vector(axis)

    # Compute twist angle from rotation matrix
    cos_angle = (np.trace(R) - 1) / 2
    cos_angle = max(min(cos_angle, 1.0), -1.0)
    angle_deg = math.degrees(math.acos(cos_angle)) * -1 # To account for direcetion

    # Compute rise along the axis
    rise = np.dot(t, axis)

    # Determine handedness
    handedness = "right" if angle_deg > 0 else "left"

    # Compute half-pitch = rise / twist_radians
    half_pitch = (abs(rise) * 180) / abs(angle_deg)

    return angle_deg, abs(rise), handedness, half_pitch


# ---------------------------------------------------------
# Load CSVs
# ---------------------------------------------------------

COM = pd.read_csv("Output/PDBs/COM_and_fibril.csv")
pairs = pd.read_csv("Output/PDBs/min_COM_distances.csv")
pairs.columns = [c.strip() for c in pairs.columns]


# ---------------------------------------------------------
# Process each chain pair and extract twist/rise
# ---------------------------------------------------------

parser = PDBParser(QUIET = True)

results = []

for pdb_id in pairs["PDB"].unique():
    pdb_file = f"Output/PDBs/published_structure/{pdb_id}.pdb"

    if not os.path.exists(pdb_file):
        print(f"WARNING: File {pdb_file} not found. Skipping.")
        continue

    structure = parser.get_structure(pdb_id, pdb_file)

    these_pairs = pairs[pairs["PDB"] == pdb_id]

    for _, row in these_pairs.iterrows():
        chain1_id = str(row["chain_1"]).strip()
        chain2_id = str(row["chain_2"]).strip()


        # Look up fibril IDs for each chain
        chain1_fibril = COM.loc[(COM["pdb_id"] == pdb_id) & (COM["chain"] == chain1_id), "fibril"]
        chain2_fibril = COM.loc[(COM["pdb_id"] == pdb_id) & (COM["chain"] == chain2_id), "fibril"]

        # Skip if fibril info is missing or fibrils differ
        if chain1_fibril.empty or chain2_fibril.empty or chain1_fibril.iloc[0] != chain2_fibril.iloc[0]:
            continue

        try:
            chain1 = structure[0][chain1_id]
            chain2 = structure[0][chain2_id]
        except KeyError:
            print(f"Missing chain {chain1_id} or {chain2_id} in {pdb_id}. Skipping.")
            continue

        try:
            twist, rise, handed, half_pitch = extract_twist_and_rise(chain1, chain2)
        except Exception as e:
            print(f"Alignment failed for {pdb_id} {chain1_id}-{chain2_id}: {e}")
            continue

        results.append({
            "pdb": pdb_id,
            "chain_1": chain1_id,
            "chain_2": chain2_id,
            "twist": twist,
            "rise": rise,
            "half_pitch": half_pitch,
            "handedness": handed
        })


def majority_hand(values):
    try:
        return mode(values)
    except:
        return values.iloc[0]  # fallback if tie

# ---------------------------------------------------------
# Save output
# ---------------------------------------------------------

results_df = pd.DataFrame(results)

# One summary row per PDB
summary_df = (
    results_df
        .groupby("pdb", as_index = False)
        .agg({
            "twist": "mean",
            "rise": "mean",
            "half_pitch": "mean",
            "handedness": majority_hand
        })
)

summary_df.to_csv(os.path.join("Output", "PDBs", "twist_rise_summary.csv"), index = False)
print(summary_df)
