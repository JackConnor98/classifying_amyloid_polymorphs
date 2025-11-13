import os
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import re

print("Running: plot_ordered_residues.py")



def get_residue_ranges(pdb_path):
    residues = set()

    with open(pdb_path, "r") as pdb_file:
        for line in pdb_file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    res_num = int(line[22:26].strip())
                    residues.add(res_num)
                except ValueError:
                    continue

    residues = sorted(residues)
    if not residues:
        return []

    ranges = []
    start = prev = residues[0]

    for res in residues[1:]:
        if res == prev + 1:
            prev = res
        else:
            ranges.append((start, prev))
            start = prev = res
    ranges.append((start, prev))

    return ranges


def alphanum_key(s):
    """Turn a string into a list of ints and strings for natural sorting."""
    return [int(text) if text.isdigit() else text.lower() for text in re.split('([0-9]+)', s)]



# =============================================================================
# Load and clean data
# =============================================================================
data = pd.read_csv(os.path.join("Output", "selected_pdbs_metadata.csv"))

# Remove NMR structures
data = data[data["Method"] != "ssNMR"]

# Keep relevant columns
residues = data[["PDB ID", "Residues Ordered"]].copy()
residues.columns = ["pdb", "ordered_residues"]



### Adding in local structures
local_pdbs = [f for f in os.listdir("Local") if f.endswith(".pdb")]
local_dir = "Local"
local_pdbs = [f for f in os.listdir(local_dir) if f.endswith(".pdb")]

# Collect the residue ranges for each PDB
local = []
for pdb_file in local_pdbs:
    pdb_path = os.path.join(local_dir, pdb_file)
    residue_ranges = get_residue_ranges(pdb_path)
    formatted_ranges = [f"{start}-{end}" if start != end else f"{start}" for start, end in residue_ranges]
    # Add to data with filename minus .pdb
    local.append({
        "pdb": os.path.splitext(pdb_file)[0],
        "ordered_residues": ", ".join(formatted_ranges)
    })

# Convert to DataFrame
local_df = pd.DataFrame(local)

residues = pd.concat([residues, local_df], ignore_index=True)




# Split by commas (indicating discontinuities)
residues["ordered_residues"] = residues["ordered_residues"].str.split(",")

# Flatten the lists into long form
residues = residues.explode("ordered_residues").reset_index(drop=True)

# Clean messy characters
residues["ordered_residues"] = (
    residues["ordered_residues"]
    .astype(str)
    .str.replace(r'[c()"()]', "", regex=True)
    .str.strip()
)



# Split comma-separated regions
residues["ordered_residues"] = residues["ordered_residues"].str.split(",")

# Expand into separate rows
residues = residues.explode("ordered_residues").reset_index(drop=True)

# Clean up stray characters
residues["ordered_residues"] = (
    residues["ordered_residues"]
    .astype(str)
    .str.replace(r'[c()"()]', "", regex=True)
    .str.strip()
)

# Split each stretch (e.g. 20-96) into start and end
split_cols = residues["ordered_residues"].str.split("-", expand=True)
split_cols.columns = ["start", "end"]

# Convert to numeric (ignore errors from NaN or empty values)
split_cols = split_cols.apply(pd.to_numeric, errors="coerce")

# Add a group ID for each continuous region per PDB
residues = pd.concat([residues[["pdb"]], split_cols], axis=1)
residues["group"] = residues.groupby("pdb").cumcount() + 1

# Expand start/end into separate rows
residues_long = pd.melt(
    residues,
    id_vars=["pdb", "group"],
    value_vars=["start", "end"],
    var_name="boundary",
    value_name="ordered_residues"
).drop(columns=["boundary"])

# Convert to numeric
residues_long["ordered_residues"] = pd.to_numeric(residues_long["ordered_residues"], errors="coerce")
residues_long = residues_long.dropna(subset=["ordered_residues"]).reset_index(drop=True)


residues_long["pdb_group"] = residues_long["pdb"].astype(str) + "_" + residues_long["group"].astype(str)

# =============================================================================
# Plotting
# =============================================================================
pdbs = sorted(residues_long["pdb"].unique(), key=alphanum_key)
num_pdbs = len(pdbs)
protein_names = data["PDB ID"].unique()

# Set base width per label (tweak this depending on how dense your labels are)
width_per_label = 0.2  # inches per PDB label
min_width = 8  # minimum plot width

# Calculate final width
plot_width = max(min_width, num_pdbs * width_per_label)
plot_width = min(plot_width, 20)  # cap maximum width to 25 inches

# Set base font size and minimum font size
base_font_size = 14
min_font_size = 8

# Scale font size down as number of labels increases
# The more labels, the smaller the text
x_axis_font_size = max(min_font_size, base_font_size - 0.05*num_pdbs)

# Create an interaction column
residues_long["pdb_group_interaction"] = residues_long["pdb"].astype(str) + "_" + residues_long["group"].astype(str)

### Creating protein specific annotations ###
protein_type = data["Protein"].str.lower()  # lowercase for easier matching

# Initialize annotation and label DataFrames
annotations = pd.DataFrame(columns=["xmin","xmax","ymin","ymax","fill"])
labels = pd.DataFrame(columns=["x","y","label"])

# Setting plot dimensions
fig, ax = plt.subplots(figsize=(max(8, num_pdbs * 0.2), 8))

# Plot each PDB as a line with dots for residues
for pdb, group_df in residues_long.groupby("pdb_group_interaction"):
    ax.plot(group_df["pdb"], group_df["ordered_residues"], lw=2, color = "black")
    ax.scatter(group_df["pdb"], group_df["ordered_residues"], s=30, color="black", zorder=3)

for pdb, group_df in residues_long.groupby("pdb"):
    # Sort residues numerically
    group_df = group_df.sort_values("ordered_residues")
    # Get unique continuous regions
    regions = group_df.groupby("pdb_group_interaction")["ordered_residues"].agg(["min", "max"]).reset_index()
    # Draw dashed lines between discontinuous regions
    for i in range(len(regions)-1):
        end_prev = regions.loc[i, "max"]
        start_next = regions.loc[i+1, "min"]
        ax.plot(
            [pdb, pdb],                  # same x-coordinate for this PDB
            [end_prev, start_next],       # vertical line connecting segments
            linestyle="--", color="black", lw=2
        )

# Tidy the axis labels
ax.set_xlabel("PDB", fontsize=24, fontweight="bold")
ax.set_ylabel("Residue Position", fontsize=24, fontweight="bold")
ax.set_xticks(range(len(pdbs)))
ax.set_xticklabels(pdbs, rotation=90, ha="center", fontsize=x_axis_font_size)
ax.tick_params(axis="y", labelsize=18)
ax.set_xlim(-0.5, num_pdbs - 0.5) 


# Determine protein type for annotations
protein_type = data["Protein"].str.lower()
ann_regions = []

if protein_type.str.contains("synuclein").any():
    ann_regions = [
        (1, 60, "darkblue", "N-Term"),
        (61, 95, "red", "NAC"),
        (96, 140, "green", "C-Term")
    ]
    rect_height = 15  
    rect_width = max(2, num_pdbs * 0.1) # width scales with number of PDBs

elif protein_type.str.contains("amyloid-").any():
    ann_regions = [
        (1, 16, "darkblue", "N-Term"),
        (17, 21, "black", "Hydrophobic\nCore"),
        (24, 27, "red", "Turn Region"),
        (28, 35, "orange", "Hydrophobic\nRegion"),
        (36, 42, "green", "C-Term")
    ]
    rect_height = 5  
    rect_width = max(2, num_pdbs * 0.15) # width scales with number of PDBs

elif protein_type.str.contains("tau").any():
    ann_regions = [
        (243, 273, "red", "R1"),
        (274, 304, "green", "R2"),
        (305, 335, "orange", "R3"),
        (336, 367, "darkblue", "R4")
    ]
    rect_height = 15  
    rect_width = max(2, num_pdbs * 0.1) # width scales with number of PDBs

# Extend y-axis to include both residues and annotations
all_y = residues_long["ordered_residues"].values
if ann_regions:
    ann_ymin = min([ymin for ymin, ymax, _, _ in ann_regions])
    ann_ymax = max([ymax for ymin, ymax, _, _ in ann_regions])
    y_min = min(all_y.min(), ann_ymin)
    y_max = max(all_y.max(), ann_ymax)
else:
    y_min, y_max = all_y.min(), all_y.max()
ax.set_ylim(y_min - 5, y_max + 5)

# Main-axis rectangles (full width)
for (ymin, ymax, color, _) in ann_regions:
    rect_main = Rectangle(
        (-0.5, ymin),          # start at left edge of x-axis
        num_pdbs - (-0.5),     # width spans all PDBs
        ymax - ymin,
        facecolor=color,
        edgecolor="black",
        alpha=0.5,
        clip_on=True           # clipped to axes
    )
    ax.add_patch(rect_main)

# Draw annotation rectangles outside the axes with labels inside
x_offset = max(1, num_pdbs * 0.01)   # 1% of the total x-range

for (ymin, ymax, color, label) in ann_regions:
    # Centre the rectangle at the midpoint of the region
    y_mid = (ymin + ymax) / 2
    rect = Rectangle(
        (num_pdbs + x_offset, y_mid - rect_height/2),  # centre rectangle
        rect_width,
        rect_height,
        facecolor=color,
        edgecolor="black",
        alpha=0.5,
        clip_on=False
    )
    ax.add_patch(rect)
    
    # Text remains centred
    ax.text(
        num_pdbs + x_offset + rect_width/2, 
        y_mid, 
        label, 
        va="center", 
        ha="center", 
        fontsize=12, 
        fontweight="bold", 
        clip_on=False
    )

# Panel Grid
ax.grid(True, linestyle="--", color="black", alpha = 0.5, lw = 0.25)

# Save plot
plt.savefig("Output/Ordered_Residues.png", dpi=300, bbox_inches='tight')
print("Plot saved to Output/Ordered_Residues.png")
