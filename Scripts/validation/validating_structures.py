import os
import pandas as pd
import numpy as np
import gzip
import xml.etree.ElementTree as ET
from plotnine import *
from scipy.stats import gaussian_kde

############################################################################################

# Locate all PDB Q-score data in folder that need to be looped through
folder_path = os.path.join("Output", "Validation")

data_path = os.path.join(folder_path, "data")

# Get a list of all .csv files in the folder
xml_files = [
    os.path.join(data_path, f)
    for f in os.listdir(data_path)
    if f.endswith(".gz")
]

############################################################################################

save_path = os.path.join(folder_path, "plots")
single_pdb_path = os.path.join(save_path, "single_pdbs")

# Creating Directories
if not os.path.exists(save_path):
    os.mkdir(save_path)
if not os.path.exists(single_pdb_path):
    os.mkdir(single_pdb_path)

############################################################################################

# Importing metadata
selected_pdbs_metadata = pd.read_csv(os.path.join("Output", "selected_pdbs_metadata.txt"), sep="\t")

# Fixing names
selected_pdbs_metadata = selected_pdbs_metadata.rename(columns={
    "PDB ID": "PDB",
    "Resol- ution (Å)": "resolution"
})

# Selecting columns
resolution_df = selected_pdbs_metadata[["PDB", "resolution"]]

############################################################################################

##################################################################################
### Creating a dataframe containing the Q-scores for every residue of each PDB ###
##################################################################################

validation_data = pd.DataFrame()

for file in xml_files:
    # Get filename
    filename = os.path.basename(file)

    # Read XML
    try:
        with gzip.open(file, 'rt', encoding='utf-8') as f:  # 'rt' = read text mode
            tree = ET.parse(f)
            root = tree.getroot()
    except ET.ParseError:
        print(f"Skipping {file}: invalid XML.")
        continue
    except OSError:
        print(f"Skipping {file}: not a valid gzip file.")
        continue
    
    # Find all ModelledSubgroup nodes
    nodes = root.findall(".//ModelledSubgroup")
    
    # Get all unique attribute names
    all_attributes = set()
    for node in nodes:
        all_attributes.update(node.attrib.keys())
    all_attributes = list(all_attributes)
    
    # Check if Q_score attribute exists
    if "Q_score" in all_attributes:
        rows = []
        for node in nodes:
            attrs = node.attrib.copy()
            # Ensure all columns are present
            for col in all_attributes:
                if col not in attrs:
                    attrs[col] = None
            rows.append(attrs)
        
        # Convert to dataframe
        df = pd.DataFrame(rows)
        
        # Add PDB column (everything before first underscore)
        df["pdb"] = filename.split("_")[0]
        
        # Append to main dataframe
        validation_data = pd.concat([validation_data, df], ignore_index=True)
    else:
        print(f"{filename.split('_')[0]} : No Q-scores found")


# Selecting key columns 
q_score_df = validation_data[["pdb", "Q_score", "chain", "resnum", "resname"]]

# Changing resnum to resno to better fit other data
q_score_df = q_score_df.rename(columns={
    "pdb": "PDB",
    "resnum": "resno",
    "chain": "published_chain"
})

# Converting columns to numeric
q_score_df[["Q_score", "resno"]] = q_score_df[["Q_score", "resno"]].apply(pd.to_numeric, errors='coerce')


# Removing NAs (This will remove ssNMR structures as Q-scores are for cryoEM)
q_score_df = q_score_df.dropna(subset=["Q_score"])

# Character vector containing the capitalized 3-letter codes for all amino acids in all capitals
amino_acids = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
]

# Removing non-amino acid components
q_score_df = q_score_df[q_score_df["resname"].isin(amino_acids)]


### Issue - renaming chains during PDB handling means they dont match up to Q-score data ###
# Need to store the original chain IDs from the PDB files and match them to the Q-score data

# Importing chain mapping data
chain_mapping_df = pd.read_csv(os.path.join("Output", "PDBs", "chain_mapping.csv"))

# Merging chain mapping with q_score_df
q_score_df = pd.merge(chain_mapping_df, q_score_df, on=["PDB", "published_chain"])

# Renaming modified_chain to chain
q_score_df = q_score_df.rename(columns={"modified_chain": "chain"})

### Adding in pdb_id ###

# Importing data
pdb_info = pd.read_csv(os.path.join("Output", "PDBs", "COM_and_fibril.csv"), sep=",")

# Selecting columns of interest
pdb_info = pdb_info[["pdb_id", "PDB", "chain", "fibril"]]

# Merge pdb_info and q_score_df on 'pdb_id' and 'chain'
q_score_df = pd.merge(pdb_info, q_score_df, on=["PDB", "chain"])

# Saving Q-Score data frame
q_score_df.to_csv(os.path.join(folder_path, "Q_Scores.csv"), index=False)

##########################################
### Filtering Out Poorly Resolved PDBs ###
##########################################

# Merging resolution with q_scores
q_score_df = pd.merge(q_score_df, resolution_df, on = "PDB")

# Calculating the mean and standard deviation for the Q-score across all structures
# Threshold for Q-score will me set to mean - 1SD
mean_Q_score = np.mean(q_score_df["Q_score"])
sd_Q_score = np.std(q_score_df["Q_score"])
Q_score_threshold = mean_Q_score - sd_Q_score


####################################################
### Calculating the mean Q-score for each pdb_id ###
####################################################
PDB_q_score = (
    q_score_df.groupby("pdb_id", as_index=False)
              .agg(mean_Q_score=("Q_score", "mean"))
)

# Putting in descending order
PDB_q_score = PDB_q_score.sort_values(by="mean_Q_score", ascending=False)

# Making pdb_id an ordered categorical variable for plotting
PDB_q_score["pdb_id"] = pd.Categorical(
    PDB_q_score["pdb_id"],
    categories=PDB_q_score["pdb_id"].tolist(),
    ordered=True
)

# Saving the mean Q-score for each PDB
PDB_q_score.to_csv(os.path.join(folder_path, "mean_PDB_Q_score.csv"), index=False)

# Plottin the mean Q-score for each PDB
p = (
    ggplot(PDB_q_score, aes(x="pdb_id", y="mean_Q_score"))
    + geom_hline(yintercept = mean_Q_score, colour = "black", size = 0.5)
    + geom_hline(yintercept=Q_score_threshold, colour="red", linetype="dashed", size=0.5)
    + geom_point(size=1)
    + labs(
        x="PDB ID",
        y="Mean Q-score",
        title=""
    )
    + theme_classic()
    + theme(
        panel_border=element_rect(linewidth = 1, fill = None),
        panel_grid_major = element_line(linewidth = 0.25, colour = "grey", linetype = "dashed"),
        axis_title=element_text(size=18, colour="black", face="bold"),
        axis_text_y = element_text(size = 14, colour = "black"),
        axis_text_x = element_text(size = 6, colour = "black", 
                                   angle = 90, vjust = 0.25)
    )
)

p.save(os.path.join(save_path, "Mean_PDB_Q_scores.png"), height=4, width=8, dpi=300)

###############################################################
### Plotting each PDB individually for quality/sense checks ###
###############################################################

pdb_names = pd.unique(q_score_df["pdb_id"])

max_y = np.max(q_score_df["Q_score"])
min_y = np.min(q_score_df["Q_score"])
max_x = np.max(q_score_df["resno"])
min_x = np.min(q_score_df["resno"])

for i in pdb_names:
    
    # Filtering current PDB
    x = q_score_df[q_score_df["pdb_id"] == i].copy()

    # Colouring outliers
    x["outliers"] = np.where(x["Q_score"] < Q_score_threshold, "outlier", "good")

    # Getting Mean Q-score for current PDB
    current_mean_Q_score = PDB_q_score["mean_Q_score"][PDB_q_score["pdb_id"] == i].iloc[0]

    # Round to 3 significant figures
    current_mean_Q_score = float(f"{current_mean_Q_score:.3g}")  

    # Get unique resolution for the current PDB
    current_resolution = q_score_df.loc[q_score_df["pdb_id"] == i, "resolution"].unique()

    # Plotting current PDB Q-Scores
    p = (
        ggplot(x, aes(x = "resno", y = "Q_score", colour = "outliers"))
        + geom_hline(yintercept = mean_Q_score, colour = "black", size = 0.5)
        + geom_hline(yintercept=Q_score_threshold, colour="red", linetype="dashed", size=0.5)
        + geom_point(size=1)
        + scale_colour_manual(values = {"good": "black", "outlier": "grey"})
        + labs(
            title = f"PDB ID = {i} | Resolution = {current_resolution[0]} Å | Mean Q-score = {current_mean_Q_score}",
            x="Residue Number",
            y="Q-score",
        )
        + xlim(min_x, max_x) 
        + ylim(min_y, max_y)
        + theme_classic()
        + theme(
            panel_border=element_rect(linewidth = 1, fill = None),
            panel_grid_major = element_line(linewidth = 0.25, colour = "grey", linetype = "dashed"),
            plot_title=element_text(size=14, colour="black", face="bold", ha="center"),
            axis_title=element_text(size=18, colour="black", face="bold"),
            axis_text = element_text(size = 14, colour = "black"),
            legend_position = "none"
        )
    )

    p.save(os.path.join(single_pdb_path, f"{i}_Q_score.png"), height=4, width=8, dpi=300)


######################################
### Removing Low Mean Q-score PDBs ###
######################################

# Selecting PDBs with a mean Q-score > Q_score_threshold
PDB_q_score["pdb_id"] = PDB_q_score["pdb_id"].astype(str)
high_resolutuion_PDBs = PDB_q_score["pdb_id"][PDB_q_score["mean_Q_score"] > Q_score_threshold].tolist()

# Saving a list of high_resolution PDBs
high_resolutuion_PDBs_df = PDB_q_score[PDB_q_score["mean_Q_score"] > Q_score_threshold]
high_resolutuion_PDBs_df.to_csv(os.path.join(folder_path, "high_resolution_pdb_ids.csv"), index=False)

#####################################
### Removing Low Q-score Residues ###
#####################################

# Removing low resolution PDBs from q_score_df
filtered_df = q_score_df[q_score_df["pdb_id"].isin(high_resolutuion_PDBs)].copy()
filtered_df = filtered_df[["pdb_id", "chain", "resno", "Q_score"]]

# Recalculating mean and standard deviation with low resolution PDBs removed
mean_Q_score = np.mean(filtered_df["Q_score"])
sd_Q_score = np.std(filtered_df["Q_score"])
Q_score_threshold = mean_Q_score - sd_Q_score

# Selecting good resolution residues
good_resolution = filtered_df[filtered_df["Q_score"] >= Q_score_threshold]

################################################
### Visualising Q-score Vairance in Residues ###
################################################

# Extract Q-scores as a numpy array
q_scores = filtered_df["Q_score"].values

# Calculate density
kde = gaussian_kde(q_scores)
x = np.linspace(q_scores.min(), q_scores.max(), 1000)  # 1000 points across range
y = kde(x)

# Create a DataFrame from the density data
density_df = pd.DataFrame({"x": x, "y": y})

# Saving density data
density_df.to_csv(os.path.join(folder_path, "Q_score_density_data.csv"), index=False)

# Split density_df into two for coloring
density_below = density_df[density_df["x"] <= Q_score_threshold]
density_above = density_df[density_df["x"] > Q_score_threshold]

# Plotting a density plot of Q-scores with areas colored based on threshold
p = (
    ggplot() +
    
    # Area for x <= threshold
    geom_area(density_below, aes(x="x", y="y"), fill="grey", alpha=0.75, colour="black", size=0.8) +
    
    # Area for x > threshold
    geom_area(density_above, aes(x="x", y="y"), fill="blue", alpha=0.75, colour="black", size=0.8) +
    
    # Vertical lines
    geom_vline(xintercept=mean_Q_score, colour="black", linetype="solid", size=1) +
    geom_vline(xintercept=Q_score_threshold, colour="red", linetype="dashed", size=1) +
    
    # Labels
    labs(x="Q Score", y="Density") +
    
    # X and Y scales
    scale_x_continuous(breaks = np.linspace(0, 1, num=11), expand=(0, 0, 0, 0)) +
    scale_y_continuous(expand=(0, 0, 0.1, 0.)) +

    # Theme
    theme_classic() +
    theme(
        panel_border=element_rect(linewidth=1, fill=None),
        axis_title=element_text(size=18, face="bold", colour="black"),
        axis_text=element_text(size=14, colour="black")
    )
)

# Save plot
#p.save(os.path.join(save_path, "Q_score_density.png"), height=4, width=6, dpi=300)

### Stacked Bar Chart ###

# Count occurrences of each resno in good_resolution
good_resolution_counts = (
    good_resolution.groupby("resno")
    .size()
    .reset_index(name="good_count")
)

# Count occurrences of each resno in filtered_df
total_counts = (
    filtered_df.groupby("resno")
    .size()
    .reset_index(name="total_count")
)

# Merging good and total counts
count_df = pd.merge(total_counts, good_resolution_counts, on="resno", how="left")

# Fill NaN values in good_count with 0 (in case some residues have no good counts)
count_df["good_count"] = count_df["good_count"].fillna(0)

print(count_df)

# Calculate % good and bad count
count_df["percentage"] = count_df["good_count"] / count_df["total_count"] * 100
count_df["bad_count"] = count_df["total_count"] - count_df["good_count"]

# Transform to long format
long_count_df = count_df.melt(
    id_vars=["resno", "percentage"],    
    value_vars=["good_count", "bad_count"],  
    var_name="count_type",               
    value_name="count"                   
)

# Saving bar chart data
long_count_df.to_csv(os.path.join(folder_path, "resolution_bar_chart_data.csv"), index=False)

# Plotting stacked bar chart
p = (
    ggplot(long_count_df, aes(x="resno", y="count", fill="count_type"))
    + geom_bar(stat = "identity", alpha=0.75, colour="black")
    + scale_fill_manual(values={"good_count": "blue", "bad_count": "grey"}, labels={"good_count": "Good", "bad_count": "Poor"})
    + labs(x="Residue Number", y="Count", fill="Resolution")
    + scale_x_continuous(breaks=np.arange(0, max_x+1, 5), expand=(0.01, 0, 0.01, 0))
    + scale_y_continuous(expand=(0, 0, 0.05, 0.))
    + theme_classic()
    + theme(
        panel_border=element_rect(linewidth=1, fill=None),
        axis_title=element_text(size=18, face="bold", colour="black"),
        axis_text=element_text(size=14, colour="black"),
        legend_title=element_text(size=14, face="bold", colour="black"),
        legend_text=element_text(size=10, colour="black")
    )
)

p.save(os.path.join(save_path, "resolution_bar_chart.png"), height=4, width=8, dpi=300)


# Important Note:
# Possible issue with chains having different number of residues
# A single residue may occur above and below the Q-score threshold in different chains of the same fibril
# Current solution - if the residue occurs at good resolution at least once in a fibril - it is included

# Adding back in other information
pdb_info_no_chains = pdb_info[["pdb_id", "PDB", "fibril"]].drop_duplicates()
tmp = pd.merge(good_resolution, pdb_info_no_chains, on="pdb_id")

# Selecting residues that are present in at least one chain of a fibril at good resolution
high_resolution_residues = tmp[['pdb_id', 'fibril', 'resno']].drop_duplicates()

# Putting in ascending residue order
high_resolution_residues = high_resolution_residues.sort_values(
    by=["pdb_id", "fibril", "resno"],
    ascending=[True, True, True]
).reset_index(drop=True)

# Saving high resolution residues
high_resolution_residues.to_csv(os.path.join(folder_path, "high_resolution_residues.csv"), index=False)