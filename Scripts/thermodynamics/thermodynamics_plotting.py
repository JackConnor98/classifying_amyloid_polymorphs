import sys
import os
import re
import pandas as pd
import numpy as np
from plotnine import *

# Read command line arguments
if len(sys.argv) > 1: 
    try: 
        remove_poorly_resolved = int(sys.argv[1]) 
    except ValueError: 
        print("Warning: invalid argument — defaulting to REMOVE poorly resolved residues") 
        remove_poorly_resolved = 1 
else: 
    remove_poorly_resolved = 1 
    print("No input given, defaulting to REMOVE poorly resolved residues")

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
metadata = pd.read_csv(os.path.join("Output", "selected_pdbs_metadata.csv"))
chain_fibril = pd.read_csv(os.path.join("Output", "PDBs", "fibrils_extended", "chain_fibril.csv"))
exterior_chains = pd.read_csv(os.path.join("Output", "PDBs", "fibrils_extended", "exterior_chains.csv"))
com_and_fibril = pd.read_csv(os.path.join("Output", "PDBs", "COM_and_fibril.csv"))
high_resolution_residues = pd.read_csv(os.path.join("Output", "Validation", "high_resolution_residues.csv"))
cluster_groups = pd.read_csv(os.path.join("Output", "RMSD", "data", "RMSD_cluster_groups.csv"))

######################################
### Removing non-cryoEM structures ###
######################################

# Extracting PDB IDs with Method "cryoEM" from metadata
cryoem_ids = metadata[metadata["Method"] == "cryoEM"]["PDB ID"].tolist()
# Removing non-cryoEM structures from df
df = df[df["PDB"].isin(cryoem_ids)].reset_index(drop=True)

#################################################
### Removing residues outside the fibril core ###
#################################################

# Select and rename columns
metadata = metadata[['PDB ID', 'Residues Ordered']].rename(columns={'PDB ID': 'PDB', 'Residues Ordered': 'ordered_residues'})

# Add new columns for core start and end
metadata['core_start'] = None
metadata['core_end'] = None

# Check which PDBs are in metadata but not in df
missing_pdbs = set(metadata['PDB']) - set(df['PDB'])
print(f"Missing PDBs: {missing_pdbs}")

# Parse residue ranges
def parse_residue_range(residue_str):
    # Split on '-' or ','
    parts = re.split(r'[-,]', residue_str)
    parts = [p.strip() for p in parts if p.strip()]  # clean up empty strings
    return int(parts[0]), int(parts[-1])

metadata[['core_start', 'core_end']] = metadata['ordered_residues'].apply(
    lambda s: pd.Series(parse_residue_range(s))
)

# Adding Metadata to df
df = pd.merge(df, metadata, on = "PDB").copy()

# Removing residues outside the fibril core
df = df[(df["Pos"] >= df["core_start"]) & (df["Pos"] <= df["core_end"])]

df.drop(["ordered_residues", "core_start", "core_end"], axis=1, inplace=True)

#######################################
### Ading fibril and polymorph data ###
#######################################

# Renaming columns
df = df.rename(columns={"Mol": "chain"})

# Merging df and chain_fibril
df = pd.merge(df, chain_fibril, on = ["PDB", "chain"])

#####################################
### Removing Head and Tail Chains ###
#####################################

mask = ~df["PDB"].astype(str).add(df["chain"].astype(str)).isin(
    exterior_chains["PDB"].astype(str).add(exterior_chains["chain"].astype(str))
)
filtered_df = df[mask]

######################################################################
### Adding in pdb_id to handle PDBs with 2+ distinct amyloid folds ###
######################################################################

# Selecting columns of interest
pdb_id_data = com_and_fibril[["PDB", "pdb_id", "fibril"]].drop_duplicates()

# Adding pdb_id to filtered_df
filtered_df = pd.merge(pdb_id_data, filtered_df, on = ["PDB", "fibril"])

###################################################################
### Filtering based on Q-scores (see validation folder/scripts) ###
###################################################################

if remove_poorly_resolved == 1: 
    '''
    Removes entire PDBs with a mean Q-score below the threshold
    Additionally removes individual residues with low resolution from an otherwise well resolved PDB
    '''
    
    # Renaming column to match filtered_df
    high_resolution_residues = high_resolution_residues.rename(columns={"resno": "Pos"})
    
    # Merging data
    filtered_df = pd.merge(filtered_df, high_resolution_residues, on = ["pdb_id", "fibril", "Pos"])
    

if remove_poorly_resolved == 0:
    
    '''
    Removes entire PDBs with a mean Q-score below the threshold
    '''
    # Extracting high resolution pdb_ids
    high_resolution_PDBs = high_resolution_residues["pdb_id"].unique()
    
    # Keeping only well resolved pdb_ids
    filtered_df = filtered_df[filtered_df["pdb_id"].isin(high_resolution_PDBs)]

# Saving filtered_df
filtered_df.to_csv(os.path.join(thermo_path, "foldx_stability_filtered.csv"), index = False)


#############################################################################################################################
#############################################################################################################################
#############################################################################################################################

##############################################
### Plotting mean total deltaG per residue ###
##############################################

# Calculating the mean delta G per residue for each pdb_id
mean_per_res = (
    filtered_df.groupby(["pdb_id", "Pos"], as_index=False)
    .agg(mean_energy=("total", "mean"))
)

# Saving mean_per_res
mean_per_res.to_csv(os.path.join(thermo_path, "mean_deltaG_per_residue.csv"), index = False)

# Calculating the threshold for residues to be considered stabilising
mean_stability = mean_per_res['mean_energy'].mean(skipna=True)
SD_stability = mean_per_res['mean_energy'].std(skipna=True)
stabilising_threshold = mean_stability - SD_stability
destabilising_threshold = mean_stability + SD_stability

# Density Plot
p = (
    ggplot(mean_per_res, aes(x = "mean_energy", group = "pdb_id"))
    + geom_density()
    + geom_vline(xintercept = 0, colour = "black", size = 1, linetype = "dashed", alpha = 0.9)
    + scale_x_continuous(breaks=range(-100, 100, 1))
    + scale_y_continuous(expand = (0, 0, 0.05, 0))
    + labs(title="", 
           x=r'$\mathbf{\Delta G^{\circ}}$ per residue (kcal·mol$^{-1}$)', 
           y="Density")
    + theme_classic()
    + theme(
        panel_border = element_rect(color = "black", size = 1.5),
        panel_grid_major = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
        panel_grid_minor = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
        axis_title = element_text(size = 20, face = "bold"),
        axis_text = element_text(size = 14, colour = "black"),
        legend_position = "none"
    )
)

# Saving density plot
p.save(os.path.join(deltaG_analysis_path, "stability_density_plot.png"), height=6, width=8, dpi=300)

# Line Plot

# Creating a group variable based on pdb_id to include line breaks if residue numbers are not continuous
mean_per_res = mean_per_res.sort_values(['pdb_id', 'Pos'])
mean_per_res['segment'] = mean_per_res.groupby('pdb_id')['Pos'].diff().ne(1).cumsum()
mean_per_res['interaction'] = mean_per_res['pdb_id'].astype(str) + "_" + mean_per_res['segment'].astype(str)

# Plotting line plot
p = (
    ggplot(mean_per_res, aes(x='Pos'))
    + geom_line(aes(y='mean_energy', group='interaction'), size=0.75, alpha=0.25)
    + geom_hline(yintercept=0, colour="blue", size=0.75, linetype="solid", alpha=0.9)
    + scale_x_continuous(
        breaks=np.arange(0, 9999999, 1),
        labels=[str(i) if i % 10 == 0 else "" for i in range(0, 9999999, 1)],
        expand=(0.01, 0.01))
    + labs(
        y=r'$\mathbf{\Delta G^{\circ}\ per\ residue\ (kcal\cdot mol^{-1})}$',
        x="Residue Number"
    )
    + theme_classic()
    + theme(
        panel_border = element_rect(color = "black", size = 1.5),
        panel_grid_major = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
        panel_grid_minor = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
        axis_title = element_text(size = 20, face = "bold"),
        axis_text = element_text(size = 14, colour = "black"),
    )
)

# Saving line plot
p.save(os.path.join(deltaG_analysis_path, "stability_line_plot.png"), height=6, width=12, dpi=300)

######################################
### Plotting each PDB individually ###
######################################

pdb_names = filtered_df["pdb_id"].unique()

for pdb in pdb_names:
    x = mean_per_res[mean_per_res["pdb_id"] == pdb]
        
    p = (
        ggplot(x, aes(x='Pos'))
        + geom_line(aes(y='mean_energy', group='interaction'), size=2)
        + geom_hline(yintercept = 0, colour = "blue", size = 1, linetype = "dashed", alpha = 0.9)
        + geom_hline(yintercept = stabilising_threshold, colour = "red", size = 1, linetype = "dashed", alpha = 0.9)
        + scale_x_continuous(
            breaks=np.arange(0, 9999999, 1),
            labels=[str(i) if i % 10 == 0 else "" for i in range(0, 9999999, 1)],
            expand=(0.01, 0.01))
        + labs(
            title = f"{pdb}",
            y=r'$\mathbf{\Delta G^{\circ}\ per\ residue\ (kcal\cdot mol^{-1})}$',
            x="Residue Number"
        )
        + theme_classic()
        + theme(
            panel_border = element_rect(color = "black", size = 1.5),
            panel_grid_major = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
            panel_grid_minor = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
            plot_title = element_text(size = 20, face = "bold", hjust = 0.5),
            axis_title = element_text(size = 20, face = "bold"),
            axis_text = element_text(size = 14, colour = "black"),
        )
    )

    p.save(os.path.join(single_pdb_path, f"{pdb}_stability.png"), height=6, width=12, dpi=300)


#############################################################################
### Calculating the Pearson Correlation Coefficient Between PDB Stability ### 
############################################################################# 

# Select only the relevant columns
mean_per_res_selected = mean_per_res[['pdb_id', 'Pos', 'mean_energy']]

# Pivot without keeping Pos as a column
wide_df = mean_per_res_selected.pivot(index='Pos', columns='pdb_id', values='mean_energy').reset_index()

# Select columns for correlation calculation (excluding the first column, which is 'Pos')
columns_to_correlate = wide_df.iloc[:, 1:]

# Calculate correlation matrix, ignoring NA values
correlation_matrix = columns_to_correlate.corr()

# Saving pearson correlation
correlation_matrix.to_csv(os.path.join(pearson_path, "Pearson_Scores.csv"))

# Convert correlation matrix to long (tidy) format
correlation_long = correlation_matrix.reset_index().melt(
    id_vars='pdb_id', 
    var_name='pdb_id_2', 
    value_name='correlation'
)

### Plotting Pearson Correlation ###

# Heatmap
p = (
    ggplot(correlation_long, aes(x = "pdb_id", y = "pdb_id_2", fill = "correlation"))
    + geom_tile(size = 0.2, colour = "black")
    + scale_fill_gradientn(colors = ["blue", "white", "red"],
                       values=[0, 0.5, 1],
                       breaks = [-1,-0.5,0,0.5,1],
                       labels = [-1,-0.5,0,0.5,1],
                       limits = [-1,1]) 
    + labs(
        title = "Pearson Correlation between PDBs per residue FoldX stability",
        x = "",
        y = "",
        fill = "Pearson\nCorrelation"
    )
    + theme_minimal()
    + theme(
        plot_title = element_text(size = 18, face = "bold", colour = "black", ha = "center"),
        axis_text_y = element_text(size = 10, face = "bold", colour = "black"),
        axis_text_x = element_text(size = 10, face = "bold", colour = "black",
                                   angle = 90, vjust = 0.25, hjust = 1),
        axis_line = element_blank(),
        legend_title = element_text(size = 12, colour = "black", face = "bold"),
        legend_text = element_text(size = 10, colour = "black")
    )
)
# Saving Heatmap
p.save(os.path.join(pearson_path, "Pearson_heatmap.png"), height=11, width=13, dpi=300)

### Boxplot ###

# Create the comparison column
correlation_long['comparison'] = correlation_long['pdb_id'] + '-' + correlation_long['pdb_id_2']
# Remove self comparisons
correlation_long = correlation_long[correlation_long['pdb_id'] != correlation_long['pdb_id_2']]

# Creating a dummy collumn
correlation_long["dummy"] = "" 

p = (
    ggplot(correlation_long, aes(x = "dummy", y = "correlation"))
    + geom_violin(fill = "grey", size = 1, width = 0.5)
    + geom_boxplot(colour = "black", size = 1, width = 0.1, fill = None, outlier_size = 1)
    + labs(
        x = "",
        y = "Pearson Correlation",
    )
    + theme_classic()
    + theme(
        panel_border = element_rect(colour = "black", fill = None, size = 1.5),
        panel_grid_major = element_line(colour = "grey", size = 0.3, linetype = "dashed"),
        panel_grid_minor = element_line(colour = "grey", size = 0.3, linetype = "dashed"),
        axis_title = element_text(size = 20, face = "bold", colour = "black"),
        axis_text = element_text(size = 14, colour = "black")
    )
)
# Saving Heatmap
p.save(os.path.join(pearson_path, "Pearson_boxplot.png"), height=6, width=5, dpi=300)

### Saving Pearson correlation summary statistics ###
with open(os.path.join(pearson_path, "average_correlation.txt"), "w") as f:
    f.write("Average Pearson Correlation Score\n\n")
    f.write(f"Mean: {correlation_long['correlation'].mean(skipna=True)}\n")
    f.write(f"Median: {correlation_long['correlation'].median(skipna=True)}\n")
    f.write(f"SD: {correlation_long['correlation'].std(skipna=True)}\n")


####################################
### Plotting individual energies ### 
####################################

# Define the columns to average
#selected_columns <- c("total", "sideHbond", "energy_VdW", "electro", "energy_SolvP", "energy_SolvH", "energy_vdwclash")
selected_columns = filtered_df.columns[10:32] 

# Calculating the mean delta G per residue for each PDB
pdb_means = (
    filtered_df
    .groupby(['pdb_id', 'Pos'])[selected_columns]
    .mean()
    .reset_index()
    .rename(columns={col: f"mean_{col}" for col in selected_columns})
)

# Reshape from wide to long format
long_df = (
    pdb_means
    .melt(id_vars=['pdb_id', 'Pos'],
          value_vars=[col for col in pdb_means.columns if col.startswith('mean_')],
          var_name='variable',
          value_name='value')
)

# Remove rows with NA values
long_df = long_df.dropna(subset=['value'])

# Remove the 'mean_' prefix from variable names
long_df['variable'] = long_df['variable'].str.replace('mean_', '', regex=False)

# Sort data
long_df = long_df.sort_values(by=['variable', 'pdb_id', 'Pos']).copy()

# Ensure the variable column is a categorical with the same order as selected_columns
long_df['variable'] = pd.Categorical(long_df['variable'], categories=selected_columns, ordered=True)

# Create a grouping column
# i.e. break the line when Pos is not consecutive
long_df['interaction'] = (
    long_df.groupby('pdb_id')['Pos']
    .transform(lambda x: np.cumsum(np.r_[0, np.diff(x) != 1]))
)

# Combine pdb_id and interaction to form a unique group identifier
long_df['line_group'] = long_df['pdb_id'].astype(str) + '_' + long_df['interaction'].astype(str)

# Create the plot
p = (
    ggplot(long_df, aes(x='Pos', y='value'))
    + geom_line(aes(group = 'line_group'), size = 1, alpha = 0.25)
    + geom_hline(yintercept = 0, size = 0.5, linetype = 'dashed')
    + facet_wrap('~variable')
    + labs(
        y = r'$\mathbf{\Delta G^{\circ}}$ per residue (kcal·mol$^{-1}$)',
        x = "Position"  
    )
    + theme_minimal()
    + theme(
        axis_line = element_blank(),
        strip_text = element_text(size = 14, face = "bold"),
        axis_title = element_text(size = 22, face = "bold", colour = "black"),
        axis_text = element_text(size = 14, colour = "black"),
        legend_position='none'
    )
)

# Saving plot
p.save(os.path.join(deltaG_analysis_path, "all_energies_separate_pdb.png"), height=12, width=15, dpi=300)

###################################
### Characterizing Residue Type ### 
###################################

# Replace "H1S" and "H2S" with "HIS"
filtered_df.loc[filtered_df['Code'] == 'H1S', 'Code'] = 'HIS'
filtered_df.loc[filtered_df['Code'] == 'H2S', 'Code'] = 'HIS'

# Count how many PDBs each residue position is found in
all_pos_counts = filtered_df.groupby('Pos').size().reset_index(name='n')

# Count how many PDBs each residue position is stabilising / destabilising
stabilising_df = filtered_df[filtered_df['total'] <= stabilising_threshold]
destabilising_df = filtered_df[filtered_df['total'] >= destabilising_threshold]

# Create amino acid property DataFrame
amino_acids = pd.DataFrame({
    "Code": ["ALA", "ARG", "ASN", "ASP", 
             "CYS", "GLN", "GLU", "GLY", 
             "HIS", "ILE", "LEU", "LYS", 
             "MET", "PHE", "PRO", "SER", 
             "THR", "TRP", "TYR", "VAL"],
    "property": ["Hydrophobic", "Positive", "Polar", "Negative", 
                 "Polar", "Polar", "Negative", "Polar", 
                 "Positive", "Hydrophobic", "Hydrophobic", "Positive", 
                 "Hydrophobic", "Aromatic", "Hydrophobic", "Polar", 
                 "Polar", "Aromatic", "Aromatic", "Hydrophobic"]
})

# Counting the total of each residue type for all PDBs
all_residue_counts = filtered_df.groupby('Code').size().reset_index(name='total_count')

# Counting the residue type from stabilising residues for all PDBs
stable_residue_counts = stabilising_df.groupby('Code').size().reset_index(name='stabilising_count')

# Counting the residue type from destabilising residues for all PDBs
destable_residue_counts = destabilising_df.groupby('Code').size().reset_index(name='destabilising_count')

# Merge the counts together
residue_count_df = pd.merge(stable_residue_counts, destable_residue_counts, on='Code', how='outer')
residue_count_df = pd.merge(residue_count_df, all_residue_counts, on='Code', how='left')

# Normalising
residue_count_df['norm_stable_count'] = (residue_count_df['stabilising_count'] / residue_count_df['total_count']) * 100
residue_count_df['norm_destable_count'] = (residue_count_df['destabilising_count'] / residue_count_df['total_count']) * 100

# Add residue properties for colouring
residue_count_df = pd.merge(residue_count_df, amino_acids, on='Code', how='left')

# Replacing NaN values with 0
residue_count_df = residue_count_df.fillna(0)

# Save residue counts to CSV
residue_count_df.to_csv(os.path.join("Output", "thermodynamics", "residue_analysis", "residue_counts.csv"), index=False)

### Plotting Stable and Destable Residue Counts ###
def plot_residue_counts(df, value_col, y_label, output_file, fill_colors=None, height=8, width=7, dpi=300):
    """
    df: DataFrame containing 'Code', value_col, and 'property'
    value_col: column to plot on y-axis (e.g., 'norm_stable_count')
    y_label: y-axis label
    output_file: path to save plot
    fill_colors: dict mapping property -> color
    """
    # Sort by value_col and reset Code categorical
    df_sorted = df.sort_values(by=value_col, ascending=False).copy()
    df_sorted['Code'] = pd.Categorical(
        df_sorted['Code'],
        categories=df_sorted['Code'].tolist(),
        ordered=True
    )

    # Create plot
    p = (
        ggplot(df_sorted, aes(x="Code", y=value_col, fill="property"))
        + geom_col()
        + labs(y=y_label, x="")
        + scale_x_discrete(expand=(0.01, 0.01))
        + scale_y_continuous(expand=(0, 0.05))
    )
    
    # Apply fill colours if provided
    if fill_colors:
        p += scale_fill_manual(values=fill_colors)

    # Save plot
    p.save(output_file, height=height, width=width, dpi=dpi)
    print(f"Saved plot to {output_file}")


# Define colours for residue properties
fill_colors = {
    "Aromatic": "sienna",
    "Hydrophobic": "goldenrod",
    "Negative": "dodgerblue",
    "Polar": "yellowgreen",
    "Positive": "crimson"
}

# Plot stable counts
plot_residue_counts(
    residue_count_df,
    value_col="norm_stable_count",
    y_label="Stabilising Occurences / Total Occurences (%)",
    output_file=os.path.join(residue_analysis_path, "stabilising_residue_type_normalised.png"),
    fill_colors=fill_colors
)

# Plot destabilising counts
plot_residue_counts(
    residue_count_df,
    value_col="norm_destable_count",
    y_label="Destabilising Occurences / Total Occurences (%)",
    output_file=os.path.join(residue_analysis_path, "destabilising_residue_type_normalised.png"),
    fill_colors=fill_colors
)

### Plotting residue type and mean dG ###

res_type_df = (
    filtered_df
    .groupby(['Pos', "Code"])["total"]
    .mean()
    .reset_index(name='mean_energy')
    )

res_type_df = pd.merge(res_type_df, amino_acids, on = "Code")

# Line Plot
p = (
    ggplot(res_type_df, aes(x = "Pos", y = "mean_energy"))
    + geom_line(size = 0.75)
    + geom_point(aes(colour = "property"), size = 1.5)
    + labs(
        y = r'Mean $\mathbf{\Delta G^{\circ}}$ per residue (kcal·mol$^{-1}$)',
        x = "Residue Number",
        colour = "Property"
    )
    + theme_classic()
    + theme(
        panel_border = element_rect(colour = "black", size = 1),
        panel_grid_major = element_line(colour = "grey", size = 0.25, linetype = "dashed"),
        panel_grid_minor = element_line(colour = "grey", size = 0.25, linetype = "dashed"),
        axis_title = element_text(size = 12, face = "bold"),
        axis_text = element_text(size = 10, colour = "black"),
        legend_title = element_text(size = 12, colour = "black", face = "bold"),
        legend_text = element_text(size = 10, colour = "black")
    )
)

# Saving Line plot
p.save(os.path.join(residue_analysis_path, "mean_dG_per_residue_all_PDBs.png"), height=4, width=10, dpi=300)

####################################
### Sum deltaG per PDB analysis ###
####################################

sum_per_pdb = (
    filtered_df
    .groupby(["pdb_id", "Pos"])["total"]
    .mean()
    .reset_index(name = "mean_energy")
    .groupby("pdb_id")["mean_energy"]
    .sum()
    .reset_index(name = "sum_mean_energy")
)

### By RMSD cluster group ###

# Selecting columns of interesst from cluster_groups
cluster_groups = cluster_groups[["pdb_id", "group"]]

# Remove all instances of duplicated rows based on the pdb_id column
cluster_groups = cluster_groups[cluster_groups['pdb_id'].map(cluster_groups['pdb_id'].value_counts() == 1)]

# Creating rmsd_cluster_df
rmsd_cluster_df = pd.merge(sum_per_pdb, cluster_groups, on = "pdb_id")
rmsd_cluster_df["group"] = rmsd_cluster_df["group"].astype(str)
rmsd_cluster_df = rmsd_cluster_df.sort_values(by='group').reset_index(drop=True)

# Saving rmsd_cluster_df
rmsd_cluster_df.to_csv(os.path.join(deltaG_comp_path, "sum_deltaG_per_PDB.csv"), index = False)

# Plotting sum_mean_energy for each RMSD cluster group
p = (
    ggplot(rmsd_cluster_df, aes(x = "group", y = "sum_mean_energy")) 
    + geom_boxplot(size = 0.75, colour = "black")
    + geom_point(size = 2.5)
    + labs(
        y = r'$\mathbf{\Delta G^{\circ}}$ per PDB (kcal·mol$^{-1}$)',
        x = "RMSD Cluster Group"
    )
    + theme_classic()
    + theme(panel_border = element_rect(colour = "black", size = 1.5),
        panel_grid_major = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
        panel_grid_minor = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
        axis_text = element_text(size = 16, colour = "black", face = "bold"),
        axis_title = element_text(size = 18, face = "bold"))
)

p.save(os.path.join(deltaG_comp_path, "sum_deltaG_per_PDB_by_RMSD_cluster.png"), height = 6, width = 8, dpi = 300)

### Comparing by fibril source ###

# Reimporting metadata
metadata = pd.read_csv(os.path.join("Output", "selected_pdbs_metadata.csv"))

# Rename column
metadata = metadata.rename(columns={"PDB ID": "PDB"})

# Select the column for pattern matching
fibril_origins = metadata['Fibril Origins']  

# Create boolean masks for each feature
cond_extracted = fibril_origins.str.contains(r"extracted|patient|case|atrophy", case=False, na=False)  
# True if the text mentions 'extracted', 'patient', 'case', or 'atrophy' (case-insensitive).  
# na=False ensures that missing values (NaN) are treated as False.

cond_seeded = fibril_origins.str.contains(r"seed|seeded", case=False, na=False)  
# True if the text mentions 'seed' or 'seeded', also case-insensitive.

# Define conditions
conditions = [
    cond_extracted & cond_seeded,        # both ex vivo terms AND seeded → "Seeded From Ex Vivo"
    cond_extracted,                      # only ex vivo terms → "Ex Vivo"
    ~cond_extracted & cond_seeded,       # only seeded but not ex vivo → "Seeded From In Vitro"
    ~cond_extracted                      # neither extracted nor ex vivo terms → "In Vitro"
]

# Corresponding labels for each condition
choices = [
    "Seeded From Ex Vivo",
    "Ex Vivo",
    "Seeded From In Vitro",
    "In Vitro"
]

# 5. Apply conditions to create a new column
metadata['condition'] = np.select(conditions, choices, default="other")  
# np.select goes through the list of conditions in order and assigns the corresponding choice.  
# If none of the conditions match, it defaults to "other".

# Selecting columns
metadata = metadata[["PDB", "condition"]]

# Creating PDB from pdb_id to merge sum_per_pdb and metadata
# If pdb_id contains "_", take the part before it; otherwise keep pdb_id
sum_per_pdb['PDB'] = sum_per_pdb['pdb_id'].str.split('_').str[0]

# Merging sum_per_pdb and metadata
sum_energy_by_condition = pd.merge(sum_per_pdb, metadata, on = "PDB")

# Saving sum_energy_by_condition
sum_energy_by_condition.to_csv(os.path.join(deltaG_comp_path, "sum_deltaG_per_PDB_by_fibril_condition.csv"), index = False)

# Plotting sum_mean_energy for each fibril formation condition
p = (
    ggplot(sum_energy_by_condition, aes(x = "condition", y = "sum_mean_energy"))
    + geom_boxplot(size = 0.75, colour = "black")
    + geom_point(size = 2.5)
    + labs(
        y = r'$\mathbf{\Delta G^{\circ}}$ per PDB (kcal·mol$^{-1}$)',
        x = "Fibril Condition"
        )
    + scale_x_discrete(labels = {"Ex vivo": "Ex Vivo",
                                 "In Vitro": "In Vitro",
                                 "Seeded From In Vitro": "Seeded From\nIn Vitro",
                                 "Seeded From Ex Vivo": "Seeded From\nEx Vivo"}
                                 )
    + theme_classic()
    + theme(
        panel_border = element_rect(colour = "black", size = 1.5),
        panel_grid_major = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
        panel_grid_minor = element_line(colour = "grey", size = 0.5, linetype = "dashed"),
        axis_text = element_text(size = 14, colour = "black", face = "bold"),
        axis_title = element_text(size = 18, face = "bold")
        )
)
# Saving plot
p.save(os.path.join(deltaG_comp_path, "sum_deltaG_per_PDB_by_fibril_condition.png"), height = 6, width = 8, dpi = 300)

###########################
### Stability Bar Chart ###
###########################

# Calculate mean total energy per residue position
tmp = (
    filtered_df
    .groupby('Pos', as_index=False)['total']
    .mean()
    .rename(columns={'total': 'mean_energy'})
)

# Create colour_category based on thresholds
tmp['colour_category'] = np.select(
    [
        tmp['mean_energy'] > 0,
        (tmp['mean_energy'] <= 0) & (tmp['mean_energy'] > stabilising_threshold),
        tmp['mean_energy'] <= stabilising_threshold
    ],
    [
        'above_0',
        'below_0_above_threshold',
        'below_threshold'
    ],
    default='other'
)

# Plotting
p = (
    ggplot(tmp, aes(x = "Pos", y = "mean_energy", fill = "colour_category")) 
    + geom_col() 
    + geom_hline(yintercept = 0, colour = "black", size = 0.5)
    + geom_hline(yintercept = stabilising_threshold, colour = "#f27272", size = 1, linetype = "dashed") 
    + scale_x_continuous(breaks = np.arange(0, 9999999, 10), expand = (0,0)) 
    + scale_fill_manual(values = {
        "above_0": "#6BAF5F",
        "below_0_above_threshold": "#fabb1b",
        "below_threshold": "#f27272"
    }
    ) 
    + labs(
        y = r'Mean $\mathbf{\Delta G^{\circ}}$ per residue (kcal·mol$^{-1}$)',
        x = "Residue Number"
        ) 
    + theme_classic() 
    + theme(axis_text = element_text(size = 12, colour = "black", face = "bold"),
            axis_title = element_text(size = 14, face = "bold"),
            legend_position = "none"
            ) 
)
# Saving plot
p.save(os.path.join(deltaG_analysis_path, "bar_chart.png"), height = 5, width = 8, dpi = 300)