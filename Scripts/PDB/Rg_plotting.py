import pandas as pd
import numpy as np
import os
from plotnine import *

# Import data
df = pd.read_csv(os.path.join("Output", "PDBs", "COM_and_fibril.csv"), sep=",")

print(df.head())

# Mean Rg for each fibril per PDB
fibril_mean_df = (
    df.groupby(["PDB", "fibril"], as_index=False)
      .agg(fibril_mean_rg=("Rg", "mean"))
)

# Calculate % difference from fibril 1
def diff_from_fibril1(group):
    # get fibril 1’s mean Rg
    fibril1_rg = group.loc[group["fibril"] == 1, "fibril_mean_rg"]
    if fibril1_rg.empty:
        group["diff_from_fibril1"] = np.nan
    else:
        fibril1_rg = fibril1_rg.iloc[0]
        group["diff_from_fibril1"] = (
            abs(group["fibril_mean_rg"] - fibril1_rg) / group["fibril_mean_rg"] * 100
        )
    return group

diff_df = fibril_mean_df.groupby("PDB", group_keys=False).apply(diff_from_fibril1)

# Order by difference magnitude
rg_diff_order = (
    diff_df.groupby("PDB", as_index=False)
           .agg(max_diff=("diff_from_fibril1", "max"))
           .sort_values("max_diff", ascending=True)
)

# Ensure PDB factor order matches rg_diff_order
diff_df["PDB"] = pd.Categorical(
    diff_df["PDB"], categories=rg_diff_order["PDB"], ordered=True
)

# Convert fibril to string for discrete colours
diff_df["fibril"] = diff_df["fibril"].astype(int).astype(str)

# --- PLOT ---

# Count distinct PDB values
num_pdbs = diff_df["PDB"].nunique()

# Scaling plot width with number of PDBs
width_per_pdb = 0.2  
dynamic_width = min(max(8, num_pdbs * width_per_pdb), 25)

p = (
    ggplot(diff_df, aes(x="PDB", y="diff_from_fibril1", colour="fibril"))
    + geom_point(size=3, alpha=0.75, stroke=1)
    + geom_hline(yintercept=5, colour="black", linetype="dashed", size=1)
    + labs(
        y="% difference in the mean Rg\nfor each fibril compared to fibril 1",
        colour="Fibril Number",
    )
    + scale_colour_discrete(limits = sorted(diff_df["fibril"].unique(), key = int))
    + scale_y_continuous(breaks = np.arange(0, diff_df["diff_from_fibril1"].max() + 5, 5))
    + theme_bw()
    + theme(
        axis_title=element_text(size=16, colour="black", weight="bold"),
        axis_text_y=element_text(size=11, colour="black"),
        axis_text_x=element_text(size = min(max(6, 11 - (num_pdbs*0.1)), 11), 
                                 angle=90, hjust=0.5, vjust=0.5, colour="black"),
        legend_title=element_text(size=16, colour="black", weight="bold"),
        legend_text=element_text(size=14, colour="black"),
    )
)

# --- SAVE PLOT ---
output_path = os.path.join("Output", "PDBs", "rg_difference_plot.png")
p.save(output_path, height=10, width=dynamic_width, dpi=300)