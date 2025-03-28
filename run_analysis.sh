#!/bin/bash

# Setting Run Parameters
scrape=0                        # 0 - Don't Web Scrape              | 1 - Web Scrape Amyloid Atlas
PDB=0                           # 0 - Don't Analyse PDBs            | 1 - Analyse PDBs
align=0                         # 0 - Do not align                  | 1 - Run alignment
validation=0                    # 0 - Do not validate               | 1 - Run validation
RMSD=0                          # 0 - Do not calculate              | 1 - Run RMSD
thermodynamics=0                # 0 - Do not run thermodynamics     | 1 - Run thermodynamic analysis
stable_regions=0                # 0 - Do not analyse stable regions | 1 - Run stable region analysis
b_sheet=0                       # 0 - Do not analyse B-Sheet        | 1 - Run B-Sheet

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

if [ $scrape -eq  1 ]; then
    # Web scraping Amyloid Atlas
    python Scripts/web_scraping/amyloid_atlas_scraper.py

    # Plotting the residues ordered in the fibril core
    Rscript Scripts/web_scraping/plot_ordered_residues.R
fi

if [ $PDB -eq 1 ]; then
    # Fetch pdbs to get .cif files
    python Scripts/PDB_analysis/fetch_pdb_isolate_chains.py

    # Create info df for each PDB for use in subsequent analysis
    python Scripts/PDB_analysis/extract_fibril_info.py

    # Asymmetric Unit Figure
    python Scripts/PDB_analysis/asymmetric_unit_png_generator.py
    python Scripts/PDB_analysis/asymmetric_unit_figure_maker.py
fi

if [ $align -eq 1 ]; then
    # Align unique chains
    python Scripts/PDB_analysis/unique_chain_alignment.py
fi

if [ $validation -eq 1 ]; then
    # Get Q-Scores
    python Scripts/Q_scores/Q_score_scraper.py

    # Selecting good resolution structures
    Rscript Scripts/Q_scores/validating_structures.R

fi

if [ $RMSD -eq 1 ]; then
    # Converting .cif alignment to a .csv file
    python Scripts/RMSD/create_df.py

    # Performing RMSD analysis
    Rscript Scripts/RMSD/RMSD_analysis.R <<EOF
plot=TRUE 
custom_cut_height=6.7
EOF
# Plot individual comparisons for each structure (TRUE/FALSE)
# On first run, set custom_cut_height=0 and it will use the mean euclidean distance for the cut height

# Tau 3R+4R and 4R Variants cut height = 5.5
# asyn cut height = 6.7 (old version = 8.3)

fi

if [ $thermodynamics -eq 1 ]; then
   # Scrape EMDB to get fibril twist and rise values
   python Scripts/thermodynamics/EMDB_scraper.py

   # Extend the asymetric units to a layer depth of 10 chains
   python Scripts/thermodynamics/extend_fibril_layers.py

   # Calculate the FoldX per residue stability
   python Scripts/thermodynamics/foldx_analysis.py

   # Plotting thermodynamic results
   Rscript Scripts/thermodynamics/thermodynamics_plotting.R

fi

if [ $stable_regions -eq 1 ]; then

    # Defining Stable Regions
    Rscript Scripts/stable_regions/defining_stable_regions.R

    # Calculate the Stable Region Distances
    Rscript Scripts/stable_regions/stable_region_distances.R

    # Colouring asymetric units by stable regions
    python Scripts/stable_regions/stable_region_colouring.py

    # Creating a figure showing each structure with coloured stable regions and grouped by RMSD
    python Scripts/stable_regions/png_maker.py
    python Scripts/stable_regions/cluster_group_and_stable_regions_figure.py

fi

if [ $b_sheet -eq  1 ]; then
    # Analyse B-sheets
    Rscript Scripts/beta_sheet_analysis/beta_sheet_analysis.R
fi