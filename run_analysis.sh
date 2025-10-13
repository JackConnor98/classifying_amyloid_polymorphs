#!/bin/bash

# Setting Run Parameters
scrape=0                        # 0 - Don't Web Scrape              | 1 - Web Scrape Amyloid Atlas
PDB=0                           # 0 - Don't Analyse PDBs            | 1 - Analyse PDBs
validation=1                    # 0 - Do not validate               | 1 - Run validation
RMSD=0                          # 0 - Do not calculate              | 1 - Run RMSD
thermodynamics=0                # 0 - Do not run thermodynamics     | 1 - Run thermodynamic analysis
stable_regions=0                # 0 - Do not analyse stable regions | 1 - Run stable region analysis
beta_sheet=0                    # 0 - Do not analyse Beta-Sheet     | 1 - Run Beta-Sheet
PNG=0                           # 0 - Do not generate PNGs          | 1 - Generate PNGs

#########################
### Optional settings ###
#########################

# Would you like to add a penalty to non-overlapping residues in the RMSD calculation?
penalty=1                       # 0 - No penalty | 1 - Mean + 1SD | 2 - Mean + 2SD | 3 - Mean + 3SD etc...

# Set a custom cut height for the dendrogram in the RMSD analysis (Reccomened to use 0 for the first run)
custom_cut_height=0           # 0 - Use the median euclidean distance as the cut height | Any other number will be used as the cut height

# Would you like to exclude poorly resolved single residues from thermodynamic analysis? 
# 0 - No | 1 - Yes
remove_poorly_resolved=1

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

if [ $scrape -eq  1 ]; then
    # Web scraping Amyloid Atlas
    python Scripts/scrape/amyloid_atlas_scraper.py

    # Plotting the residues ordered in the fibril core
    Rscript Scripts/scrape/plot_ordered_residues.py
fi

if [ $PDB -eq 1 ]; then
    # Fetch pdbs to get .cif files
    python Scripts/PDB/fetch_pdb_isolate_chains.py

    # Calculating the maximum Rg within a PDB
    python Scripts/PDB/Rg_plotting.py

fi

if [ $validation -eq 1 ]; then
    # Get Q-Scores
    python Scripts/validation/Q_score_scraper.py

    # Selecting good resolution structures
    python Scripts/validation/validating_structures.py

fi

if [ $RMSD -eq 1 ]; then

    # Align unique chains
    python Scripts/RMSD/unique_chain_alignment.py $penalty

    # Performing RMSD analysis
    Rscript Scripts/RMSD/RMSD_analysis.R <<EOF
custom_cut_height=$custom_cut_height
EOF

fi


if [ $thermodynamics -eq 1 ]; then

   # Getting fibril twist and rise
   python Scripts/thermodynamics/EMDB_scraper.py

   # Extend the asymetric units to a layer depth of 10 chains
   python Scripts/thermodynamics/extend_fibril_layers.py

   # Calculate the FoldX per residue stability
   python Scripts/thermodynamics/foldx_analysis.py

   # Plotting thermodynamic results
   Rscript Scripts/thermodynamics/thermodynamics_plotting.R <<EOF
remove_poorly_resolved=$remove_poorly_resolved
EOF
fi

if [ $stable_regions -eq 1 ]; then

    # Defining Stable Regions
    Rscript Scripts/stable_regions/defining_stable_regions.R

    # Calculate the Stable Region Distances
    Rscript Scripts/stable_regions/stable_region_distances.R

fi

if [ $beta_sheet -eq  1 ]; then
    # Analyse B-sheets
    Rscript Scripts/beta_sheet/beta_sheet_analysis.R
fi

if [ $PNG -eq 1 ]; then

    # Asymmetric Unit Figure
    python Scripts/PNG/asymmetric_unit_png_generator.py
    python Scripts/PNG/asymmetric_unit_figure_maker.py

    # Colouring asymetric units by stable regions
    python Scripts/PNG/stable_region_colouring.py

    # Creating a figure showing each structure with coloured stable regions and grouped by RMSD
    python Scripts/PNG/stable_region_png_maker.py
    python Scripts/PNG/cluster_group_and_stable_regions_figure.py

fi