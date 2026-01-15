#!/bin/bash

# TODO
# Calculate intersheet spacing
# Give option: manually calculate twist and rise only when EMDB data is not available
# Fix plot size scaling to handle large numbers of PDBs better
# Test Local PDB functionality fully

# Setting Run Parameters
scrape=0                        # 0 - Don't Web Scrape              | 1 - Web Scrape Amyloid Atlas
PDB=0                           # 0 - Don't Analyse PDBs            | 1 - Analyse PDBs
validation=0                    # 0 - Do not validate               | 1 - Run validation
RMSD=0                          # 0 - Do not calculate              | 1 - Run RMSD
thermodynamics=1                # 0 - Do not run thermodynamics     | 1 - Run thermodynamic analysis
stable_regions=0                # 0 - Do not analyse stable regions | 1 - Run stable region analysis
beta_strand=0                   # 0 - Do not analyse Beta-Sheet     | 1 - Run Beta-Sheet
PNG=0                           # 0 - Do not generate PNGs          | 1 - Generate PNGs

#########################
### Optional settings ###
#########################

# Would you like to use Local PDBs?
use_local=0                    # 0 = No | 1 = Yes

# Web scraping: would you like to use the GUI (1) or command line (0) version
scrape_version=1

# Specify the Q-score threshold [0-1] (default is set to the mean - SD Q-score accross all PDBs)
q_score_threshold="automatic"

# Would you like to add a penalty to non-overlapping residues in the RMSD calculation?
penalty=0                     # 0 - No penalty | 1 - Mean + 1SD | 2 - Mean + 2SD | 3 - Mean + 3SD etc...

# Set a custom cut height for the dendrogram in the RMSD analysis (Reccomened to use 0 for the first run)
custom_cut_height=0          # 0 - Use the median euclidean distance as the cut height | Any other number will be used as the cut height

# Would you like to use publised twist and rise values or calulate them locally?
twist_rise_source=0     # 0 - Use EMDB values where available | 1 - Calculate locally

# Specify the number of layers to extend the fibril to (default is 10)
num_layers=10

# Would you like to exclude poorly resolved single residues from thermodynamic analysis? 
# 0 - No | 1 - Yes
remove_poorly_resolved=0

# Specify FoldX parameters - if not specified it will use the default [pH = 7.4, temp = 298K, ionstrength = 0.150]
# The specified parameters will be applied to all structures
# I plan to update this in the future to allow different parameters for each PDB to better capture the fibril formation conditions
ph=NA              
temp=NA # Temperature in Kelvin                
ionstrength=NA          

# Specify a number of parallel process to be used for side chain relaxation (Bottleneck of the analysis is FoldX's RepairPDB)
num_processes=4

# Set window size for rolling average of FoldX data (3 is reccomended for most proteins)
window=3
min_stable_region_size=0

# Specify an intersheet distance threshold for stable region contacts (If not provided it will use the default of 10.8A)
distance_threshold=NA

# Specify the minimum number of residues required to be considered a B-strand (default = 4)
min_length=NA

# Specify how you want the PNGs of fibril units to be coloured
# 0 = Single Colour
# 1 = Protein Domain
# 2 = Protofilament
# 3 = Amyloid Fold
png_colouring=0
# If you selected Single Colour (0) Protofilament (2) or Amyloid Fold (3) you can enter a list of colours to use otherwise it will use a default palette
png_palette="tv_blue, tv_red, tv_green, tv_orange"

#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################
#############################################################################################################################################

if [ $scrape -eq  1 ]; then

    # Web scraping Amyloid Atlas
    if [ $scrape_version -eq 1 ]; then
        python Scripts/scrape/amyloid_atlas_scraper_gui.py
    else
        python Scripts/scrape/amyloid_atlas_scraper.py
    fi
    # Plotting the residues ordered in the fibril core
    python Scripts/scrape/plot_ordered_residues.py $use_local
fi

if [ $PDB -eq 1 ]; then
    # Fetch pdbs to get .cif files
    python Scripts/PDB/fetch_pdb_isolate_chains.py $use_local

    # Calculating the maximum Rg within a PDB
    python Scripts/PDB/Rg_plotting.py

    # Sequence aligning PDBs
    python Scripts/PDB/sequence_alignment.py

fi

if [ $validation -eq 1 ]; then
    # Get Q-Scores
    python Scripts/validation/Q_score_scraper.py

    # Selecting good resolution structures
    python Scripts/validation/validating_structures.py $q_score_threshold $use_local

fi

if [ $RMSD -eq 1 ]; then

    # Align unique chains
    python Scripts/RMSD/unique_chain_alignment.py $penalty

    # Performing RMSD analysis
    python Scripts/RMSD/RMSD_analysis.py $custom_cut_height

fi

if [ $thermodynamics -eq 1 ]; then

   # Getting fibril twist and rise
   python Scripts/thermodynamics/EMDB_scraper.py 

   # Calculating twist and rise locally
   python Scripts/thermodynamics/calculate_twist_rise.py

   # Extend the asymetric units to a layer depth of 10 chains
   python Scripts/thermodynamics/extend_fibril_layers.py $twist_rise_source $num_layers

   # Calculate the FoldX per residue stability
   python Scripts/thermodynamics/foldx_analysis_multithread.py $ph $temp $ionstrength $num_processes

   # Plotting thermodynamic results
   python Scripts/thermodynamics/thermodynamics_plotting.py $remove_poorly_resolved

fi

if [ $stable_regions -eq 1 ]; then

    # Defining Stable Regions
    python Scripts/stable_regions/defining_stable_regions.py $window $min_stable_region_size

    # Calculate the Stable Region Distances
    python Scripts/stable_regions/stable_region_distances.py $distance_threshold
fi

if [ $beta_strand -eq  1 ]; then
    # Analyse B-sheets
    python Scripts/beta_strand/beta_strand_analysis.py $min_length
fi

if [ $PNG -eq 1 ]; then

    # Generate PNGs
    python Scripts/PNG/generate_pngs.py $png_colouring $png_palette
    
    # Scrape the polarity maps from Amyloid Atlas
    python Scripts/PNG/polarity_map_scraper.py

fi