import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import pandas as pd
from Bio.PDB import PDBParser
import requests
from io import StringIO
import numpy as np
import pymol
from pymol import cmd

'''

### TODO ###
- Add contact threshold to legend
- add colouring by stable regions
- Add the ability to specify residues to colour
- Make seperate chains more obvious (transparent or different colour backbone)
- Add N and C term labels
- make compare function to show two PDBs and only show unique contacts to each fold (or shared contacts)

'''

### --- Global Variables

# Mapping of 3-letter amino acid codes to 1-letter codes
AA_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
    'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
    'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
    'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
    'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
}

AA_LOOKUP = pd.DataFrame({
    "Code": ["A", "R", "N", "D", 
            "C", "Q", "E", "G", 
            "H", "I", "L", "K", 
            "M", "F", "P", "S", 
            "T", "W", "Y", "V"],
    "property": ["Hydrophobic", "Positive", "Polar", "Negative", 
                "Polar", "Polar", "Negative", "Polar", 
                "Positive", "Hydrophobic", "Hydrophobic", "Positive", 
                "Hydrophobic", "Aromatic", "Hydrophobic", "Polar", 
                "Polar", "Aromatic", "Aromatic", "Hydrophobic"],
    "colour": ["goldenrod", "crimson", "yellowgreen", "dodgerblue",
            "yellowgreen" ,"yellowgreen", "dodgerblue", "yellowgreen",
            "crimson", "goldenrod", "goldenrod", "crimson",
            "goldenrod", "sienna", "goldenrod", "yellowgreen",
            "yellowgreen", "sienna", "sienna", "goldenrod"]
})

PROP_TO_COLOUR = dict(zip(AA_LOOKUP["property"], AA_LOOKUP["colour"]))
CODE_TO_PROP = dict(zip(AA_LOOKUP["Code"], AA_LOOKUP["property"]))
CODE_TO_COLOUR = dict(zip(AA_LOOKUP["Code"], AA_LOOKUP["colour"]))

### --- Functions ---

######################################################################
### 1) Reading in local or remote PDB files and extracting CA atom ###
######################################################################

# Fetch PDB structure from RCSB given a PDB code
def fetch_pdb_structure(pdb_code):

    parser = PDBParser(QUIET = True)

    #url = f"https://files.rcsb.org/download/{pdb_code}.pdb"
    url = f"https://www.rcsb.org/pdb/download/downloadFile.do?fileFormat=pdb&compression=NO&structureId={pdb_code}"
    
    response = requests.get(url)
    
    structure = None
    
    if response.status_code != 200:
        
        structure = fetch_structure_pymol(pdb_code)
        
        if not structure:
            raise ValueError(f"Could not fetch PDB: {pdb_code}")
        
    if not structure:
        
        pdb_data = StringIO(response.text)
        
        structure = parser.get_structure(pdb_code, pdb_data)
    
    return structure

# Read PDB structure from a local file
def read_local_pdb(pdb_file):
    parser = PDBParser(QUIET=True)

    try:

        pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
        pdb_name = pdb_name.replace("_asym_unit", "")
        
        structure = parser.get_structure(pdb_name, pdb_file)
        
        return structure
          
    except FileNotFoundError:
        print(f"PDB not found at {pdb_file}")
        return None

    except Exception as e:
        print(f"Error reading {pdb_file}: {e}")
        return None

def launch_pymol():
    # Set PyMOL to run without GUI
    pymol.pymol_argv = ['pymol', '-qc']
    pymol.finish_launching()
    
def close_pymol():
    # Close PyMol
    cmd.quit()
    
def remove_cif():

    # Deleting the .cif files that automatically save
        for filename in os.listdir():
            if filename.endswith((".cif", ".pdb1")):
                file_path = os.path.join(filename)
                os.remove(file_path)
    
# Failsafe function to fetch PDB using python if fetching from website fails
def fetch_structure_pymol(pdb_code):

    # Reinitialize PyMol
    cmd.reinitialize()
    
    # Fetch the pdb
    cmd.fetch(pdb_code, type="pdb1")

    # check whether PyMOL created the object
    if not cmd.get_object_list(pdb_code):
        print(f"Biological assembly failed for {pdb_code}, retrying...")
        cmd.reinitialize()
        cmd.fetch(pdb_code)
        check_states = False
    else:
        check_states = True


    # Handling PDBs with multiple states (e.g. 8azs)
    if check_states:
        # Count how many states (models) are present
        n_states = cmd.count_states(pdb_code)
        if n_states > 1:
            print(f"{pdb_code} has {n_states} state(s).")
            print(f"→ Splitting {n_states} states and recombining into one object...")
            # Split all states into separate objects
            cmd.split_states(pdb_code)
            # Delete the original object
            cmd.delete(pdb_code)
            # Create an empty object to hold the merged structure
            cmd.create(pdb_code, "none")
            # Loop over each split state and copy it into the main object
            for i in range(1, n_states + 1):
                state_obj = f"{pdb_code}_{i:04d}"  # e.g. 8AZS_0001
                cmd.copy_to(pdb_code, state_obj)
                cmd.delete(state_obj)  # cleanup
                
                
    # Save the structure as a .pdb file in the published_structure directory
    output_path = os.path.join("Output", "PNG", "contact_maps", f"{pdb_code}.pdb")
    cmd.save(output_path, pdb_code)
    
    structure = read_local_pdb(output_path)
    
    # Remove temp file
    try:
        os.remove(output_path)
        remove_cif()
    except FileNotFoundError:
        pass
    
    # Delete all loaded structures            
    pymol.cmd.delete("all")
    
    return structure

# Extract atom coordinates from PDB
def get_atom_coordinates(structure, pdb_name = None):

    # --- Extract atom info ---
    atom_data = []
    for atom in structure.get_atoms():
        parent = atom.get_parent()
        
        resname = parent.get_resname()
        one_letter = AA_3TO1.get(resname)
        
        if one_letter is None:
            continue # Skip non-protein residues
        
        if atom.get_name() == "CA":  # Only CA atoms
            atom_data.append({
                "PDB": pdb_name,
                "fibril": atom.bfactor,   # fibril number is stored as the B-factor
                "chain": parent.get_parent().id,
                "resno": parent.get_id()[1],
                "resname": AA_3TO1.get(parent.get_resname(), 'X'),
                "x": atom.coord[0],
                "y": atom.coord[1],
                "z": atom.coord[2]
            })

    df = pd.DataFrame(atom_data)
    
    # Removing non-protein components
    

    return df

# Main function to parse PDB file (local or remote) and extract CA coordinates into a DataFrame
def parse_pdb(pdb_file, pdb_name = None):
    
    parser = PDBParser(QUIET=True)

    # Check if PDB is a local file or a PDB code
    if os.path.isfile(pdb_file):
    
        structure = read_local_pdb(pdb_file)

        if structure is None:
            raise ValueError(f"Could not read PDB from {pdb_file}")

        df = get_atom_coordinates(structure, pdb_name)
            
    else:
        pdb_name = pdb_file.upper()

        structure = fetch_pdb_structure(pdb_name)

        df = get_atom_coordinates(structure, pdb_name)

    return df 

#######################################
### 2) Finding the asymmetric units ###
#######################################

# Function to calculate center of mass for a set of coordinates
def calculate_center_of_mass(coords):
    x_com = np.mean(coords['x'])
    y_com = np.mean(coords['y'])
    z_com = np.mean(coords['z'])
    return x_com, y_com, z_com

# Function to calculate Euclidean distance between two points in 3D space
def calculate_distance(point1, point2):
    return np.sqrt((point1[0] - point2[0])**2 + (point1[1] - point2[1])**2 + (point1[2] - point2[2])**2)

# Function to group chains by connectivity using Union-Find algorithm
def group_chains_by_connectivity(df_distances, distance_cutoff=8.0):
    '''
    Group chains that are connected through neighbors using Union-Find algorithm.
    Chains within distance_cutoff are considered neighbors. Chains connected through
    intermediate chains are placed in the same group.
    
    Example: A-B: 4.8Å, B-C: 4.8Å, A-C: 9Å → A, B, C all in same group
    '''
    from collections import defaultdict
    
    # Get all unique chains
    all_chains = set(df_distances['chain1'].unique()) | set(df_distances['chain2'].unique())
    
    # Union-Find data structure
    parent = {chain: chain for chain in all_chains}
    
    def find(x):
        if parent[x] != x:
            parent[x] = find(parent[x])
        return parent[x]
    
    def union(x, y):
        root_x = find(x)
        root_y = find(y)
        if root_x != root_y:
            parent[root_x] = root_y
    
    # Union chains that are within distance cutoff
    for _, row in df_distances.iterrows():
        if row['distance'] < distance_cutoff:
            union(row['chain1'], row['chain2'])
    
    # Group chains by their root parent
    groups = defaultdict(list)
    for chain in all_chains:
        root = find(chain)
        groups[root].append(chain)
    
    return dict(groups)

# Main function to find asymmetric unit chains based on COM and distance between COMs
def find_asymmetric_unit(df):
    ''' 
    Logic: for identical residues on different chains, calculate the COM and the distance betwen each chains COM
    Those with a significantly large distance between COMs will be different chains
    Once each chain has its on fibril number, the lowest z-position chain will be selected
    '''

    for chain in df["chain"].unique():
        chain_coords = df[df["chain"] == chain][["x", "y", "z"]]
        x_com, y_com, z_com = calculate_center_of_mass(chain_coords)
        df.loc[df["chain"] == chain, "x_com"] = x_com
        df.loc[df["chain"] == chain, "y_com"] = y_com
        df.loc[df["chain"] == chain, "z_com"] = z_com

    # Get unique chain COM data
    df_com = df.drop_duplicates(subset=['chain'], keep='first')[['chain', 'x_com', 'y_com', 'z_com']]

    # Calculate and store pairwise distances
    distance_data = []
    for i, row1 in df_com.iterrows():
        for j, row2 in df_com.iterrows():
            if i < j:
                distance = calculate_distance((row1['x_com'], row1['y_com'], row1['z_com']),
                                              (row2['x_com'], row2['y_com'], row2['z_com']))
                distance_data.append({
                    'chain1': row1['chain'],
                    'chain2': row2['chain'],
                    'distance': distance
                })

    df_distances = pd.DataFrame(distance_data)
    
    if len(df["chain"].unique()) > 1:    
        # Group chains by connectivity (distance < 8 Angstroms)
        chain_groups = group_chains_by_connectivity(df_distances, distance_cutoff=8)
        
        print("\nChain groupings (same stacking chain):")
        for group_id, chains in chain_groups.items():
            print(f"  Group: {sorted(chains)}")
    
        # Find lowest z_com for each group and assign fibril number
        group_fibril_numbers = {}
        chains_to_keep = []
        for group_id, chains in chain_groups.items():
            min_z = min(df_com[df_com['chain'].isin(chains)]['z_com'])
            group_fibril_numbers[group_id] = min_z
            
            # Find the chain(s) with the lowest z_com in this group
            lowest_chain = df_com[df_com['chain'].isin(chains) & (df_com['z_com'] == min_z)]['chain'].values
            if len(lowest_chain) > 0:
                chains_to_keep.append(lowest_chain[0])  # Take first if there are ties
    else:
        # Single chain case: keep the only chain
        chains_to_keep = df["chain"].unique().tolist()
    
    # Filter original dataframe to keep only the selected chains
    df_asymmetric_unit = df[df['chain'].isin(chains_to_keep)][['x', 'y', 'z', 'chain', 'resno', 'resname']].copy()

    return df_asymmetric_unit

##########################################################
### 3) calculating pairwise distances between CA atoms ###
##########################################################

def find_min_ca_ca_distance(df):

    # --- Calculate pairwise distances between CA atoms ---
    coords = df[["x", "y", "z"]].to_numpy()
    diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
    CA_distances = np.linalg.norm(diff, axis=2)

    # Convert to long format
    n = len(df)
    CA_distance_df = pd.DataFrame({
        "resno1": np.repeat(df["resno"].values, n),
        "chain1": np.repeat(df["chain"].values, n),
        "resno2": np.tile(df["resno"].values, n),
        "chain2": np.tile(df["chain"].values, n),
        "distance": CA_distances.flatten()
    })

    filtered_distances = CA_distance_df[
        ~((abs(CA_distance_df["resno1"] - CA_distance_df["resno2"]) == 0) &
          (CA_distance_df["chain1"] == CA_distance_df["chain2"]))
    ]
    min_distance = filtered_distances["distance"].min()

    # Further filter to exclude residues within 3 positions of each other in same chain
    CA_distance_df = CA_distance_df[
        ~((abs(CA_distance_df["resno1"] - CA_distance_df["resno2"]) <= 3) &
          (CA_distance_df["chain1"] == CA_distance_df["chain2"]))
    ]

    
    return min_distance, CA_distance_df

###############################################################
### 4) Plotting the PDB coordinates and annotating contacts ###
###############################################################

def calculate_plot_size(df):
    # Extract x and y values to determine scale
    x_vals = df['x'].tolist()
    y_vals = df['y'].tolist()
    
    if not x_vals:
        print("No coordinates to plot")
        return
    
    # Calculate range for scaling
    x_min, x_max = min(x_vals), max(x_vals)
    y_min, y_max = min(y_vals), max(y_vals)
    x_range = x_max - x_min
    y_range = y_max - y_min
    
    # Add small epsilon to avoid division by zero
    epsilon = 1e-6
    x_range = x_range if x_range > epsilon else 1
    y_range = y_range if y_range > epsilon else 1
    
    # Preserve aspect ratio: scale to fixed width, and let height adjust
    plot_width = x_range * 1.2  # Add some padding to width
    aspect_ratio = y_range / x_range
    plot_height = plot_width * aspect_ratio

    return x_min, y_min, x_range, y_range, aspect_ratio, plot_width, plot_height

def add_backbone(df, ax, label = "number", rescol=None):
    
    # Get circle size by setting radius to one third of the minimum CA-CA distance
    # find minimum CA-CA distance
    min_dist, _ = find_min_ca_ca_distance(df)
    print(f"Minimum CA-CA distance: {min_dist:.2f} Å")
    radius = min_dist / 3
    print(f"Circle radius set to: {radius:.2f} Å")
    
    for pos in range(len(df)):
        row = df.iloc[pos]
        resno = row['resno']
        resname = row['resname']
        x = row['x_norm']
        y = row['y_norm']
        
        # Draw lines connecting consecutive CA atoms
        if pos > 0:
            prev_row = df.iloc[pos-1]
            prev_resno = prev_row['resno']
            prev_chain = prev_row['chain']
            prev_x = prev_row['x_norm']
            prev_y = prev_row['y_norm']
            chain = row['chain']
            
            if prev_chain == chain:
                if int(resno) - int(prev_resno) != 1:  # Only connect sequential residues
                    ax.plot([prev_x, x], [prev_y, y], color='black', lw=3, linestyle='--', zorder=1)
                else:
                    ax.plot([prev_x, x], [prev_y, y], color='black', lw=4.5, zorder=1)

        # Setting residue colour
        if not rescol:
            circle_col = "white"
        if rescol in ["type", "property"]:
            if rescol in ["property", "type"]:
                circle_col = CODE_TO_COLOUR.get(resname, "white")
        else: 
            circle_col = "white"
        
        # Draw circle at (x, y)
        circle = plt.Circle((x, y), radius=radius, facecolor=circle_col, edgecolor='black', linewidth=3, zorder=2)
        ax.add_patch(circle)
        
        # Setting residue labels
        if label in ["number", "resno", "pos", "position"]:
            # Add residue label (showing residue number)
            ax.text(x, y, str(resno), ha='center', va='center', fontsize=radius*8, weight='bold', color='black', zorder=3)
        
        elif label in ["letter", "resname", "residue"]:
            # Add residue label (showing residue name)
            ax.text(x, y, str(resname), ha='center', va='center', fontsize=radius*8, weight='bold', color='black', zorder=3)


    # Finished drawing all residues; nothing to return
    return None
    
def add_contacts(dist_data, pdb_data, ax, threshold = None):

    if threshold:

        distance_threshold = float(threshold)

        df = dist_data[dist_data["distance"] < distance_threshold].copy()
        
    else:
        df = dist_data.copy()

    # Drawing contacts
    for i, row in df.iterrows():
        
        resno1 = row["resno1"]
        chain1 = row["chain1"]
        resno2 = row["resno2"]
        chain2 = row["chain2"]

        # Getting coordinates for residue 1
        resno1_xpos = pdb_data[(pdb_data["resno"] == resno1) & (pdb_data["chain"] == chain1)]["x_norm"].values
        resno1_ypos = pdb_data[(pdb_data["resno"] == resno1) & (pdb_data["chain"] == chain1)]["y_norm"].values
        # Getting coordinates for residue 2
        resno2_xpos = pdb_data[(pdb_data["resno"] == resno2) & (pdb_data["chain"] == chain2)]["x_norm"].values
        resno2_ypos = pdb_data[(pdb_data["resno"] == resno2) & (pdb_data["chain"] == chain2)]["y_norm"].values

        if chain1 == chain2:
            ax.plot([resno1_xpos, resno2_xpos], [resno1_ypos, resno2_ypos], color='blue', lw=1.5, zorder=1)
        if chain1 != chain2:
            ax.plot([resno1_xpos, resno2_xpos], [resno1_ypos, resno2_ypos], color='red', lw=1.5, zorder=1)

    return None

def add_legend(ax, rescol, draw_contacts):
    '''
    Add a legend to the plot indicating the meaning of line colors.
    
    Red lines: inter-protofilament contacts
    Black lines: intra-protofilament contacts
    '''
    
    legend_handles = []

    # --- Residue type legend ---
    if rescol in ["type", "property"]:

        unique_props = AA_LOOKUP[["property", "colour"]].drop_duplicates()

        for _, row in unique_props.iterrows():
            circle = mlines.Line2D(
                [0], [0],
                marker = "o",
                color = "black",
                markerfacecolor = row["colour"],
                markersize = 12,
                linestyle = "None",
                label = row["property"]
            )
            legend_handles.append(circle)

    # --- Contact legend ---
    if draw_contacts:
        inter_legend = mlines.Line2D([], [], color = "red", linewidth = 1.5,
                                    label = "Inter-protofilament contacts")
        intra_legend = mlines.Line2D([], [], color = "blue", linewidth = 1.5,
                                    label = "Intra-protofilament contacts")

        legend_handles.extend([inter_legend, intra_legend])

    # --- Draw legend ---
    if legend_handles:
        ax.legend(
            handles = legend_handles,
            loc = "upper left",
            bbox_to_anchor = (1.02, 1),  # push outside right
            borderaxespad = 0,
            fontsize = 14,
            frameon = True
        )
        plt.subplots_adjust(right = 1.2)

### Master function ###

def plot_pdb_coordinates(pdb_file, outfile=None, label=None, rescol=None, draw_contacts=False, threshold=None):
    
    '''
    Master function for generating amyloid folds annotated by contacts
    
    ### Required Arguments ###
    - pdb_file = the location of the PDB file you wish to generate a map for
    - outfile = the save location of the generated map
    
    ### Optional Arguments ###
    - label = How should the residues be labelled; None, "number" for residue number, "letter" for 1 letter AA code
    - rescol = Set how the residues should be coloured: "property" or "type" colours by AA property
    - draw_contacts = Whether to draw contacts on the plot
    - threshold = Maximum distance in Angstroms for Ca-Ca to be considered a contact
    
    '''
    
    print(f"\nProcessing: {pdb_file}")
    
    df = find_asymmetric_unit(parse_pdb(pdb_file))
    
    x_min, y_min, x_range, y_range, aspect_ratio, plot_width, plot_height = calculate_plot_size(df)
    
    # Adjust figure size based on aspect ratio
    fig_width = 0.15*plot_width
    fig_height = fig_width * aspect_ratio
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    
    # Create normalized coordinates
    df_normalized = df.copy()
    df_normalized['x_norm'] = (df['x'] - x_min) / x_range * plot_width
    df_normalized['y_norm'] = (df['y'] - y_min) / y_range * plot_height
    
    add_backbone(df_normalized, ax, label, rescol)
    
    CA_distance_df = pd.DataFrame()

    if draw_contacts:
        # Calculate Ca-Ca distances
        _, CA_distance_df = find_min_ca_ca_distance(df)    
        add_contacts(CA_distance_df, df_normalized, ax, threshold = threshold)
    
    add_legend(ax, rescol, draw_contacts)
        
    
    # Set axis limits with padding
    padding = 2
    ax.set_xlim(-padding, plot_width + padding)
    ax.set_ylim(-padding, plot_height + padding)
    
    ax.set_aspect('equal')
    ax.axis('off')
    
    if outfile:
        if outfile and not outfile.lower().endswith((".png", ".pdf", ".jpg")):
            print("File extension not provided or unaccepted. Defaulting to .png")
            outfile += ".png"
            
        plt.savefig(outfile, dpi=150, bbox_inches='tight')
        print(f"Plot saved to {outfile}")
        plt.close(fig)
    else:
        plt.show()
        plt.close(fig)
    
    return CA_distance_df

#######################################################################################################################

if __name__ == "__main__":

    # # Residue contacts 
    # residue_contacts = os.path.join("Output", "stable_regions", "residue_contacts", "residue_distances_by_cluster_group.csv")
    # # regions contacts
    # region_contacts = os.path.join("Output", "stable_regions", "region_contacts", "region_distances_by_cluster_group.csv")

    ### Asyn Tests ###
    #pdb_file = os.path.join("Output", "PDBs", "asymetric_unit", "6cu7_asym_unit.pdb")
    #pdb_file = os.path.join("Output", "PDBs", "asymetric_unit", "6ssx_asym_unit.pdb")
    #pdb_file = os.path.join("Output", "PDBs", "asymetric_unit", "8a9l_asym_unit.pdb")
    #pdb_file = os.path.join("Output", "PDBs", "unique_chains", "8y2p.pdb")
    #pdb_file = os.path.join("Output", "PDBs", "unique_chains", "6cu7.pdb")
    pdb_file = "7lc9"
    
    ### IAPP Tests ###
    #pdb_file = os.path.join("Output", "PDBs", "asymetric_unit", "8awt_asym_unit.pdb")
    #pdb_file = os.path.join("Output", "PDBs", "asymetric_unit", "9gz6_asym_unit.pdb")
    #pdb_file = "8awt"
    
    plot_pdb_coordinates(pdb_file, outfile="Output/PNG/test.png", label="number", rescol="type", 
                         draw_contacts=False, threshold=10)