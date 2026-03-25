from annotate_contacts import *
import os

# Creating output dir
save_dir = os.path.join("Output", "PNG", "contact_maps")

if not os.path.exists(save_dir):
    os.makedirs(save_dir)

file_path = os.path.join("Output", "pdb_names.txt")

# PDB list
with open(file_path, 'r') as file:
    lines = file.readlines()
    # Skip the first line (header) and strip whitespace
    pdb_names = [line.strip() for line in lines[1:]]    
    
launch_pymol()    
    
all_dist_data = pd.DataFrame()
    
for pdb in pdb_names:
    
    distance_df = plot_pdb_coordinates(pdb, outfile=os.path.join(save_dir, f"{pdb}.png"), label="number", rescol="type",
                                       draw_contacts=False, threshold=10)
    
    all_dist_data = pd.concat([all_dist_data, distance_df], ignore_index=True)


if not all_dist_data.empty:
    all_dist_data.to_csv(os.path.join(save_dir, "dist_data.csv"), index=False)   

close_pymol()
