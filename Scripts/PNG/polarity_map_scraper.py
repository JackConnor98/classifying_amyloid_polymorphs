import csv
import os
import pandas as pd
import requests

########################################
### Downlaoading polarity map images ###
########################################

# Specifying save directory
save_folder = os.path.join("Output", "PNG", "polarity_maps")

# Create the save folder if it doesn't exist
os.makedirs(save_folder, exist_ok=True)

# Get PDB names
df = pd.read_csv('Output/pdb_names.txt', delimiter='\t')
pdb_names = df["PDB"]

# Setting the prefix and suffix for the polarity map url
url_prefix = "https://people.mbi.ucla.edu/sawaya/amyloidatlas/residuefull/"
url_suffix = "_polarimap.png"


for pdb in pdb_names:
    
    # Create the image url
    img_url = url_prefix + pdb + url_suffix

    try:
        img_data = requests.get(img_url).content

        # Save the image to the specified folder
        img_name = os.path.join(save_folder, pdb + ".png")
        with open(img_name, 'wb') as img_file:
            img_file.write(img_data)

        print(f"{pdb} Polarity Map Downloaded Successfully.")
    
    except Exception as e:
        print(f"Error downloading image from '{img_url}': {e}")
