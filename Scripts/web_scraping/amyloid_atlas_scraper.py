#%%

from bs4 import BeautifulSoup
import requests 
import csv
import os
import tkinter as tk
from tkinter import simpledialog
import pandas as pd

# alpha-syn = synuclein
# tau = Tau
# abeta = Amyloid-
# IAPP = IAPP
# Light chain antibody = light chain
# TDP-43 = tdp


# Reading in a HTML from a website
url = "https://people.mbi.ucla.edu/sawaya/amyloidatlas/"
result = requests.get(url)
doc = BeautifulSoup(result.text, "html.parser")

# Accessing Table Data = <td><big><centre> data of interest </centre></big></td> in html terms
table_values = doc.find_all("td")

# initializing headers
headers = []

# Initializing lists
row_data = []
current_row = []


# The first 9 values in table_values are the column headers
# Of which 6 and 7 we do not want (they contain images)
# for loop extracts information from every row except for index 5 and 6 because of negative indexing
   
for i in range(0,len(table_values)):
    if i < 9 and i % 9 != 5 and i % 9 != 6:
        headers.append(table_values[i].get_text(strip=True))
    
    if i >= 9 and i % 9 != 5 and i % 9 != 6:
        x = table_values[i].get_text(strip=True)

        current_row.append(x)
        #current_row.append(table_values[i].get_text(strip=True))

    if i > 8 and i % 9 == 8:
        row_data.append(current_row)
        current_row = []    

# adding headers to the beginning of the list
row_data.insert(0,headers) 


# Creating Output Directory if it does not exist
if not os.path.exists("Output"):
    os.makedirs("Output")

# Function to show the protein column and get user input
def show_protein_column_and_get_input(protein_data):
    # Create the main window
    root = tk.Tk()
    root.withdraw()  # Hide the main window

    # Create a popup window
    popup = tk.Toplevel(root)
    popup.title("Protein Column")

    # Create a scrollbar
    scrollbar = tk.Scrollbar(popup)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    # Create a text widget
    text_widget = tk.Text(popup, wrap=tk.NONE, yscrollcommand=scrollbar.set)
    text_widget.pack(expand=True, fill=tk.BOTH)

    # Insert protein data into the text widget
    for protein in protein_data:
        text_widget.insert(tk.END, protein + '\n')

    # Configure the scrollbar
    scrollbar.config(command=text_widget.yview)

    # Create a simple dialog to get user input
    user_input = simpledialog.askstring("Selecting PDBs", "Please enter a pattern to match your protein of interest:", parent=popup)

    # Close the popup window
    popup.destroy()
    root.quit()

    return user_input

# Convert row_data to a DataFrame
headers = row_data[0]
data = row_data[1:]
df = pd.DataFrame(data, columns=headers)

# Extract the 'Protein' column from row_data
protein_index = row_data[0].index('Protein')
protein_data = [row[protein_index] for row in row_data[1:]]

# Show the protein column and get user input
pattern = show_protein_column_and_get_input(protein_data)

# Convert row_data to a DataFrame
headers = row_data[0]
data = row_data[1:]
df = pd.DataFrame(data, columns=headers)


# Selecting only cryoEM structures
cryoEM_df = df[df['Method'].str.contains('cryoEM', case=False)]

# filter DataFrame by Protein column matching pattern "syn"
selected_rows = cryoEM_df[cryoEM_df["Protein"].str.contains(str(pattern), case = False)]

# Saving filtered DataFrame to a text file
selected_rows.to_csv("Output/selected_pdbs_metadata.txt", index = False, header = True, sep = "\t")

# Extracting PDB names from selected rows
pdb_names = selected_rows["PDB ID"].tolist()

# Printing number of PDBs
print("\n You Selected {} PDBs".format(len(pdb_names)))

# Adding column name to start of the list
pdb_names.insert(0, "PDB")

# Saving PDB names to a text file
with open('Output/pdb_names.txt', 'w') as file:
    # write each value to a new line
    for value in pdb_names:
        file.write(value + '\n')

# %%
