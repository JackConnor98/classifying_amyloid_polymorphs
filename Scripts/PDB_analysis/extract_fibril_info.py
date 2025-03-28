import os
import pandas as pd
import numpy as np

# Specify the directory you want to set as the working directory
working_directory = '/mnt/c/Users/Jack_/Desktop/Leeds/Papers/structure_thermodynamics_polymorphism/RSMD_Thermodynamics_Analysis'

# Change the current working directory
os.chdir(working_directory)

# Creating pdb_info_df for use in further analysis
df = pd.read_csv('Output/PDBs/COM_and_fibril.csv')

current_pdb_df = df[df["PDB"] == "6l1u"]

current_pdb_df["polymorph"] = np.nan

# Calculate mean Rg for each fibril
fibril_mean_Rg = current_pdb_df.groupby('fibril')['Rg'].mean().reset_index()
fibril_mean_Rg.columns = ['fibril', 'mean_Rg']

print("Fibril Mean Rg = \n", fibril_mean_Rg)

# Calculate 10% confidence intervals for each fibril's mean Rg
fibril_mean_Rg['lower_bound'] = fibril_mean_Rg['mean_Rg'] * 0.9
fibril_mean_Rg['upper_bound'] = fibril_mean_Rg['mean_Rg'] * 1.1

# Check for distinct mean Rg values among fibrils
distinct_fibrils = []
similar_pairs = []
for i, fibril in fibril_mean_Rg.iterrows():
    distinct = True
    similar_to = []
    for j, compare_fibril in fibril_mean_Rg.iterrows():
        if i != j:
            if fibril['lower_bound'] <= compare_fibril['mean_Rg'] <= fibril['upper_bound']:
                distinct = False
                similar_to.append(compare_fibril['fibril'])
    if distinct:
        distinct_fibrils.append(fibril['fibril'])
    else:
        if similar_to:
            similar_pairs.append((fibril['fibril'], similar_to))



# Assign unique polymorph numbers to distinct fibrils

current_pdb_df["polymorph"] = np.nan

polymorph_counter = 1
for fibril in distinct_fibrils:
    current_pdb_df.loc[current_pdb_df['fibril'] == fibril, 'polymorph'] = polymorph_counter
    polymorph_counter += 1

# If there is only one polymorph - set all to 1
if current_pdb_df['polymorph'].isnull().all():
    current_pdb_df['polymorph'] = 1

# Ensure that 'polymorph' is treated as a string
current_pdb_df['polymorph'] = current_pdb_df['polymorph'].astype(int).astype(str)

# Create the 'pdb_id' column based on the number of polymorphs
if distinct_fibrils:
    current_pdb_df['pdb_id'] = current_pdb_df['PDB'] + '_' + current_pdb_df['polymorph']
else:
    current_pdb_df['pdb_id'] = current_pdb_df['PDB']



print(similar_pairs)
print("\n\n")
print(distinct_fibrils)
print("\n\n")
print(current_pdb_df)