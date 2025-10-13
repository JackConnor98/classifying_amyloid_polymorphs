from bs4 import BeautifulSoup
import requests 
import os
import pandas as pd

output_dir = "Output"
url = "https://people.mbi.ucla.edu/sawaya/amyloidatlas/"

# Creating Output Directory if it does not exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# Reading in a HTML from a website
def fetch_html(url):
    response = requests.get(url)
    doc = BeautifulSoup(response.text, "html.parser")
    
    return doc

def parse_table_data(doc):

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

    return headers, row_data 

def get_user_protein_pattern(df, silent=False):
    protein_values = df["Protein"].dropna().unique()
    protein_values.sort()

    if not silent:
        print("\nUnique Proteins Found:\n")
        for idx, protein in enumerate(protein_values, 1):
            print(f"{idx:2}: {protein}")
        
    print("\nEnter a pattern matching the protein name you're interested in or type 'q' to quit.")

    if not silent:
        print("\ne.g. for α-synuclein, enter 'synuclein'")
        print("\ne.g. for Amyloid-β(1-42), enter 'amyloid-'")

    pattern = input("\n>>> ").strip()

    if pattern.lower() == "q":
        return None

    return pattern

def pattern_match(pattern, df, continue_selection=True):

    # filter DataFrame by Protein column matching pattern 
    selected_rows = df[df["Protein"].str.contains(str(pattern), case = False)]

    while selected_rows.empty:

        if pattern is None:
            print("👋 Exiting program. Bye.")
            return
        
        print(f"\n🚫 No PDBs found matching the pattern '{pattern}'. Please try again.")
        pattern = get_user_protein_pattern(df, silent=True)
        selected_rows = df[df["Protein"].str.contains(str(pattern), case = False)]

    # Extracting PDB names from selected rows
    pdb_names = selected_rows["PDB ID"].tolist()

    # Adding column name to start of the list
    pdb_names.insert(0, "PDB")

    # Showing the selected PDBs and their metadata
    print("\nSelected PDBs Metadata:")
    print(selected_rows)

    # Printing number of PDBs
    print("\nYou Selected {} PDBs".format(len(pdb_names)))

    while True:
        happy = input("Are you sure you want to continue with these PDBs? (y/n) >>> ").strip().lower()
        if happy.lower() == 'y':
            return selected_rows, pdb_names, False
        elif happy.lower() == 'n':
            print("Please try again with a different pattern.")
            return selected_rows, pdb_names, True
        else:
            print("🚫 Invalid input. Please enter 'y' or 'n'.")

def save_data(selected_rows, pdb_names):

    # Saving filtered DataFrame to a text file
        selected_rows.to_csv(os.path.join(output_dir, "selected_pdbs_metadata.txt"), index = False, header = True, sep = "\t")

        # Saving PDB names to a text file
        with open(os.path.join(output_dir, "pdb_names.txt"), 'w') as file:
            # write each value to a new line
            for value in pdb_names:
                file.write(value + '\n')

def main():
    try:
        doc = fetch_html(url)
        headers, rows = parse_table_data(doc)
        df = pd.DataFrame(rows, columns=headers)

        # Selecting only cryoEM structures
        cryoEM_df = df[df['Method'].str.contains('cryoEM', case=False)]
        
        pattern = get_user_protein_pattern(df)

        continue_selection = True

        while continue_selection and pattern is not None:

            selected_rows, pdb_names, continue_selection = pattern_match(pattern, cryoEM_df)

            if not continue_selection:

                # Save the selected data
                save_data(selected_rows, pdb_names)

                break

            # If user wants to try again, get a new pattern
            pattern = get_user_protein_pattern(df)

        if pattern is None:
            print("👋 Exiting program. Bye.")
            return

    
    except Exception as e:
        print(f"💥 Something went wrong: {e}")

if __name__ == "__main__":
    main()
