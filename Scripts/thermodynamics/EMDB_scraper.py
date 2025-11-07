from bs4 import BeautifulSoup
import requests 
import csv
import os
import PyPDF2
import re
import gzip
import shutil
import pandas as pd

##### STEP 1 - Converting PDB ids into EMDB numbers #####

# Getting EMDB names from PDB names
url = "https://www.rcsb.org/structure/"

# List of PDBs to find EMDB data for (Only using those that passed validation)
pdb_df = pd.read_csv(os.path.join("Output", "Validation", "high_resolution_pdb_ids.csv"))
pdb_ids = pdb_df['pdb_id'].tolist()
pdb = list(set([x.split("_")[0] for x in pdb_ids]))

emdb_results = []

for i in pdb:
    pdb_query = url + i
    response = requests.get(pdb_query)
    doc = BeautifulSoup(response.text, "html.parser")
    try:
        header_emdb = doc.find("li", id="header_emdb")
        strong_element = header_emdb.find("strong")
        emdb_text = strong_element.get_text(strip=True)
        emdb_number = emdb_text.split("-")[1]
        emdb_number = emdb_number[:-1]
        emdb_results.append((i, emdb_number))
        print("PDB: " + i + "   EMDB: " + emdb_number + "\n")
    except AttributeError:
        emdb_results.append((i, "NA"))
        print("PDB: " + i + "   EMDB: Not Found\n")

##### Step 2 - Using the EMDB numbers to scrape helical twist and helical rise from EMDB and calculate half pitch #####

emdb_numbers = [emdb for pdb, emdb in emdb_results]
results = []
missing_data = []

for i in range(len(emdb_numbers)):
    if emdb_numbers[i] != "NA":
        url = f"https://files.wwpdb.org/pub/emdb/validation_reports/EMD-{emdb_numbers[i]}/emd_{emdb_numbers[i]}_full_validation.pdf.gz"
        response = requests.get(url, stream=True)
        gz_path = "downloaded.pdf.gz"
        pdf_path = "downloaded.pdf"
        second_page_text = ""
        if response.status_code == 200:
            with open(gz_path, 'wb') as f:
                for chunk in response.iter_content(1024):
                    f.write(chunk)
            # Check if file is gzipped
            with open(gz_path, 'rb') as f:
                magic = f.read(2)
            is_gzipped = magic == b'\x1f\x8b'
            try:
                if is_gzipped:
                    with gzip.open(gz_path, 'rb') as f_in:
                        with open(pdf_path, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
                else:
                    os.rename(gz_path, pdf_path)
                try:
                    pdf = PyPDF2.PdfReader(pdf_path)
                except Exception as e:
                    print(f"Could not read PDF for {emdb_numbers[i]}: {e}")
                    results.append((pdb[i], emdb_numbers[i], "NA", "NA", "NA", "NA"))
                    continue
            except Exception as e:
                print(f"Could not extract PDF for {emdb_numbers[i]}: {e}")
                results.append((pdb[i], emdb_numbers[i], "NA", "NA", "NA", "NA"))
                continue

            if len(pdf.pages) >= 2:
                second_page = pdf.pages[1]
                second_page_text = second_page.extract_text()
            else:
                print("The PDF does not have at least 2 pages.")
                results.append((pdb[i], emdb_numbers[i], "NA", "NA", "NA", "NA"))
                continue
        else:
            print("Failed to download the PDF.")
            results.append((pdb[i], emdb_numbers[i], "NA", "NA", "NA", "NA"))
            continue

        twist_pattern = r"wist=(-?\d+\.\d+)"
        rise_pattern = r"rise=([\d.]+)"
        twist_match = re.search(twist_pattern, second_page_text)
        rise_match = re.search(rise_pattern, second_page_text)

        if twist_match and rise_match:
            extracted_twist = float(twist_match.group(1))
            extracted_rise = float(rise_match.group(1))
            if abs(extracted_twist) > 100:
                twist_converted = False
                if extracted_twist < 0:
                    extracted_twist = round(180 - abs(extracted_twist), 2)
                    twist_converted = True
                if extracted_twist > 0 and twist_converted == False:
                    extracted_twist = round(extracted_twist - 180, 2)
            half_pitch = round(abs((180/extracted_twist)*extracted_rise), 2)
            handed = "right" if extracted_twist > 0 else "left"
            results.append((pdb[i], emdb_numbers[i], extracted_twist, extracted_rise, half_pitch, handed))
            print("Helical Twist: ", extracted_twist)
            print("\nHelical Rise: ", extracted_rise)
            print("\nHalf Pitch: ", half_pitch)
            print("\nHandedness: ", handed)
        else:
            results.append((pdb[i], emdb_numbers[i], "NA", "NA", "NA", "NA"))
            missing_data.append(pdb[i])
            print("No match found.")

        # Clean up files
        if os.path.exists(pdf_path):
            os.remove(pdf_path)
        if os.path.exists(gz_path):
            os.remove(gz_path)
    else:
        results.append((pdb[i], emdb_numbers[i], "NA", "NA", "NA", "NA"))

# Create and write the results to a text file
with open(os.path.join("Output", "emdb_data.txt"), "w") as file:
    file.write("pdb\temdb\ttwist\trise\thalf_pitch\thandedness\n")
    for pdb_value, emdb_value, twist_value, rise_value, half_pitch, handed in results:
        file.write(f"{pdb_value}\t{emdb_value}\t{twist_value}\t{rise_value}\t{half_pitch}\t{handed}\n")

print("Results saved to 'Output/emdb_data.txt'")

#######################################################
### Adding PDBs missing EMDB data to exclusion list ###
#######################################################

excluded_file = os.path.join("Output", "excluded_list.txt")
file_exists = os.path.exists(excluded_file) # Check if file exists

# Getting a list of pdbs already added to exclusion list to avoid duplicate entries
existing_pdbs = set()
if file_exists:
    with open(excluded_file, "r") as f:
        for line in f:
            if line.strip() and not line.startswith("PDB ID"):
                existing_pdbs.add(line.split("\t")[0])  # Get PDB ID before first tab

# Adding missing PDBs with missing Q-scores to excluded list
if missing_data:
    with open(excluded_file, "a") as f:
        # Write header only if the file didn't exist before
        if not file_exists:
            f.write("PDB ID\tReason\n")

        # Write your excluded entries
        for pdb in missing_data:
            if pdb not in existing_pdbs:  
                line = f"{pdb}\tEMDB data not found\n"
                f.write(line)

    print(f"\nEMDB data not found for {len(missing_data)} PDB(s)")
else:
    print("\n🎉 EMDB data found for all PDBs.")

