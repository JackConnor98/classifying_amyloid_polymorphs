import subprocess
import customtkinter as ctk
import tkinter as tk
from tkinter import messagebox

# -----------------------------------------------------------
# Utility
# -----------------------------------------------------------

def safe(val, default = "NA"):
    if val is None:
        return default
    s = str(val).strip()
    return s if s != "" else default


# -----------------------------------------------------------
# Pipeline execution
# -----------------------------------------------------------

def run_pipeline():
    selected_jobs = []
    scrape = scrape_var.get()
    pdb = pdb_var.get()
    validation = validation_var.get()
    rmsd = rmsd_var.get()
    thermodynamics = thermodynamics_var.get()
    stable_regions = stable_regions_var.get()
    b_strand = b_strand_var.get()
    png = png_var.get()

    penalty = safe(penalty_var.get(), "0")
    custom_cut_height = safe(custom_cut_height_var.get(), "0")
    remove_poorly_resolved = str(remove_poorly_resolved_var.get())
    ph = safe(ph_var.get(), "7.4")
    temp = safe(temp_var.get(), "298")
    ionstrength = safe(ion_var.get(), "0.15")
    window = safe(window_var.get(), "3")
    min_stable = safe(min_stable_region_size_var.get(), "0")
    distance_threshold = safe(distance_threshold_var.get(), "10.8")
    min_length = safe(min_length_var.get(), "0")
    png_colouring = safe(png_colouring_var.get(), "0")
    png_palette = safe(png_palette_var.get(), "tv_blue, tv_red, tv_green, tv_orange")

    try:
        if scrape == 1:
            subprocess.run(["python", "Scripts/scrape/amyloid_atlas_scraper_gui.py"], check = True)
            subprocess.run(["python", "Scripts/scrape/plot_ordered_residues.py"], check = True)
            selected_jobs.append("1.Web Scrape")

        if pdb == 1:
            subprocess.run(["python", "Scripts/PDB/fetch_pdb_isolate_chains.py"], check = True)
            subprocess.run(["python", "Scripts/PDB/Rg_plotting.py"], check = True)
            subprocess.run(["python", "Scripts/PDB/sequence_alignment.py"], check = True)
            selected_jobs.append("2.Analyse PDBs")

        if validation == 1:
            subprocess.run(["python", "Scripts/validation/Q_score_scraper.py"], check = True)
            subprocess.run(["python", "Scripts/validation/validating_structures.py"], check = True)
            selected_jobs.append("3.Validation")

        if rmsd == 1:
            subprocess.run(["python", "Scripts/RMSD/unique_chain_alignment.py", penalty], check = True)
            subprocess.run(["python", "Scripts/RMSD/RMSD_analysis.py", custom_cut_height], check = True)
            selected_jobs.append("4.RMSD")

        if thermodynamics == 1:
            subprocess.run(["python", "Scripts/thermodynamics/EMDB_scraper.py"], check = True)
            subprocess.run(["python", "Scripts/thermodynamics/extend_fibril_layers.py"], check = True)
            subprocess.run([
                "python",
                "Scripts/thermodynamics/foldx_analysis.py",
                ph,
                temp,
                ionstrength
            ], check = True)
            subprocess.run(["python", "Scripts/thermodynamics/thermodynamics_plotting.py", remove_poorly_resolved], check = True)
            selected_jobs.append("5.Thermodynamics")

        if stable_regions == 1:
            subprocess.run(["python", "Scripts/stable_regions/defining_stable_regions.py", window, min_stable], check = True)
            subprocess.run(["python", "Scripts/stable_regions/stable_region_distances.py", distance_threshold], check = True)
            selected_jobs.append("6.Stable Regions")

        if b_strand == 1:
            subprocess.run(["python", "Scripts/beta_strand/beta_strand_analysis.py", min_length], check = True)
            selected_jobs.append("7.Beta-Strands")

        if png == 1:
            subprocess.run([
                "python",
                "Scripts/PNG/asymmetric_unit_png_generator.py",
                png_colouring,
                png_palette
            ], check = True)
            subprocess.run(["python", "Scripts/PNG/asymmetric_unit_figure_maker.py"], check = True)
            subprocess.run(["python", "Scripts/PNG/stable_region_colouring.py"], check = True)
            subprocess.run(["python", "Scripts/PNG/stable_region_png_maker.py"], check = True)
            subprocess.run(["python", "Scripts/PNG/cluster_group_and_stable_regions_figure.py"], check = True)
            subprocess.run(["python", "Scripts/PNG/polarity_map_scraper.py"], check = True)
            selected_jobs.append("7.PNG")

        messagebox.showinfo(
            "Done",
            "Selected pipeline steps completed successfully:\n\n" + "\n".join(selected_jobs)
        )


    except subprocess.CalledProcessError as e:
        messagebox.showerror("Error", f"A step failed: {e}\nSee console for details.")


# -----------------------------------------------------------
# GUI Build
# -----------------------------------------------------------

ctk.set_appearance_mode("dark")
ctk.set_default_color_theme("blue")

app = ctk.CTk()
app.title("Calypso: Structural and Thermodynamic Analysis of Amyloid Fibrils")
app.geometry("820x820")

scroll = ctk.CTkScrollableFrame(app, width = 780, height = 780)
scroll.pack(padx = 10, pady = 10, fill = "both", expand = True)

# -----------------------------------------------------------
# Helper functions
# -----------------------------------------------------------

def section(title):
    frame = ctk.CTkFrame(scroll)
    label = ctk.CTkLabel(frame, text = title, font = ("Arial", 26, "bold"))
    label.pack(anchor = "w", pady = (4, 8))
    frame.pack(fill = "both", padx = 10, pady = 10)
    return frame

def make_label(parent, text, size = 16, weight = "normal", side = "left", pad_y = (4, 4)):
    font = ("Arial", size, weight)
    label = ctk.CTkLabel(parent, text = text, font = font)
    label.pack(side = side, pady = pad_y)
    return label

def make_label_new_line(parent, text, size = 16, weight = "normal", pad_y = (4, 4), pad_x = 6):
    label = ctk.CTkLabel(
        parent,
        text = text,
        justify = "left",
        anchor = "w",
        font = ("Arial", size, weight)
    )
    label.pack(anchor = "w", padx = pad_x, pady = pad_y)
    return label

def make_checkbox(parent, text, variable, size = 16, weight = "normal", pad_y = 4):
    checkbox = ctk.CTkCheckBox(parent, text = text, variable = variable, font = ("Arial", size, weight))
    checkbox.pack(anchor = "w", pady = pad_y)
    return checkbox

def make_radiobutton(parent, text, variable, value, size = 16, weight = "normal", side = "left", pad_x = 6):
    radiobutton = ctk.CTkRadioButton(parent, text = text, variable = variable, value = value, font = ("Arial", size, weight))
    radiobutton.pack(side = side, padx = pad_x)
    return radiobutton

def make_entry(parent, variable, size = 16, weight = "normal", side = "left", width = 120, pad_x = 6, pad_y = 4):
    entry = ctk.CTkEntry(parent, textvariable = variable, width = width, font = ("Arial", size, weight))
    entry.pack(side = side, padx = pad_x, pady = pad_y)
    return entry

def make_button(parent, text, command, size = 20, weight = "bold", height = 40, corner_radius = 8):
    button = ctk.CTkButton(parent, text = text, command = command, height = height, corner_radius = corner_radius, 
                           font = ("Arial", size, weight))
    button.pack(pady = 10)
    return button

# -----------------------------------------------------------
# 1. Web scrape
# -----------------------------------------------------------

scrape_frame = section("1. Web Scrape PDB Names From Amyloid Atlas")
scrape_var = ctk.IntVar(value = 1)
make_checkbox(scrape_frame, "Run web scraping", scrape_var)


# -----------------------------------------------------------
# 2. Analyse PDBs
# -----------------------------------------------------------

pdb_frame = section("2. Analyse PDBs")
pdb_var = ctk.IntVar(value = 1)
make_checkbox(pdb_frame, "Fetch and analyse PDBs", pdb_var)

# -----------------------------------------------------------
# 3. Validation
# -----------------------------------------------------------

validation_frame = section("3. Validation")
validation_var = ctk.IntVar(value = 1)
make_checkbox(validation_frame, "Run validation (Q-scores and selection)", validation_var)

# -----------------------------------------------------------
# 4. RMSD
# -----------------------------------------------------------

rmsd_frame = section("4. RMSD")
rmsd_var = ctk.IntVar(value = 1)
make_checkbox(rmsd_frame, "Run RMSD alignment and clustering", rmsd_var)

penalty_var = ctk.StringVar(value = "0")
custom_cut_height_var = ctk.StringVar(value = "0")

row = ctk.CTkFrame(rmsd_frame)
row.pack(fill = "x", pady = 4)
make_label_new_line(row, "Set a penalty for non-overlapping residues in RMSD calculations.\nPenalty = Mean RMSD + (penalty * SD)")
make_entry(row, penalty_var, width = 50)

row2 = ctk.CTkFrame(rmsd_frame)
row2.pack(fill = "x", pady = 4)
make_label_new_line(row2, "Specify the Euclidean Distance to be used for as the cut height to determine RMSD cluster groups.\nFor the first run, 0 is reccomended as this will default to using the mean Euclidean distance.")
make_entry(row2, custom_cut_height_var, width = 50)

# -----------------------------------------------------------
# 5. Thermodynamics
# -----------------------------------------------------------

thermo_frame = section("5. Thermodynamics")
thermodynamics_var = ctk.IntVar(value = 1)
make_checkbox(thermo_frame, "Run thermodynamic analysis (FoldX)", thermodynamics_var)

remove_poorly_resolved_var = ctk.IntVar(value = 1)

make_label_new_line(thermo_frame, "Remove individual poorly resolved residues from otherwise well resolved PDBs?")

remove_poorly_resolved_var = ctk.StringVar(value = "Yes")  # default choice

radio_frame = ctk.CTkFrame(thermo_frame)
radio_frame.pack(anchor = "w", pady = 2)
make_radiobutton(radio_frame, "Yes", remove_poorly_resolved_var, value = "Yes")
make_radiobutton(radio_frame, "No", remove_poorly_resolved_var, value = "No")

ph_var = ctk.StringVar(value = "7.4")
temp_var = ctk.StringVar(value = "298")
ion_var = ctk.StringVar(value = "0.15")

grid = ctk.CTkFrame(thermo_frame)
grid.pack(fill = "x", pady = 4)

make_label_new_line(grid, "Here you can control the parameters used for FoldX calculations")

for label_text, var in [
    ("pH", ph_var),
    ("Temperature (K)", temp_var),
    ("Ionic Strength", ion_var)
]:
    make_label(grid, label_text)
    make_entry(grid, var, width = 80, side = "left")
    
# -----------------------------------------------------------
# 6. Stable regions
# -----------------------------------------------------------

stable_frame = section("6. Stable Regions")
stable_regions_var = ctk.IntVar(value = 1)
make_checkbox(stable_frame, text = "Define and calculate stable regions", variable = stable_regions_var)

grid2 = ctk.CTkFrame(stable_frame)
grid2.pack(fill = "x", pady = 4)

window_var = ctk.StringVar(value = "3")
min_stable_region_size_var = ctk.StringVar(value = "0")

make_label_new_line(grid2, text = "Specify the window size for the sliding window rolling average used for stable region detection",)
make_entry(grid2, window_var, width = 50)

distance_threshold_var = ctk.StringVar(value = "10.8")
row3 = ctk.CTkFrame(stable_frame)
row3.pack(fill = "x", pady = 4)
make_label_new_line(row3, text = "Distance Threshold for Stabilising Regions/Residues Contacts (Å)")
make_entry(row3, distance_threshold_var, width = 50)

# -----------------------------------------------------------
# 7. Beta strand
# -----------------------------------------------------------

beta_frame = section("7. Beta Strand Analysis")
b_strand_var = ctk.IntVar(value = 1)
make_checkbox(beta_frame, text = "Analyse beta sheets", variable = b_strand_var)

min_length_var = ctk.StringVar(value = "0")
row4 = ctk.CTkFrame(beta_frame)
row4.pack(fill = "x", pady = 4)
make_label_new_line(row4, text = "Set the minimum number of residues required to be considered a beta-strand")
make_entry(row4, min_length_var, width = 50)


# -----------------------------------------------------------
# 8. PNG generation
# -----------------------------------------------------------

png_frame = section("8. PNG Generation")
png_var = ctk.IntVar(value = 1)
make_checkbox(png_frame, text = "Generate PNGs and figures", variable = png_var)

png_colouring_var = ctk.StringVar(value = "0")
selected_colour_var = ctk.StringVar(value = "tv_blue")
png_palette_var = ctk.StringVar(value = "tv_blue, tv_red, tv_green, tv_orange")

row5 = ctk.CTkFrame(png_frame)
row5.pack(fill = "x", pady = 4)


label_entry_frame = ctk.CTkFrame(row5)
label_entry_frame.pack(fill="x", pady = 2)
make_label(label_entry_frame, text = "Please select a colour scheme (0-3)", side = "left")
make_entry(label_entry_frame, png_colouring_var, width = 50, side = "left")

bullet_frame = ctk.CTkFrame(row5)
bullet_frame.pack(fill = "x", pady = 2, anchor = "w")

make_label_new_line(bullet_frame, text = "0 = Single Colour")
make_label_new_line(bullet_frame, text = "1 = Protein Domain (Tau, Asyn or A-beta Only)")
make_label_new_line(bullet_frame, text = "2 = Protofilament")
make_label_new_line(bullet_frame, text = "3 = Amyloid Fold")


row6 = ctk.CTkFrame(png_frame)
row6.pack(fill = "x", pady = 4)
make_label_new_line(row6, text = "Please provide a comma seperated colour palette for Pymol:")
make_entry(row6, png_palette_var, width = 400)

# -----------------------------------------------------------
# Run button
# -----------------------------------------------------------

make_button(scroll, text = "Run Selected Jobs", command = run_pipeline)

# -----------------------------------------------------------
# Main loop
# -----------------------------------------------------------

app.mainloop()
