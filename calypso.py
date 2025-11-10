import subprocess
import customtkinter as ctk
from tkinter import messagebox
from PIL import Image, ImageTk

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

# # ----------------------------------------------------------
# # ------------------ ADD BACKGROUND IMAGE ------------------
# # ----------------------------------------------------------

# bg_image = Image.open("Figures/Calypso.png")
# bg_image = bg_image.resize((780, 780))  # resize to match window
# bg_photo = ImageTk.PhotoImage(bg_image)

# bg_label = ctk.CTkLabel(scroll, image=bg_photo, text="")
# bg_label.place(x=0, y=0, relwidth=1, relheight=1)
# bg_label.lower()  # send it behind all other widgets in scroll

# -----------------------------------------------------------
# Helper for section frames
# -----------------------------------------------------------

def section(title):
    frame = ctk.CTkFrame(scroll)
    label = ctk.CTkLabel(frame, text = title, font = ("Arial", 20, "bold"))
    label.pack(anchor = "w", pady = (4, 8))
    frame.pack(fill = "both", padx = 10, pady = 10)
    return frame


# -----------------------------------------------------------
# 1. Web scrape
# -----------------------------------------------------------

scrape_frame = section("1. Web Scrape PDB Names From Amyloid Atlas")
scrape_var = ctk.IntVar(value = 1)
ctk.CTkCheckBox(scrape_frame, text = "Run web scraping", variable = scrape_var).pack(anchor = "w", pady = 4)


# -----------------------------------------------------------
# 2. Analyse PDBs
# -----------------------------------------------------------

pdb_frame = section("2. Analyse PDBs")
pdb_var = ctk.IntVar(value = 1)
ctk.CTkCheckBox(pdb_frame, text = "Fetch and analyse PDBs", variable = pdb_var).pack(anchor = "w", pady = 4)


# -----------------------------------------------------------
# 3. Validation
# -----------------------------------------------------------

validation_frame = section("3. Validation")
validation_var = ctk.IntVar(value = 1)
ctk.CTkCheckBox(validation_frame, text = "Run validation (Q-scores and selection)", variable = validation_var).pack(anchor = "w", pady = 4)


# -----------------------------------------------------------
# 4. RMSD
# -----------------------------------------------------------

rmsd_frame = section("4. RMSD")
rmsd_var = ctk.IntVar(value = 1)
ctk.CTkCheckBox(rmsd_frame, text = "Run RMSD alignment and clustering", variable = rmsd_var).pack(anchor = "w", pady = 4)

penalty_var = ctk.StringVar(value = "0")
custom_cut_height_var = ctk.StringVar(value = "0")

row = ctk.CTkFrame(rmsd_frame)
row.pack(fill = "x", pady = 4)

ctk.CTkLabel(
    row,
    text = "Set a penalty for non-overlapping residues in RMSD calculations.\nPenalty = Mean RMSD + (penalty * SD)",
    justify = "left",
    anchor = "w"
    ).pack(anchor = "w", padx = 6)
ctk.CTkEntry(row, textvariable = penalty_var, width = 120).pack(anchor = "w", padx = 6, pady = 4)


row2 = ctk.CTkFrame(rmsd_frame)
row2.pack(fill = "x", pady = 4)
ctk.CTkLabel(row2, 
             text = "Specify the Euclidean Distance to be used for as the cut height to determine RMSD cluster groups.\nFor the first run. 0 is reccomended as this will default to using the mean Euclidean distance.",
             justify = "left",
             anchor = "w"
             ).pack(anchor = "w", padx = 6)
ctk.CTkEntry(row2, textvariable = custom_cut_height_var, width = 80).pack(side = "left", padx = 6)


# -----------------------------------------------------------
# 5. Thermodynamics
# -----------------------------------------------------------

thermo_frame = section("5. Thermodynamics")
thermodynamics_var = ctk.IntVar(value = 1)
ctk.CTkCheckBox(thermo_frame, text = "Run thermodynamic analysis (FoldX)", variable = thermodynamics_var).pack(anchor = "w", pady = 4)

remove_poorly_resolved_var = ctk.IntVar(value = 1)
ph_var = ctk.StringVar(value = "7.4")
temp_var = ctk.StringVar(value = "298")
ion_var = ctk.StringVar(value = "0.15")

# ctk.CTkCheckBox(thermo_frame, text = "Remove individual poorly resolved residues from otherwise well resolved PDBs?", variable = remove_poorly_resolved_var).pack(anchor = "w", pady = 4)

ctk.CTkLabel(thermo_frame, 
             text = "Remove individual poorly resolved residues from otherwise well resolved PDBs?",
             anchor = "w",
             justify = "left"
             ).pack(anchor = "w", pady = 4)

remove_poorly_resolved_var = ctk.StringVar(value = "Yes")  # default choice

radio_frame = ctk.CTkFrame(thermo_frame)
radio_frame.pack(anchor = "w", pady = 2)

ctk.CTkRadioButton(radio_frame, text = "Yes", variable = remove_poorly_resolved_var, value = "Yes").pack(side = "left", padx = 6)
ctk.CTkRadioButton(radio_frame, text = "No", variable = remove_poorly_resolved_var, value = "No").pack(side = "left", padx = 6)


grid = ctk.CTkFrame(thermo_frame)
grid.pack(fill = "x", pady = 4)

ctk.CTkLabel(grid, 
             text = "Here you can control the parameters used for FoldX calculations",
             justify = "left",
             anchor = "w"
             ).pack(anchor = "w", padx = 6)

for label_text, var in [
    ("pH", ph_var),
    ("Temperature (K)", temp_var),
    ("Ionic Strength", ion_var)
]:
    ctk.CTkLabel(grid, text = label_text).pack(side = "left", padx = 6)
    ctk.CTkEntry(grid, textvariable = var, width = 80).pack(side = "left", padx = 6)

# -----------------------------------------------------------
# 6. Stable regions
# -----------------------------------------------------------

stable_frame = section("6. Stable Regions")
stable_regions_var = ctk.IntVar(value = 1)
ctk.CTkCheckBox(stable_frame, text = "Define and calculate stable regions", variable = stable_regions_var).pack(anchor = "w", pady = 4)

grid2 = ctk.CTkFrame(stable_frame)
grid2.pack(fill = "x", pady = 4)

window_var = ctk.StringVar(value = "3")
min_stable_region_size_var = ctk.StringVar(value = "0")

ctk.CTkLabel(grid2, 
             text = "Specify the window size for the sliding window rolling average used for stable region detection",
             justify = "left",
             anchor = "w"
             ).pack(side = "left", padx = 6)
ctk.CTkEntry(grid2, textvariable = window_var, width = 80).pack(side = "left", padx = 6)

ctk.CTkLabel(grid2, 
             text = "Specify the minimum number of residues required to be considered a stabilising region",
             justify = "left",
             anchor = "w"
             ).pack(side = "left", padx = 6)
ctk.CTkEntry(grid2, textvariable = min_stable_region_size_var, width = 80).pack(side = "left", padx = 6)


distance_threshold_var = ctk.StringVar(value = "10.8")
row3 = ctk.CTkFrame(stable_frame)
row3.pack(fill = "x", pady = 4)
ctk.CTkLabel(row3, text = "Distance Threshold for Stabilising Regions/Residues Contacts (Å)").pack(side = "left", padx = 6)
ctk.CTkEntry(row3, textvariable = distance_threshold_var, width = 80).pack(side = "left", padx = 6)


# -----------------------------------------------------------
# 7. Beta strand
# -----------------------------------------------------------

beta_frame = section("7. Beta Strand Analysis")
b_strand_var = ctk.IntVar(value = 1)
ctk.CTkCheckBox(beta_frame, text = "Analyse beta sheets", variable = b_strand_var).pack(anchor = "w", pady = 4)

min_length_var = ctk.StringVar(value = "0")
row4 = ctk.CTkFrame(beta_frame)
row4.pack(fill = "x", pady = 4)
ctk.CTkLabel(row4, 
             text = "Set the minimum number of residues required to be considered a beta-strand",
             justify = "left",
             anchor = "w"
             ).pack(side = "left", padx = 6)
ctk.CTkEntry(row4, textvariable = min_length_var, width = 80).pack(side = "left", padx = 6)


# -----------------------------------------------------------
# 8. PNG generation
# -----------------------------------------------------------

png_frame = section("8. PNG Generation")
png_var = ctk.IntVar(value = 1)
ctk.CTkCheckBox(png_frame, text = "Generate PNGs and figures", variable = png_var).pack(anchor = "w", pady = 4)

png_colouring_var = ctk.StringVar(value = "0")
selected_colour_var = ctk.StringVar(value = "tv_blue")
png_palette_var = ctk.StringVar(value = "tv_blue, tv_red, tv_green, tv_orange")

row5 = ctk.CTkFrame(png_frame)
row5.pack(fill = "x", pady = 4)

ctk.CTkLabel(
    row5,
    text = "Please select a colour scheme (0–3)",
    anchor = "w",
    justify = "left"
).pack(anchor = "w", padx = 6)

bullet_frame = ctk.CTkFrame(row5)
bullet_frame.pack(anchor = "w", padx = 12)

ctk.CTkLabel(bullet_frame, text = "• 0 = Single Colour").pack(anchor = "w")
ctk.CTkLabel(bullet_frame, text = "• 1 = Protein Domain (Tau, Asyn or A-beta Only)").pack(anchor = "w")
ctk.CTkLabel(bullet_frame, text = "• 2 = Protofilament").pack(anchor = "w")
ctk.CTkLabel(bullet_frame, text = "• 3 = Amyloid Fold").pack(anchor = "w")

ctk.CTkEntry(row5, textvariable = png_colouring_var, width = 80).pack(anchor = "w", padx = 6, pady = 4)

row6 = ctk.CTkFrame(png_frame)
row6.pack(fill = "x", pady = 4)
ctk.CTkLabel(row6, 
             text = "Please provide a comma seperated colour palette for Pymol:",
             justify = "left",
             anchor = "w"
             ).pack(anchor = "w", padx = 6)
ctk.CTkEntry(row6, textvariable = png_palette_var, width = 300).pack(anchor = "w", padx = 6, pady = 4)

# -----------------------------------------------------------
# Run button
# -----------------------------------------------------------

run_button = ctk.CTkButton(scroll, text = "Run Selected Jobs", command = run_pipeline, height = 40, corner_radius = 8)
run_button.pack(pady = 10)


# -----------------------------------------------------------
# Main loop
# -----------------------------------------------------------

app.mainloop()
