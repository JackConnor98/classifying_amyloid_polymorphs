import pymol
from pymol import cmd
import os
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
import math
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

def color_white():
    cmd.color("white")


def validate_color(color):
    valid_colors = {name for name, _ in cmd.get_color_indices()}
    if color not in valid_colors:
        print(f"Warning: '{color}' is not a recognised PyMOL color. Reverting to 'white'.")
        return "white"
    return color


def set_coloring(png_coloring, user_palette, metadata=None, fibril_info=None, stable_regions=None):
    """Return a `color_residues(object_name, pdb_name)` callable based on settings.

    - `metadata` and `fibril_info` should be pandas DataFrames when required.
    - `stable_regions` path is required for png_coloring==4.
    - The returned function should be called with the current `object_name` and `pdb_name`.
    """

    # Normalize palette
    user_palette = user_palette or []

    # Helper to return a no-op colorer
    def _noop(object_name, pdb_name):
        color_white()

    # png_coloring == 0: single user color
    if png_coloring == 0:
        if not user_palette:
            selected_color = "white"
        else:
            selected_color = validate_color(user_palette[0])

        def _color_single(object_name, pdb_name):
            cmd.color(selected_color, object_name)

        return _color_single

    # png_coloring == 1: domain-based coloring driven by `metadata`
    if png_coloring == 1 and metadata is not None:
        protein_column = metadata["Protein"].str.lower()
        asyn = protein_column.str.contains("-synuclein").any()
        tau = protein_column.str.contains("tau").any()
        abeta = protein_column.str.contains("amyloid-").any()

        # Asyn Colouring
        if asyn:
            def _color_domains(object_name, pdb_name):
                # N-terminus
                cmd.color("tv_blue", f"{object_name} and resi 1-60")
                # NAC
                cmd.color("tv_red", f"{object_name} and resi 61-95")
                # C-terminus
                cmd.color("tv_green", f"{object_name} and resi 96-140")

            return _color_domains

        # Tau Colouring
        if tau:
            def _color_domains(object_name, pdb_name):
                # R1
                cmd.color("tv_red", f"{object_name} and resi 243-273")
                # R2
                cmd.color("tv_green", f"{object_name} and resi 274-304")
                # R3
                cmd.color("tv_orange", f"{object_name} and resi 305-335")
                # R4
                cmd.color("tv_blue", f"{object_name} and resi 336-367")

            return _color_domains

        # A-beta Colouring
        if abeta:
            def _color_domains(object_name, pdb_name):
                # N-term
                cmd.color("tv_blue", f"{object_name} and resi 1-16")
                # Central Hydrophobic Core
                cmd.color("black", f"{object_name} and resi 17-21")
                # Turn Region
                cmd.color("tv_red", f"{object_name} and resi 24-27")
                # Second Hydrophobic Region
                cmd.color("tv_orange", f"{object_name} and resi 28-35")
                # C-term
                cmd.color("tv_green", f"{object_name} and resi 36-42")

            return _color_domains

        print("Domain coloring is not supported for this protein — falling back to white.")
        return _noop

    # png_coloring in [2,3]: group-based coloring using fibril_info
    if png_coloring in [2, 3] and fibril_info is not None:
        target_col = "fibril" if png_coloring == 2 else "polymorph"
        fibril_groups = fibril_info[target_col].unique()

        default_colors = [
            "tv_red", "tv_blue", "tv_green", "tv_orange", "yellow",
            "purple", "cyan", "salmon", "lime", "deepblue"
        ]

        validated_palette = [validate_color(c) for c in user_palette]
        filtered_defaults = [c for c in default_colors if c not in validated_palette]
        pymol_colors = validated_palette + filtered_defaults

        color_map = {grp: pymol_colors[i % len(pymol_colors)] for i, grp in enumerate(fibril_groups)}

        def _color_groups(object_name, pdb_name):
            curr = fibril_info[fibril_info["PDB"] == pdb_name]
            if curr.empty:
                print(f"No fibril_info entry found for {pdb_name}, nothing coloured.")
                return
            for _, row in curr.iterrows():
                chain = row["chain"]
                group_val = row[target_col]
                colour = color_map.get(group_val, "white")
                cmd.color(colour, f"{object_name} and chain {chain}")

        return _color_groups

    # png_coloring == 4: stable region-based coloring
    if png_coloring == 4 and stable_regions is not None:
        
        df = pd.read_csv(stable_regions)

        if "residues" not in df.columns:
            raise ValueError(
                "stable_regions CSV must contain a 'residues' column"
            )

        stable_regions = df["residues"].tolist()

        # Load the qualitative palettes
        set1 = plt.get_cmap("Set1")
        dark2 = plt.get_cmap("Dark2")

        # Convert to hex
        set1_colours = [mcolors.to_hex(set1(i)) for i in range(set1.N)]
        dark2_colours = [mcolors.to_hex(dark2(i)) for i in range(dark2.N)]

        # Combine them
        base_colours = set1_colours + dark2_colours

        # Expand to match number of stable regions
        region_colours = np.tile(base_colours, int(np.ceil(len(stable_regions) / len(base_colours))))[:len(stable_regions)]

        # Register colours in PyMOL
        for idx, hex_colour in enumerate(region_colours):
            rgb = mcolors.to_rgb(hex_colour)
            cmd.set_color(f"region_colour_{idx}", rgb)

        def _color_stable_regions(object_name, pdb_name):
            # Iterate over each region and color it
            for j, row in df.iterrows():
                region_range = row["residues"]  # Getting the residues in the current stable region
                cmd.color(f"region_colour_{j}", f"{object_name} and resi {region_range}")

        return _color_stable_regions

    # Fallback
    return _noop


def process_pdb_directory(pdb_path, save_path, color_residues, apply_transparency=False, fibril_info=None, reference_pdb=0):
    """Process all PDBs in `pdb_path`, saving PNGs to `save_path`.

    Args:
        pdb_path: Directory containing PDB files
        save_path: Directory to save PNG images to
        color_residues: Callable accepting (object_name, pdb_name) to apply coloring
        apply_transparency: If True, apply transparency to non-matching chains (requires `fibril_info`)
        fibril_info: DataFrame with chain mapping information (required if apply_transparency=True)
        reference_pdb: Index into `pdb_files` to select the reference structure (default 0)
    """
    if not os.path.exists(save_path):
        os.makedirs(save_path, exist_ok=True)

    pdb_files = [file for file in os.listdir(pdb_path) if file.endswith(".pdb")]

    # If there are no PDB files, nothing to do
    if not pdb_files:
        return None

    # Determine reference structure path from available PDB files.
    # `reference_pdb` is an index into the `pdb_files` list (default 0).
    try:
        ref_idx = int(reference_pdb)
    except Exception:
        ref_idx = 0
    if ref_idx < 0 or ref_idx >= len(pdb_files):
        ref_idx = 0
    reference_path = os.path.join(pdb_path, pdb_files[ref_idx])

    for i in pdb_files:
        current_path = os.path.join(pdb_path, i)
        object_name = i.replace(".pdb", "").replace("_asym_unit", "")
        pdb_name = object_name.split("_")[0]

        cmd.load(current_path, object_name)
        cmd.show("cartoon")
        cmd.set("cartoon_loop_radius", 1)
        # Removing dashed lines between disconnected residues
        cmd.set("cartoon_gap_cutoff", 0)

        # Apply transparency logic if enabled
        if apply_transparency and fibril_info is not None:
            # Get the fibril information for the current structure
            current_fibril_info = fibril_info[fibril_info['PDB'].str.contains(pdb_name, na=False)]
            unique_pdb_ids = current_fibril_info['pdb_id'].unique().tolist()

            # Get list of chains in the loaded structure
            chains = cmd.get_chains(object_name)

            if len(unique_pdb_ids) > 1:
                # Multi-chain case: process each pdb_id separately
                for pdb_id in unique_pdb_ids:
                    
                    # Clear previous object
                    cmd.delete(object_name) 
                    
                    pdbid_name = pdb_id
                    cmd.load(current_path, pdbid_name)

                    # Filter fibril_info to get the chains for the current pdb_id
                    expected_chains = current_fibril_info[current_fibril_info['pdb_id'] == pdb_id]['chain'].unique().tolist()

                    # Get list of chains in the loaded structure
                    chains = cmd.get_chains(pdbid_name)

                    # Identify matching and non-matching chains
                    non_matching_chains = [chain for chain in chains if chain not in expected_chains]

                    # Set non-matching chains to 50% transparency
                    for chain in non_matching_chains:
                        cmd.set("cartoon_transparency", 0.75, f"{pdbid_name} and chain {chain}")

                    # Load the reference structure and align
                    cmd.load(reference_path, "reference")
                    cmd.align(f"{pdbid_name} and name CA", "reference and name CA")

                    # Set all residues to white before coloring
                    color_white()

                    # Color specific residues
                    color_residues(pdbid_name, pdb_name)

                    # Set the background color to white
                    cmd.bg_color("white")
                    cmd.set("ray_opaque_background", 0)

                    # Orienting camera to focus on the reference so it is consistent
                    cmd.orient(pdbid_name)

                    # Deleting the reference structure
                    cmd.delete("reference")

                    # Refresh makes sure everything is updated before saving
                    cmd.refresh()

                    # Save the image as a PNG with a transparent background
                    png_path = os.path.join(save_path, f"{pdbid_name}.png")
                    cmd.png(png_path, dpi=300, ray=1, quiet=1)

                    # Clear the current object for the next iteration
                    cmd.delete(pdbid_name)
                    
            else:
                # Single chain case
                # Load the reference structure and align
                cmd.load(reference_path, "reference")
                cmd.align(f"{object_name} and name CA", "reference and name CA")

                # Set all residues to white before coloring
                color_white()

                # Color specific residues
                color_residues(object_name, pdb_name)

                # Set the background color to white
                cmd.bg_color("white")
                cmd.set("ray_opaque_background", 0)

                # Orienting camera to focus on the reference so it is consistent
                cmd.orient(object_name)

                # Deleting the reference structure
                cmd.delete("reference")

                # Refresh makes sure everything is updated before saving
                cmd.refresh()

                # Save the image as a PNG with a transparent background
                png_path = os.path.join(save_path, f"{object_name}.png")
                cmd.png(png_path, dpi=300, ray=1, quiet=1)

                # Clear the current object for the next iteration
                cmd.delete(object_name)
        else:
            # Standard processing without transparency
            # Set all residues to white before coloring
            color_white()

            # Color specific residues
            color_residues(object_name, pdb_name)

            # Set the background color to white
            cmd.bg_color("white")
            cmd.set("ray_opaque_background", 0)

            # Refresh makes sure everything is updated before saving
            cmd.refresh()

            # Save the image as a PNG with a transparent background
            png_path = os.path.join(save_path, f"{object_name}.png")
            cmd.png(png_path, dpi=300, ray=1, quiet=1)

            # Clear the current object for the next iteration
            cmd.delete(object_name)


def merge_pngs(png_dir):
    png_files = [file for file in os.listdir(png_dir) if file.endswith(".png")]
    png_files.sort()
    pdb_names = [file.replace(".png", "").replace("_asym_unit", "") for file in png_files]
    images = [Image.open(os.path.join(png_dir, file)) for file in png_files]
    img_width, img_height = images[0].size

    try:
        # Try common Windows path first
        font_path = "C:/Windows/Fonts/Arial.ttf"
        font_size = 120
        font = ImageFont.truetype(font_path, font_size)
    except Exception:
        font = ImageFont.load_default()

    sample_text = "Sample"
    text_bbox = font.getbbox(sample_text)
    text_height = text_bbox[3] - text_bbox[1] + 100

    images_per_row = math.isqrt(len(images)) or 1
    num_rows = (len(images) + images_per_row - 1) // images_per_row
    total_width = img_width * images_per_row
    total_height = int(((img_height + text_height) * num_rows) + (2 * text_height))

    combined_image = Image.new("RGB", (total_width, total_height), (255, 255, 255))
    draw = ImageDraw.Draw(combined_image)

    x_offset = 0
    y_offset = text_height
    for index, img in enumerate(images):
        combined_image.paste(img, (x_offset, y_offset))
        text_x = x_offset + (img_width // 2)
        text_y = y_offset + img_height + int(text_height / 2)
        draw.text((text_x, text_y), pdb_names[index], font=font, fill="black", anchor="mm")
        x_offset += img_width
        if (index + 1) % images_per_row == 0:
            x_offset = 0
            y_offset += img_height + text_height

    combined_image.save(os.path.join(png_dir, "combined_grid.jpg"))


def create_grouped_png_grid(
    pdb_path,
    png_dir,
    cluster_path,
    output_name = "combined_grid.png",
    transparent = False,
    font_path = "/mnt/c/Windows/Fonts/Arial.ttf",
    font_size = 140,
    title_size = 140,
    line_thickness = 10,
    top_padding = 50
):
    # ------------------------------------------------------------------
    
    # Getting RMSD cluster groups
    cluster_df = pd.read_csv(cluster_path)


    # Get PDB codes from PDB files
    pdb_files = [f for f in os.listdir(pdb_path) if f.endswith(".pdb")]
    pdb_codes = sorted({f[:4] for f in pdb_files})

    png_files = [f for f in os.listdir(png_dir) if f.endswith(".png")]

    matched_pngs = [
    f for f in png_files
    if f[:4] in pdb_codes]

    # Load images
    images = [Image.open(os.path.join(png_dir, f)) for f in matched_pngs]

    # Map pdb_name -> image
    png_names = [f.replace(".png", "").replace("_asym_unit", "") for f in matched_pngs]
    
    pdb_to_img = dict(zip(png_names, images))

    # Order according to dataframe
    ordered_pdbs = cluster_df["pdb_id"].tolist()
    images_sorted = [pdb_to_img[pdb] for pdb in ordered_pdbs]

    # ------------------------------------------------------------------
    # Fonts
    try:
        font = ImageFont.truetype(font_path, font_size)
        title_font = ImageFont.truetype(font_path, title_size)
    except IOError:
        font = ImageFont.load_default()
        title_font = font

    text_bbox = font.getbbox("Sample")
    text_height = text_bbox[3] - text_bbox[1]

    title_bbox = title_font.getbbox("Group X")
    title_height = title_bbox[3] - title_bbox[1]

    # ------------------------------------------------------------------
    # Layout
    img_width, img_height = images_sorted[0].size
    images_per_row = round(math.sqrt(len(images_sorted))) + 2

    num_images = len(cluster_df)
    rows_per_group = (
    cluster_df
    .groupby("group")
    .size()
    .apply(lambda n: math.ceil(n / images_per_row))
    )


    row_height = img_height + text_height + top_padding
    
    group_header_height = (title_height + text_height + line_thickness + text_height)

    total_height = 0

    for n_rows in rows_per_group:
        total_height += group_header_height
        total_height += n_rows * row_height

    total_width = img_width * images_per_row

    mode = "RGBA" if transparent else "RGB"
    bg_colour = (255, 255, 255, 0) if transparent else (255, 255, 255)

    combined_image = Image.new(mode, (total_width, total_height), bg_colour)
    draw = ImageDraw.Draw(combined_image)

    # ------------------------------------------------------------------
    # Drawing
    x_offset = 0
    y_offset = 0
    current_group = None

    for idx, row in cluster_df.iterrows():
        pdb_id = row["pdb_id"]
        group = row["group"]

        if group != current_group:
            current_group = group
            y_offset += text_height

            draw.text(
                (total_width // 2, y_offset),
                f"Group {group}",
                font = title_font,
                fill = "black",
                anchor = "mm",
                stroke_width = 2
            )

            y_offset += text_height
            draw.line(
                (0, y_offset, total_width, y_offset),
                fill = "black",
                width = line_thickness
            )
            y_offset += text_height
            x_offset = 0

        img = pdb_to_img[pdb_id]
        combined_image.paste(img, (x_offset, y_offset))

        draw.text(
            (x_offset + img_width // 2, y_offset + img_height + top_padding),
            pdb_id,
            font = font,
            fill = "black",
            anchor = "mm"
        )

        x_offset += img_width
        if x_offset >= total_width:
            x_offset = 0
            y_offset += img_height + text_height

        if idx < num_images - 1 and cluster_df.iloc[idx + 1]["group"] != group:
            if x_offset > 0:
                y_offset += img_height + text_height
            x_offset = 0

    # ------------------------------------------------------------------
    combined_image.save(os.path.join(png_dir, output_name))

def create_png_grid(
    pdb_path,
    png_dir,
    output_name = "combined_grid.png",
    transparent = False,
    font_path = "/mnt/c/Windows/Fonts/Arial.ttf",
    font_size = 140
):
    # ------------------------------------------------------------------
    # Get PDB codes from PDB files
    pdb_files = [f for f in os.listdir(pdb_path) if f.endswith(".pdb")]
    pdb_codes = sorted({f[:4] for f in pdb_files})

    png_files = [f for f in os.listdir(png_dir) if f.endswith(".png")]

    matched_pngs = [
    f for f in png_files
    if f[:4] in pdb_codes
]

    # Load images
    images = [Image.open(os.path.join(png_dir, f)) for f in matched_pngs]

    # Map pdb_name -> image
    png_names = [f.replace(".png", "").replace("_asym_unit", "") for f in matched_pngs]

    # ------------------------------------------------------------------
    # Fonts
    try:
        font = ImageFont.truetype(font_path, font_size)
    except IOError:
        font = ImageFont.load_default()

    text_bbox = font.getbbox("Sample")
    text_height = text_bbox[3] - text_bbox[1]
    text_padding = int(0.25 * text_height)

    # ------------------------------------------------------------------
    # Layout
    img_width, img_height = images[0].size
    images_per_row = round(math.sqrt(len(images)))

    num_images = len(images)
    num_rows = math.ceil(num_images / images_per_row)

    total_width = img_width * images_per_row
    total_height = num_rows * (img_height + text_height + text_padding)

    mode = "RGBA" if transparent else "RGB"
    bg_colour = (255, 255, 255, 0) if transparent else (255, 255, 255)

    combined_image = Image.new(mode, (int(total_width), int(total_height)), bg_colour)
    draw = ImageDraw.Draw(combined_image)

    # ------------------------------------------------------------------
    # Drawing
    x_offset = 0
    y_offset = 0
    count = 0

    for pdb_name in png_names:

        img = Image.open(os.path.join(png_dir, pdb_name + ".png"))

        combined_image.paste(img, (x_offset, y_offset))

        draw.text(
            (x_offset + img_width // 2, y_offset + img_height + text_padding),
            pdb_name,
            font = font,
            fill = "black",
            anchor = "mm"
        )

        count += 1

        x_offset += img_width
        if x_offset >= total_width:
            x_offset = 0
            y_offset += img_height + text_height

        if count == images_per_row:
            if x_offset > 0:
                y_offset += img_height + text_height
            x_offset = 0

    # ------------------------------------------------------------------
    combined_image.save(os.path.join(png_dir, output_name))


__all__ = [
    "color_white",
    "validate_color",
    "set_coloring",
    "process_pdb_directory",
    "merge_pngs",
    "create_grouped_png_grid"
]