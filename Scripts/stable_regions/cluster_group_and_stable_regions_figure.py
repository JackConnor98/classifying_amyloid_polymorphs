import os
from PIL import Image, ImageDraw, ImageFont
import pandas as pd
import math

# Getting list of RMSD cluster groups
df = pd.read_csv("Output/RMSD/RMSD_cluster_groups.csv")

# Directory containing PNGs
png_dir = "Output/stable_regions/pngs"

# Get list of PNG files
png_files = [file for file in os.listdir(png_dir) if file.endswith(".png")]

# Removing the .png to get the pdb names
pdb_names = [file.split('.', 1)[0] for file in png_files]

# Load all images
images = [Image.open(os.path.join(png_dir, file)) for file in png_files]

# Sort pdb_names based on the order in df
pdb_names_sorted = df['pdb_id'].tolist()

# Create a mapping of pdb_name to image and reorder images
pdb_to_img = dict(zip(pdb_names, images))
images_sorted = [pdb_to_img[pdb_name] for pdb_name in pdb_names_sorted]

# Load a font
font_path = "/mnt/c/Windows/Fonts/Arial.ttf"
font_size = 140
title_size = 140

try:
    font = ImageFont.truetype(font_path, font_size)
    title_font = ImageFont.truetype(font_path, title_size)
except IOError:
    font = ImageFont.load_default()

# Calculate space needed for text
sample_text = "Sample"
text_bbox = font.getbbox(sample_text)
text_height = text_bbox[3] - text_bbox[1]# + 20  # Add a bit of padding
top_padding = 25 # setting this so I can make the padding smaller at the top

# Number of images per row
images_per_row = round(math.sqrt(len(images)))# + 6
print(f"Number of images per row = {images_per_row}")

# Calculate number of rows needed, accounting for group titles
num_images = len(df)
num_groups = df['group'].nunique()
num_rows = num_images // images_per_row + num_groups


# Determine the size of the combined image (including space for text and group titles)
img_width, img_height = images[0].size
total_height = num_rows * (img_height + text_height) + (num_groups * text_height * 2)
total_width = img_width * images_per_row

# Create a new blank image with the appropriate size
combined_image = Image.new("RGB", (total_width, total_height), (255, 255, 255))  # Use RGBA for transparency

# Draw object to add text
draw = ImageDraw.Draw(combined_image)

# Line color and thickness
line_color = "black"
line_thickness = 10

# Initialize offsets
x_offset = 0
y_offset = 0
current_group = None

# Iterate through the sorted dataframe and arrange the images with titles and separators
for index, row in df.iterrows():
    pdb_id = row['pdb_id']
    group = row['group']

    # Add group title if we're in a new group
    if group != current_group:
        current_group = group
        y_offset += text_height 

        # Draw the group title centered in the row
        group_title = f"Group {group}"
        text_x = total_width // 2
        draw.text((text_x, y_offset), group_title, font=title_font, fill="black", anchor="mm", stroke_width=2)
        y_offset += text_height  # Move down after writing the title
        
        # Draw the separating line across the width of the image
        draw.line((0, y_offset, total_width, y_offset), fill=line_color, width=line_thickness)
        y_offset += text_height  # Move down after drawing the line  

        x_offset = 0  # Reset x_offset for new group

    # Paste the image in the sorted order
    img = pdb_to_img[pdb_id]
    combined_image.paste(img, (x_offset, y_offset))
    
    # Calculate position for text
    text_x = x_offset + (img_width // 2)
    #text_y = y_offset + img_height + int(text_height / 20)
    text_y = y_offset + img_height + top_padding
    
    # Draw the pdb name below the image
    draw.text((text_x, text_y), pdb_id, font=font, fill="black", anchor="mm")
    
    # Update offsets
    x_offset += img_width
    if x_offset >= total_width:  # Wrap to next line if at the end of a row
        x_offset = 0
        y_offset += img_height + text_height

    # Move to next row after the last image in a group
    if index < num_images - 1 and df.iloc[index + 1]['group'] != group:
        if x_offset > 0: 
            y_offset += img_height + text_height
        x_offset = 0


# Save the combined image
combined_image.save(os.path.join(png_dir, "combined_grid.jpg"))
