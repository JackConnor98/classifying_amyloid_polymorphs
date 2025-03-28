import os
from PIL import Image, ImageDraw, ImageFont
import math

# Directory containing PNGs
png_dir = "Output/PDBs/asymetric_unit/pngs"

# Get list of PNG files
png_files = [file for file in os.listdir(png_dir) if file.endswith(".png")]
png_files.sort()  # Sort to maintain order

# Extract the part before the first underscore
pdb_names = [file.split('_', 1)[0] for file in png_files]

# Load all images
images = [Image.open(os.path.join(png_dir, file)) for file in png_files]

# Determine the size of each image
img_width, img_height = images[0].size

# Load a font
# Adjust the path to the font as needed, or use a default PIL font
font_path = "/mnt/c/Windows/Fonts/Arial.ttf"
font_size = 120 #70

try:
    font = ImageFont.truetype(font_path, font_size)
except IOError:
    font = ImageFont.load_default()

# Calculate space needed for text
sample_text = "Sample"
text_bbox = font.getbbox(sample_text)
text_height = text_bbox[3] - text_bbox[1] + 100  # Add a bit of padding

# Number of images per row
images_per_row = math.isqrt(len(images))
print(f"Images Per Row : {images_per_row}")

# Calculate number of rows needed
num_rows = (len(images) + images_per_row - 1) // images_per_row

# Determine the size of the combined image (including space for text)
total_width = img_width * images_per_row
total_height = int(((img_height + text_height) * num_rows) + (2 * text_height)) # text_height buffer space at start and end

# Create a new blank image with the appropriate size
combined_image = Image.new("RGB", (total_width, total_height), (255, 255, 255)) # Use RGBA for transparency

# Draw object to add text
draw = ImageDraw.Draw(combined_image)

# Paste each image into the combined image in a grid and add text
x_offset = 0
y_offset = text_height
for index, img in enumerate(images):
    combined_image.paste(img, (x_offset, y_offset))
    # Calculate position for text
    text_x = x_offset + (img_width // 2)
    text_y = y_offset + img_height + int(text_height/2)
    # Draw the text
    draw.text((text_x, text_y), pdb_names[index], font=font, fill="black", anchor="mm")
    x_offset += img_width
    if (index + 1) % images_per_row == 0:
        x_offset = 0
        y_offset += img_height + text_height

# Save the combined image
combined_image.save(os.path.join(png_dir, "combined_grid.jpg"))
