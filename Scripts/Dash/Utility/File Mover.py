import os
import shutil

topdir = "Z:\\Group Share\\Dash\\Current Data\\Synapt 96 Wells"
blanks_dir = os.path.join(topdir, "Blanks")

# Create the "Blanks" directory if it doesn't exist
if not os.path.exists(blanks_dir):
    os.makedirs(blanks_dir)

# Iterate through all items in the top directory
for item in os.listdir(topdir):
    item_path = os.path.join(topdir, item)
    # Check if the item is a directory and contains the word "Blank"
    if os.path.isdir(item_path) and "Blank" in item:
        # Move the directory to the "Blanks" folder
        shutil.move(item_path, blanks_dir)
        print(f"Moved: {item_path} to {blanks_dir}")