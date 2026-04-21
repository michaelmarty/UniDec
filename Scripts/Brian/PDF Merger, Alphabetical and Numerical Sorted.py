from PyPDF2 import PdfMerger
from pathlib import Path
import re

# Ask user for folder path
folder_input = input("Enter the folder path containing the PDFs: ").strip()

# Remove surrounding quotes if present (e.g., "C:\path\to\folder")
if (folder_input.startswith('"') and folder_input.endswith('"')) or \
   (folder_input.startswith("'") and folder_input.endswith("'")):
    folder_input = folder_input[1:-1]

pdf_folder = Path(folder_input).expanduser().resolve()

if not pdf_folder.is_dir():
    raise NotADirectoryError(f"Not a valid directory: {pdf_folder}")

# Output PDF in the same folder
output_pdf = pdf_folder / "merged.pdf"

def sort_key(path: Path):
    """
    Split filename into (text_part, numeric_part) for sorting.
    Examples:
      '1.pdf'        -> ('', 1)
      '2.pdf'        -> ('', 2)
      '10.pdf'       -> ('', 10)
      'A1.pdf'       -> ('A', 1)
      'Chapter2.pdf' -> ('Chapter', 2)
    """
    name = path.stem  # filename without extension
    match = re.match(r"^(.*?)(\d+)$", name)
    if match:
        text_part, num_part = match.groups()
        return (text_part, int(num_part))
    else:
        # No trailing number; sort by full name, numeric part = 0
        return (name, 0)

# Get all PDFs in the folder, sorted with the custom key
pdf_files = sorted(pdf_folder.glob("*.pdf"), key=sort_key)

if not pdf_files:
    raise FileNotFoundError(f"No PDF files found in {pdf_folder}")

merger = PdfMerger()
for pdf_path in pdf_files:
    merger.append(str(pdf_path))

merger.write(str(output_pdf))
merger.close()

print(f"Merged PDF saved to: {output_pdf}")