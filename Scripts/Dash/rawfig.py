import numpy as np

def convert_npz_to_txt(npz_file, output_dir):
    """
    Reads a .npz file and converts its contents into .txt files.

    Parameters:
        npz_file (str): Path to the input .npz file.
        output_dir (str): Directory to save the output .txt files.
    """
    # Load the .npz file
    data = np.load(npz_file)

    # Iterate through each key in the .npz file
    for key in data.files:
        # Construct the output .txt file path
        output_file = f"{output_dir}/{key}.txt"

        # Save the array to a .txt file
        np.savetxt(output_file, data[key], fmt='%g')

        print(f"Saved {key} to {output_file}")

# Example usage
npz_file = "C:\\Xcalibur\\data\\TD\\Mike Sussman\\20250409 Native\\20240409 TD Native SEC Peak 1_CDMS_2025-04-09-08-29-05_unidecfiles\\20240409 TD Native SEC Peak 1_CDMS_2025-04-09-08-29-05_rawdata.npz"  # Replace with your .npz file path
output_dir = "C:\\Xcalibur\\data\\TD\\Mike Sussman\\20250409 Native\\20240409 TD Native SEC Peak 1_CDMS_2025-04-09-08-29-05_unidecfiles"  # Replace with your desired output directory
convert_npz_to_txt(npz_file, output_dir)