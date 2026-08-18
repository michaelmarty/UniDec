# Define file paths
file1 = r"Z:\Group Share\BHT\Q Exactive HF Data\RPLC-MS\Agilent LC\20251113\UniDec_Figures_and_Files\20251113_BHT_10XDiluted_ecolilysate_C4_ensemble_3_peaks.txt"
file2 = r"Z:\Group Share\BHT\Q Exactive HF Data\RPLC-MS\Agilent LC\20251113\20251114_BHT_10XDiluted_ecolilysate_C4_12INJ_CDMS_1_unidecfiles\20251114_BHT_10XDiluted_ecolilysate_C4_12INJ_CDMS_1_peaks.dat"

# Read the first file and extract the first column
with open(file1, 'r') as f1:
    data1 = [line.split(sep=',')[0] for line in f1 if line.strip()]  # Split by ',' and take the first column

# Read the second file and extract the first column
with open(file2, 'r') as f2:
    data2 = [line.split()[0] for line in f2 if line.strip()]  # Split by whitespace and take the first column

# Convert the first dataset to integers
data1 = [int(float(value)) for value in data1]

# Convert the second dataset to integers
data2 = [int(float(value)) for value in data2]

# Print the first few values of each file (optional)
print("First file first column preview:")
print(data1[:5])  # Print first 5 values of the first column

print("\nSecond file first column preview:")
print(data2[:5])  # Print first 5 values of the first column

# Find numbers within 10 of each other between the datasets
matches = [(val1, val2) for val1 in data1 for val2 in data2 if abs(val1 - val2) <= 10]

# Print the matches
print("Numbers within 10 of each other:")
for match in matches:
    print(match)