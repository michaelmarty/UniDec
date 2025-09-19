import pandas as pd
import numpy as np
import pprint


EXCLUSION_TIME_DIFFERENCE = 0.25

# Lab file location.
# excel_file = 'Z:Group Share\\James Rohrbough\\FPOP\\MSMS raw data\\ubiquitin\\FPOP peptides\\JGR_2319383_2Ox-calib_Peptides.xlsx'
# Local file location.

excel_file = "output.xlsx"
peptides = pd.read_excel(excel_file)

# sorted_by_protein = peptides.sort(['Protein Name'])
# sorted_by_RT = peptides.sort.value(['Scan Retention'])

# Grab filtered data.
filtered_data = peptides.loc[peptides["Protein Name"] == "Polyubiquitin-B|"]
filtered_data = filtered_data.fillna("")

# Summarize results.
print("{} peptides extracted from a total of {}".format(
    len(filtered_data), len(peptides)
))

# Store the counts for each base sequence. Maps sequence to count.
base_sequence_counts = dict()

for base_sequence in filtered_data["Base Sequence"]:
    # Create dictionary entry if it doesn't exist.
    base_sequence_counts.setdefault(base_sequence, 0)
    # Update the count.
    base_sequence_counts[base_sequence] += 1

# Maps the base_sequence to the retention time.
reference_retention_times = dict()
for key in base_sequence_counts:
    # Get group of matching base sequences
    group = peptides.loc[peptides["Base Sequence"] == key]
    total_sum = 0

    # Find the retention times for each base sequence.
    for index, row in group.iterrows():
        # Append unmodified entries to calculate reference retention time.
        if pd.isnull(row["Mods"]):
            reference_retention_times.setdefault(key, list())
            reference_retention_times[key].append(row["Scan Retention Time"])

    # Error Handling: When a reference time couldn't be found.
    if reference_retention_times.get(key) is None:
        print(
            "EXCLUSION: Could not find base sequence '{}' without"
            " mods to use as reference time.".format(key))
        continue

    # Reference retention time = mean(unmodified_entries)
    reference_retention_times[key] = np.mean(reference_retention_times[key])

    # Notify user of retention time.
    print("'{}' reference retention time: {}".format(
        key, reference_retention_times[key]))

    oxidized_data = dict()

    for index, row in group.iterrows():
        if (
            row["Scan Retention Time"] < 
            (reference_retention_times[key] - EXCLUSION_TIME_DIFFERENCE) or
            row["Scan Retention Time"] > 
            (reference_retention_times[key] + EXCLUSION_TIME_DIFFERENCE)
        ):
            essential_sequence = row["Essential Sequence"]
            parsed_sequence = essential_sequence

            # Grab start residue for current peptide.
            frag_len = int(
                row["Start and End Residues In Protein"].split(" ")[0][1:]) - 1

            # Include modified fragments only.
            if pd.isnull(row["Mods"]) is False:
                # Find oxidation based on brackets.
                while "]" in parsed_sequence:
                    total_ion_current = row["Total Ion Current"]
                    total_sum += total_ion_current
                    parsed_sequence = parsed_sequence.split("[", 1)
                    frag_len += len(parsed_sequence[0])
                    parsed_sequence = parsed_sequence[1].split("]", 1)
                    parsed_sequence = parsed_sequence[1]

                    oxidized_data.setdefault(frag_len, 0)
                    oxidized_data[frag_len] += float(total_ion_current)

    for key, item in oxidized_data.items():
        oxidized_data[key] = item / total_sum
    print(pprint.pformat(oxidized_data))
