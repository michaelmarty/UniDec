import os

def fix_pg_in_msp(neglib):
    insert = "152.9958	50"
    compound_class = None

    # Read the MSP file and modify PG entries to add the missing fragment
    # Will need to do this for EtherPG and LysoPG as well
    # Need to edit num peaks line as well
    # Check if fragment already exists to avoid duplicates
    with open(neglib, 'r') as f:
        lines = f.readlines()
    with open(neglib.replace(".msp", "_Fixed.msp"), 'w') as f:
        i = 0
        while i < len(lines):
            line = lines[i]
            if line.startswith("COMPOUNDCLASS:"):
                compound_class = line.split(":", 1)[1].strip()
                f.write(line)
                i += 1
            elif line.startswith("Num Peaks:"):
                num_peaks = int(line.split(":", 1)[1].strip())
                # Check if the compound class is PG or related
                if compound_class in ["PG"]:
                    # Peek ahead to see if the fragment already exists
                    fragment_exists = False
                    for j in range(i + 1, i + 1 + num_peaks):
                        if insert.split('\t')[0] in lines[j]:
                            fragment_exists = True
                    if not fragment_exists:
                        # Increment the number of peaks
                        num_peaks += 1
                        f.write(f"Num Peaks: {num_peaks}\n")
                        # Write the existing peaks
                        for j in range(i + 1, i + 1 + num_peaks - 1):
                            mz = float(lines[j].split('\t')[0])
                            if mz > float(insert.split('\t')[0]) and not fragment_exists:
                                # Insert the missing fragment before writing this peak
                                f.write(f"{insert}\n")
                                fragment_exists = True  # Ensure we only insert once
                            f.write(lines[j])
                        i += num_peaks  # Skip the peaks we've already written
                    else:
                        f.write(f"Num Peaks: {num_peaks}\n")
                        # Just write the existing peaks without modification
                        for j in range(i + 1, i + 1 + num_peaks):
                            f.write(lines[j])
                        i += num_peaks + 1  # Skip the peaks we've already written
                else:
                    # Just write the existing peaks without modification
                    f.write(f"Num Peaks: {num_peaks}\n")
                    for j in range(i + 1, i + 1 + num_peaks):
                        f.write(lines[j])
                    i += num_peaks + 1  # Skip the peaks we've already written
            else:
                f.write(line)
                i += 1


if __name__ == "__main__":
    os.chdir(r"Z:\Group Share\Annika\Stellar\Untargeted DDA\HEK")
    neglib = "20260112_HEK_UntDDA_Neg_New_Corrected.msp"
    fix_pg_in_msp(neglib)
