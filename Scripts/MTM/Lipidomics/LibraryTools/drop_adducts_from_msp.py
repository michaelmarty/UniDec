import os


def split_msp_entries(lines):
    """Split MSP file lines into blank-line separated entry blocks."""
    entries = []
    current = []

    for line in lines:
        if line.strip() == "":
            if current:
                entries.append(current)
                current = []
            continue
        current.append(line)

    if current:
        entries.append(current)

    return entries


def filter_msp_entries(input_path, output_path, excluded_precursor="[M+HCOO]-"):
    target_line = f"PRECURSORTYPE: {excluded_precursor}"

    with open(input_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    entries = split_msp_entries(lines)
    kept_entries = [
        entry for entry in entries if not any(line.strip() == target_line for line in entry)
    ]

    with open(output_path, "w", encoding="utf-8") as f:
        for i, entry in enumerate(kept_entries):
            f.writelines(entry)
            if i < len(kept_entries) - 1:
                f.write("\n")

    removed = len(entries) - len(kept_entries)
    print(f"Input entries: {len(entries)}")
    print(f"Removed entries: {removed}")
    print(f"Retained entries: {len(kept_entries)}")
    print(f"Wrote filtered MSP: {output_path}")


if __name__ == "__main__":
    os.chdir(r"Z:\Group Share\Annika\DDA_Compare")
    in_file = "MSDIAL-TandemMassSpectralAtlas-VS69-Neg-RTP2_Fixed.msp"
    in_file = "IQX_neg_OTOT_RTP2.msp"
    in_file = "Apex_neg_OTOT_RTP2.msp"
    out_file = os.path.splitext(in_file)[0] + "_no_hcoo.msp"
    filter_msp_entries(in_file, out_file)

