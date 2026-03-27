import re

def clean_msp(input_path, output_path):
    with open(input_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    entries = []
    current_entry = []

    # Split entries by NAME field (case-insensitive)
    for line in lines:
        if line.strip() == "":
            continue
        if line.upper().startswith("NAME:") and current_entry:
            entries.append(current_entry)
            current_entry = []
        current_entry.append(line)
    if current_entry:
        entries.append(current_entry)

    print(f"Total entries detected: {len(entries)}")

    cleaned_entries = []

    for i, entry_lines in enumerate(entries):
        # Parse key: value pairs
        entry = {}
        for line in entry_lines:
            if ":" in line:
                key, val = line.split(":", 1)
                entry[key.strip().upper()] = val.strip()

        # Clean/normalize fields
        name = entry.get("NAME", f"Unknown_{i}")
        name = name.replace(" ", "_")  # replace spaces
        entry["NAME"] = name

        # Precursor type
        prec_type = entry.get("PRECURSORTYPE", "M+")
        prec_type = re.sub(r"[\[\]]", "", prec_type).upper()
        entry["PRECURSORTYPE"] = prec_type

        # Decimal precision for precursor m/z
        precursor = entry.get("PRECURSORMZ")
        if precursor:
            entry["PRECURSORMZ"] = f"{float(precursor):.5f}"

        # Retention time
        rt = entry.get("RETENTIONTIME")
        if rt:
            entry["RETENTIONTIME"] = f"{float(rt):.5f}"

        # Parse peaks
        num_peaks = entry.get("NUM PEAKS")
        peak_lines = [l for l in entry_lines if re.match(r"^\d+\.\d+\s+\d+", l.strip())]

        # Fix peak decimal precision
        cleaned_peaks = []
        for pline in peak_lines:
            mz, intensity = pline.strip().split()
            cleaned_peaks.append(f"{float(mz):.5f}  {intensity}")

        # Fix Num Peaks
        entry["NUM PEAKS"] = str(len(cleaned_peaks))

        # Build cleaned entry text
        text_lines = []
        for field in ["NAME","PRECURSORMZ","PRECURSORTYPE","SMILES","INCHIKEY",
                      "FORMULA","RETENTIONTIME","CCS","IONMODE","COMPOUNDCLASS","COMMENT","NUM PEAKS"]:
            if field in entry:
                text_lines.append(f"{field}: {entry[field]}")
        text_lines.extend(cleaned_peaks)
        text_lines.append("")  # blank line between entries

        cleaned_entries.append("\n".join(text_lines))

    # Write cleaned MSP
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write("\n".join(cleaned_entries))

    print(f"Cleaned MSP written to {output_path}")


# ---- USAGE ----
input_msp = r"Z:\Group Share\Shivam\IQX Stellar 3\MSDial\MSDIAL-TandemMassSpectralAtlas-VS69-Pos-StellarRTP2.msp"
output_msp = r"Z:\Group Share\Shivam\IQX Stellar 3\MSDial\Cleaned_MSDIAL.msp"

clean_msp(input_msp, output_msp)