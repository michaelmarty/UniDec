import re
from collections import defaultdict

chunk_file = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\HEK Feb 2026\MSP_chunks\MSP_chunk_7.msp"
output_log = r"C:\Temp\chunk7_potential_issues.txt"


def parse_entries(file_path):
    with open(file_path, "r", errors="ignore") as f:
        content = f.read()
    entries = content.split("NAME:")
    return ["NAME:" + e for e in entries if e.strip()]


def inspect_entry(entry, seen_mz):
    issues = []

    # RETENTIONTIME
    rt_line = next((l for l in entry.split("\n") if l.startswith("RETENTIONTIME:")), None)
    if not rt_line:
        issues.append("Missing RETENTIONTIME")
    else:
        try:
            rt = float(rt_line.split(":")[1].strip())
            if rt <= 0:
                issues.append(f"Non-positive RETENTIONTIME: {rt}")
        except:
            issues.append("Non-numeric RETENTIONTIME")

    # PRECURSORMZ duplicates
    mz_line = next((l for l in entry.split("\n") if l.startswith("PRECURSORMZ:")), None)
    if mz_line:
        try:
            mz = float(mz_line.split(":")[1].strip())
            if mz in seen_mz:
                issues.append(f"Duplicate PRECURSORMZ: {mz}")
            else:
                seen_mz.add(mz)
        except:
            issues.append("Non-numeric PRECURSORMZ")
    else:
        issues.append("Missing PRECURSORMZ")

    # SMILES / NAME length
    smiles_line = next((l for l in entry.split("\n") if l.startswith("SMILES:")), "")
    name_line = next((l for l in entry.split("\n") if l.startswith("NAME:")), "")
    if len(smiles_line) > 1000:
        issues.append(f"Very long SMILES ({len(smiles_line)} chars)")
    if len(name_line) > 200:
        issues.append(f"Very long NAME ({len(name_line)} chars)")

    # Num Peaks mismatch
    num_line = next((l for l in entry.split("\n") if l.startswith("Num Peaks:")), None)
    if num_line:
        try:
            num_peaks = int(num_line.split(":")[1].strip())
            peak_lines = [l for l in entry.split("Num Peaks:")[-1].strip().split("\n")[1:] if l.strip()]
            if len(peak_lines) != num_peaks:
                issues.append(f"Num Peaks mismatch ({num_peaks} vs {len(peak_lines)})")
        except:
            issues.append("Invalid Num Peaks format")
    else:
        issues.append("Missing Num Peaks")

    return issues


entries = parse_entries(chunk_file)
seen_mz = set()
flagged_entries = []
import re

# ===== SETTINGS =====
input_file = r"Z:\Group Share\Annika\Stellar\Untargeted DDA\HEK Feb 2026\MSP_chunks\MSP_chunk_7.msp"
output_file = r"C:\Temp\MSP_chunk7_cleaned_for_MSDIAL.msp"

max_name_len = 200
max_smiles_len = 1000


# ===== FUNCTIONS =====
def parse_entries(file_path):
    with open(file_path, "r", errors="ignore") as f:
        content = f.read()
    entries = content.split("NAME:")
    return ["NAME:" + e for e in entries if e.strip()]


def sanitize_entry(entry, seen_mz):
    lines = entry.split("\n")
    new_lines = []

    for i, line in enumerate(lines):
        # ===== NAME / SMILES / INCHIKEY sanitization =====
        if line.startswith(("NAME:", "SMILES:", "INCHIKEY:")):
            key, val = line.split(":", 1)
            val = val.encode("ascii", errors="ignore").decode()
            val = val.replace("\\", "/")
            val = val.strip()
            if key == "NAME" and len(val) > max_name_len:
                val = val[:max_name_len]
            if key == "SMILES" and len(val) > max_smiles_len:
                val = val[:max_smiles_len]
            new_lines.append(f"{key}:{val}")

        # ===== RETENTIONTIME =====
        elif line.startswith("RETENTIONTIME:"):
            try:
                rt = float(line.split(":")[1].strip())
                if rt <= 0:
                    rt = 0.1
            except:
                rt = 0.1
            new_lines.append(f"RETENTIONTIME: {rt}")

        # ===== PRECURSORMZ =====
        elif line.startswith("PRECURSORMZ:"):
            try:
                mz = float(line.split(":")[1].strip())
                while mz in seen_mz:
                    mz += 0.0001  # offset duplicates slightly
                seen_mz.add(mz)
                new_lines.append(f"PRECURSORMZ: {mz}")
            except:
                new_lines.append(line)  # leave as-is if parsing fails

        # ===== Num Peaks =====
        elif line.startswith("Num Peaks:"):
            try:
                num_peaks = int(line.split(":")[1].strip())
                peak_lines = [l for l in lines[i + 1:] if l.strip()]
                if len(peak_lines) != num_peaks:
                    num_peaks = len(peak_lines)
                new_lines.append(f"Num Peaks: {num_peaks}")
            except:
                new_lines.append(line)
        else:
            new_lines.append(line)

    return "\n".join(new_lines)


# ===== PROCESS =====
entries = parse_entries(input_file)
seen_mz = set()
cleaned_entries = [sanitize_entry(e, seen_mz) for e in entries]

with open(output_file, "w", encoding="utf-8") as f:
    f.write("\n\n".join(cleaned_entries))

print(f"Processed {len(entries)} entries.")
print(f"Cleaned MSP written to {output_file}")
for e in entries:
    issues = inspect_entry(e, seen_mz)
    if issues:
        flagged_entries.append((e, issues))

# Write flagged entries
with open(output_log, "w", encoding="utf-8") as f:
    for entry, issues in flagged_entries:
        f.write(f"# Issues: {', '.join(issues)}\n{entry}\n\n")

print(f"Scanned {len(entries)} entries in Chunk 7.")
print(f"Flagged {len(flagged_entries)} potential problem entries written to {output_log}")