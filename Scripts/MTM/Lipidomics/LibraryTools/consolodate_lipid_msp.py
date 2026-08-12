import os
import numpy as np
import PublicScripts.Lipids.LipidFunctions as lf


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


def get_new_name(name):
    # classname = name.split(" ", 1)[0]
    # tname, tl, tu = lf.count_tails(name, verbose=True)
    # newtname = str(tl) + ":" + str(tu)
    # prefix = ""
    # postfix = ""
    # if "O-" in name:
    #     prefix = "O-"
    # if ";" in name:
    #     print(name)
    newname = lf.make_sum_comp(name)

    return newname, newname.split(" ", 1)[1]


def convert_to_sumcomp(input_path, output_path, drop_old_entries=False, write_output=True,
                       change_classes=["LPC", "PC", "SM"], verbose=False):
    with open(input_path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    entries = split_msp_entries(lines)
    modified_entries = []
    non_pc_entries = []
    old_pc_entries = []
    sum_comps = []

    for entry in entries:
        compound_class = next((e for e in entry if e.startswith("COMPOUNDCLASS:")), None)
        compound_class = compound_class.split(":", 1)[1].strip() if compound_class else "Unknown"
        # print("Class:", compound_class)
        if compound_class in change_classes:
            old_pc_entries.append(entry)
            # Get NAME:
            name_line = next((e for e in entry if e.startswith("NAME:")), None)
            if name_line:
                name = name_line.split(":", 1)[1].strip()
                # Get Sum Comps
                newname, newtname = get_new_name(name)

                if newtname in sum_comps:
                    continue
                else:
                    sum_comps.append(newtname)
                    if verbose:
                        print(f"Renamed: {name} -> {newname}")

                    # Replace NAME line
                    new_entry = []
                    for line in entry:
                        if line.startswith("NAME:"):
                            new_entry.append(f"NAME: {newname}\n")
                        else:
                            new_entry.append(line)
                    modified_entries.append(new_entry)

        else:
            non_pc_entries.append(entry)

    print(f"Processed {len(entries)} entries.")
    print(f"\tModified PC entries: {len(modified_entries)}")
    print(f"\tNon-PC entries: {len(non_pc_entries)}")
    print(f"\tOld PC entries: {len(old_pc_entries)}")

    if write_output:
        if not drop_old_entries:
            modified_entries = old_pc_entries + modified_entries
        modified_entries = non_pc_entries + modified_entries
        print(f"\t Output entries: {len(modified_entries)}")

        with open(output_path, "w", encoding="utf-8") as f:
            for i, entry in enumerate(modified_entries):
                f.writelines(entry)
                if i < len(modified_entries) - 1:
                    f.write("\n")


        print(f"Wrote modified MSP: {output_path}")


if __name__ == "__main__":
    os.chdir(r"Z:\Group Share\Annika\DDA_Compare")
    in_file = "MSDIAL-TandemMassSpectralAtlas-VS69-Neg-RTP2_Fixed.msp"

    in_file = "MSDIAL-TandemMassSpectralAtlas-VS69-Pos-StellarRTP2-Apex2_no_na.msp"
    in_file = "MSDIAL-TandemMassSpectralAtlas-VS69-Pos-StellarRTP2-IQX_AqX_no_na.msp"
    in_file = "MSDIAL-TandemMassSpectralAtlas-VS69-Pos-StellarRTP2_no_na.msp"
    # in_file = "IQX_neg_OTOT_RTP2.msp"
    # in_file = "Apex_neg_OTOT_RTP2.msp"
    # out_file = os.path.splitext(in_file)[0] + "_no_hcoo.msp"
    # filter_msp_entries(in_file, out_file)

    out_file = os.path.splitext(in_file)[0] + "_sumcomp.msp"
    convert_to_sumcomp(in_file, out_file, drop_old_entries=True, write_output=True)
