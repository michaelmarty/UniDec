import os

base_library = r"Z:\Group Share\JGP\Lipidomics\MSDIAL-TandemMassSpectralAtlas-VS69-Neg.msp"
os.chdir(os.path.dirname(base_library))

classes_to_keep = ["PG", "PE", "CL"]

lines_to_keep = []

current_entry = []
keep_current_entry = False
with open(base_library, "r") as file:
    for line in file:
        if "NAME" in line and len(current_entry) > 0:
            if keep_current_entry:
                lines_to_keep.extend(current_entry)
                current_entry = [line]
                continue
            else:
                current_entry = [line]
                continue

        if "COMPOUNDCLASS" in line:
            class_split = line.split(": ")
            lipid_class = class_split[1].strip()
            if lipid_class not in classes_to_keep:
                keep_current_entry = False
            else:
                keep_current_entry = True

        current_entry.append(line)

with open('filtered_library.msp', 'w') as file:
    for line in lines_to_keep:
        file.write(line)




