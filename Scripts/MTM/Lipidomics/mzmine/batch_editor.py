import os
from pathlib import Path
import xml.etree.ElementTree as ET

def read_mzmine_batch_file(batchfile):
    tree = ET.parse(batchfile)
    batch = tree.getroot()
    steps = batch.findall("batchstep")
    for step in steps:
        module_name = step.attrib["module_name"]
        print(f"Step: {module_name}")
        for param in step.findall("parameter"):
            name = param.attrib["name"]
            value = param.text
            print(f"  {name}: {value}")
            if param.findall("parameter"):
                print(f"    Sub-parameters for {name}:")
                for subparam in param.findall("parameter"):
                    subname = subparam.attrib["name"]
                    subvalue = subparam.text
                    print(f"      {subname}: {subvalue}")
            for fileparam in param.findall("file"):
                filepath = fileparam.text
                print(f"      File parameter {filepath}")

    return tree

def write_mzmine_batch_file(tree, outputfile):
    tree.write(outputfile, encoding="utf-8", xml_declaration=True)
    print(f"Batch file saved to {outputfile}")

def set_spectral_library(tree, new_library_path):
    steps = tree.getroot().findall("batchstep")
    for step in steps:
        if step.attrib["module_name"] == "Import MS data":
            for param in step.findall("parameter"):
                if param.attrib["name"] == "Spectral library files":
                    for fileparam in param.findall("file"):
                        fileparam.text = new_library_path
                        print(f"Updated spectral library path to {new_library_path}")

        if step.attrib["module_name"] == "Spectral library search":
            for param in step.findall("parameter"):
                if param.attrib["name"] == "Spectral libraries":
                    for fileparam in param.findall("specific_file"):
                        fileparam.text = new_library_path
                        print(f"Updated spectral library path to {new_library_path}")

    return tree

def get_spectral_library(tree):
    steps = tree.getroot().findall("batchstep")
    for step in steps:
        if step.attrib["module_name"] == "Spectral library search":
            for param in step.findall("parameter"):
                if param.attrib["name"] == "Spectral libraries":
                    for fileparam in param.findall("specific_file"):
                        libfile = fileparam.text
                        print(f"Current spectral library path: {libfile}")
                        return libfile
    return None

def create_duplicate_copy(batchfile, copy_suffix="_copy"):
    path = Path(batchfile)
    new_path = path.with_name(path.stem + copy_suffix + path.suffix)
    new_path.write_text(path.read_text())
    print(f"Created duplicate batch file: {new_path}")
    return str(new_path)


batchfile = r"\\marty-nas.cm.utexas.edu\MartyLab\Group Share\Annika\DDA_Compare\mzmine\Stellar_Pos_batch.mzbatch"
outfile = batchfile.replace(".mzbatch", "_no_na.mzbatch")
#
# create_duplicate_copy(batchfile)



tree = read_mzmine_batch_file(batchfile)
oldfile = get_spectral_library(tree)
newfile = oldfile.replace(".msp", "_no_na.msp")
tree = set_spectral_library(tree, newfile)
write_mzmine_batch_file(tree, outfile)