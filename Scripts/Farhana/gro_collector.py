import os
import re
import shutil

lipid = "POPG"

top_dir = "Z:\\Group Share\\Group\\VMD\\AqpZMutants\\WT_AA_PE_PG_CL_Dropbox"
grodir = "Z:\\Group Share\\Group\\VMD\\AqpZMutants\\WT_AA_PE_PG_CL_Dropbox\\GroAll_" + lipid + "\\cluster"

index=0
# Walk through directory and find folders that match Bound_Pose_ + lipid
for dirpath, dirnames, filenames in os.walk(top_dir):
    for dirname in dirnames:
        if re.match("Bound_Poses_" + lipid, dirname):
            print(dirname)
            # Walk through the directory and find all .gro files
            for dirpath2, dirnames2, filenames2 in os.walk(os.path.join(dirpath, dirname)):
                for filename in filenames2:
                    if filename.endswith(".gro") and "cluster" in filename:
                        print(filename)
                        # Copy filename to grodir
                        shutil.copy(os.path.join(dirpath2, filename), grodir)
                        # Rename file to include the index
                        os.rename(os.path.join(grodir, filename), os.path.join(grodir, filename.split(".")[0] + "_" + str(index) + ".gro"))
                        index += 1
