import os
from HoleFunctions import *


dir = "/mnt/z/Group Share/Group/VMD/AqpZMutants/WT_AA_PEPG_1/namd"
os.chdir(dir)

file1 = "step5_input.pdb"
# file2 = "step7.1_production.dcd"
file2 = "wca_s5.dcd"

run_hole(file1, file2)
stitch_holes()

