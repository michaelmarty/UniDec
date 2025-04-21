import os
from HoleFunctions import run_hole, stitch_holes
import pandas as pd

dir = "/mnt/z/Group Share/Group/VMD/AqpZMutants/WT_AA_PE_PG_CL_2/namd"
os.chdir(dir)

file1 = "step5_input.pdb"
# file2 = "step7.1_production.dcd"
file2 = "wca_s5.dcd"

run_hole(file1, file2)

stitch_holes()

