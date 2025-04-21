import os
from HoleFunctions import run_hole
import pandas as pd

dir = "/mnt/z/Group Share/Group/VMD/AqpZMutants/[location]/namd"
os.chdir(dir)

file1 = "step5_input.pdb"
# file2 = "step7.1_production.dcd"
file2 = "[dcdfile].dcd"

cps = [
    [0.51, 15.83, 4.69],
    [-15.76, 0.714, 5.576],
    [16.0175, -0.5727, 4.883],
    [-3.83, -17.9667, 6.909]
]

run_hole(file1, file2, cps)

