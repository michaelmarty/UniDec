import numpy as np

cell = 1e9 #lipids
lipids_per_tube = cell * 150e6

print("Lipids per tube:", lipids_per_tube)
print("Moles of lipid per tube:", lipids_per_tube / 6.022e23 * 1e9, "nmol")

lipids_per_msp = 40

target_msp_per_tube = lipids_per_tube / lipids_per_msp

print("Moles of MSP per tube:", target_msp_per_tube / 6.022e23 * 1e9, "nmol")



vol_det = 0.02
ddm = 0.04 # M
chaps = 0.016 # M
chol = 0.015 # M
tot_det = vol_det * (ddm + chaps + chol)

print("Detergent per tube:", tot_det * 1e3, "mmol")

