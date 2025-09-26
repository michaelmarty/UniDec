import numpy as np
import molmass as mm

omass = mm.Formula("O").monoisotopic_mass


def calc_mass_loss(c):
    h = 2 * c
    f = mm.Formula("C" + str(c) + "H" + str(h))
    mass = f.monoisotopic_mass
    return mass + omass, mass + 2 * omass


n = range(4, 13)
parent = 775.46-80

for nc in n:
    masses = calc_mass_loss(nc)
    print("Mass loss for " + str(nc) + " carbon atoms: " + str(masses))
    o3_masses = np.array([parent - m for m in masses]) + 3*omass
    o3_masses = np.round(o3_masses, 2)
    print("O3:", o3_masses)
    o6_masses = o3_masses + 3*omass
    o6_masses = np.round(o6_masses, 2)
    print("O6:", o6_masses)
    df = []
    for pm in o6_masses:
        for n2 in n:
            masses = calc_mass_loss(n2)
            df.append([pm-m for m in masses])
    df = np.array(df)
    df = np.round(df, 2)
    #print(df)

intact  = np.array([parent + 3*omass, parent + 6*omass, parent + 6*omass - 34, parent + 6*omass - 34*2, parent + 3*omass -34])
intact = np.round(intact, 2)
print("Intact:", intact)
# Calculate the molecular weight of the compound

