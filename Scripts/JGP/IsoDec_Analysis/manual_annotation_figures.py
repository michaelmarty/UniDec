import matplotlib.pyplot as plt
import matplotlib as mpl

masses = []
charges = []

uvpd_manual_annotation_file = r"Z:\Group Share\JGP\JPST001885\UVPD_Test\manual_annotations.txt"
etd_manual_annotation_file = "Z:\\Group Share\\JGP\\js8b05641_si_001\\ETD Manual Annotations.txt"
tof_manual_annotation_file = r"Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\apomyoglobin_tof_manual_annotations.txt"

with open(uvpd_manual_annotation_file, "r") as file:
    for line in file:
        try:
            split = line.split("\t")
            mz = float(split[0].strip())
            z = int(split[1].strip())
            mass = mz * z
            masses.append(mass)
            charges.append(z)
        except Exception as e:
            print(e)
            continue

with open(etd_manual_annotation_file, "r") as file:
    for line in file:
        try:
            split = line.split(" ")
            mz = float(split[0].strip())
            z = int(split[1].strip())
            mass = mz * z
            masses.append(mass)
            charges.append(z)
        except Exception as e:
            print(e)
            continue

with open(tof_manual_annotation_file, "r") as file:
    for line in file:
        try:
            split = line.split("\t")
            mz = float(split[0].strip())
            z = int(split[1].strip())
            mass = mz * z
            masses.append(mass)
            charges.append(z)
        except Exception as e:
            print(e)
            continue

# Set backend to Agg
mpl.use('WxAgg')

fig, axs = plt.subplots(1,2)

axs[0].hist(masses, bins=25)
axs[0].set_xlabel("Mass (Da)")
axs[0].spines['top'].set_visible(False)
axs[0].spines['right'].set_visible(False)

min_z = min(charges)
max_z = max(charges)
range = max_z - min_z
axs[1].hist(charges, bins=range + 1)
axs[1].set_xlim(0, max_z + 1)
axs[1].set_xlabel("Charge")
axs[1].spines['top'].set_visible(False)
axs[1].spines['right'].set_visible(False)
plt.show()