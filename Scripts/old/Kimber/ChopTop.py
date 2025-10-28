import numpy as np

file = "Z:\\Group Share\\Kimber\\PastedSpectra\\test2.txt"
outfile = "Z:\\Group Share\\Kimber\\PastedSpectra\\test2_chopped.txt"

data = np.loadtxt(file, delimiter="\t")

# Chop all intensities above value

max = np.max(data[:,1])

chop= 0.5e6
magfactor = max/ chop
print("Mag Factor: ", magfactor)

b1 = data[:,1] > chop
data[b1,1] = chop



np.savetxt(outfile, data, delimiter="\t", fmt="%f")