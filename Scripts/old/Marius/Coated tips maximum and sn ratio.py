import numpy
import numpy as np
import os

dir = "Z:\\Group Share\\Marius Kostelic\\Tiny tips and teflon\\08032021 BSA a b c all same size triplicate"
os.chdir(dir)


data = np.loadtxt(fname="20210803_MK_BSA_C_3_1min_Aqu_2.txt", delimiter=" ", skiprows=(7))
#load data and skip first 7 rows

x = (data[:, 0])
#pull x data for column 1
y = (data[:, 1])

noise = 176.06
#pull y fata for column 2
maximum = numpy.max(data[:, 1])
StoN = maximum/noise
print(maximum)
print(StoN)
