import numpy
import os

#set the directory
dir = "Z:\\Group Share\\Marius Kostelic\\EVEN More mixed lipids\\Mitochondrial\\DPPC SOPE TSCL"
os.chdir(dir)
#load csv file with numpy
masses = numpy.loadtxt(fname="DPPC SOPE TSCL all.csv", delimiter = ",")
pos = masses[0:3]
#print(POPC_SOPE_pos)
neg = masses[3:6]
#print(POPC_SOPE_neg)
imid = masses[6:9]
#print(POPC_SOPE_imid)
#methyl = masses[9:12]
#print(POPC_SOPE_methyl)
#x4 is an input for what to print
#calculate average and standard deviations for input x4
x4 = imid
print (x4.mean(axis=0))
print (x4.std(axis=0))

file1 = "x values for figure 2.txt"
methyl = x4.mean(axis=0)
methyl_deviation = x4.std(axis=0)
#file4 ="POPC SOPE POPS 2mimid y values.txt"
#file5 ="POPC SOPE POPS 2mimid error.txt"
#file6 ="POPC SOPE POPS SSM 2mimid y values.txt"
#file7 ="POPC SOPE POPS SSM 2mimid y errors.txt"
file8 = "POPC line.txt"
file9 = "POPC SOPE line.txt"
file10 = "POPC SOPE POPS line.txt"
file11 = "POPC SOPE POPS SSM line.txt"
