import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import os

rcParams['ps.useafm']=True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42
rcParams['lines.linewidth'] = 1
rcParams['axes.linewidth'] = 1
rcParams['font.size'] = 16

path = "Z:\\Group Share\\Larry Walker\\MassSpecData\\All Peptides\\Lipid_Loss\\"
os.chdir(path)

ydata1 = np.loadtxt("Bactenecin_DMPC_Averages.txt")
yerr1 = np.loadtxt("Bactenecin_DMPC_Errors.txt")
ydata2 = np.loadtxt("Bactenecin_DTT_DMPC_Averages.txt")
yerr2 = np.loadtxt("Bactenecin_DTT_DMPC_Errors.txt")
ydata3 = np.loadtxt("LactoferricinB_DMPC_Averages.txt")
yerr3 = np.loadtxt("LactoferricinB_DMPC_Errors.txt")
ydata4 = np.loadtxt("LactoferricinB_DTT_DMPC_Averages.txt")
yerr4 = np.loadtxt("LactoferricinB_DTT_DMPC_Errors.txt")
ydata5 = np.loadtxt("Magainin2_DMPC_Averages.txt")
yerr5 = np.loadtxt("Magainin2_DMPC_Errors.txt")
ydata6 = np.loadtxt("TemporinA_DMPC_Averages.txt")
yerr6 = np.loadtxt("TemporinA_DMPC_Errors.txt")
ydata7 = np.loadtxt("TemporinL_DMPC_Averages.txt")
yerr7 = np.loadtxt("TemporinL_DMPC_Errors.txt")
ydata8 = np.loadtxt("KR12_DMPC_Averages.txt")
yerr8 = np.loadtxt("KR12_DMPC_Errors.txt")
ydata9 = np.loadtxt("Indolicidin_DMPC_Averages.txt")
yerr9 = np.loadtxt("Indolicidin_DMPC_Errors.txt")

xdata = np.loadtxt("Lipid_Loss_Xvalues.txt")
xdata2 = np.loadtxt("LactoferricinB_Xvalues.txt")
xerr = np.loadtxt("Lipid_Loss_XError.txt")


plt.plot(xdata, ydata1, '-o', color='blue', linewidth=3.0, markersize=7.0)
plt.errorbar(xdata, ydata1, yerr=yerr1, color='blue', capsize=5.0, linewidth=3.0, capthick=3.0)
plt.plot(xdata, ydata2, '-o', color='steelblue', linewidth=3.0, markersize=7.0)
plt.errorbar(xdata, ydata2, yerr=yerr2, color='steelblue', capsize=5.0, linewidth=3.0, capthick=3.0)
plt.plot(xdata2, ydata3, '-o', color='limegreen', linewidth=3.0, markersize=7.0)
plt.errorbar(xdata2, ydata3, yerr=yerr3, color='limegreen', capsize=5.0, linewidth=3.0, capthick=3.0)
plt.plot(xdata2, ydata4, '-o', color='darkgreen', linewidth=3.0, markersize=7.0)
plt.errorbar(xdata2, ydata4, yerr=yerr4, color='darkgreen', capsize=5.0, linewidth=3.0, capthick=3.0)
plt.plot(xdata, ydata5, '-o', color='indigo', linewidth=3.0, markersize=7.0)
plt.errorbar(xdata, ydata5, yerr=yerr5, color='indigo', capsize=5.0, linewidth=3.0, capthick=3.0)
plt.plot(xdata, ydata6, '-o', color='black', linewidth=3.0, markersize=7.0)
plt.errorbar(xdata, ydata6, yerr=yerr6, color='black', capsize=5.0, linewidth=3.0, capthick=3.0)
plt.plot(xdata, ydata7, '-o', color='gray', linewidth=3.0, markersize=7.0)
plt.errorbar(xdata, ydata7, yerr=yerr7, color='gray', capsize=5.0, linewidth=3.0, capthick=3.0)
plt.plot(xdata, ydata8, '-o', color='sienna', linewidth=3.0, markersize=7.0)
plt.errorbar(xdata, ydata8, yerr=yerr8, color='sienna', capsize=5.0, linewidth=3.0, capthick=3.0)
plt.plot(xdata, ydata9, '-o', color='darkorange', linewidth=3.0, markersize=7.0)
plt.errorbar(xdata, ydata9, yerr=yerr9, color='darkorange', capsize=5.0, linewidth=3.0, capthick=3.0)

plt.xticks(size=20)
plt.xlim(-0.25,19)
plt.ylim(top=180)
plt.ylim(bottom=120)
plt.xticks(np.arange(0,20, step=3))
plt.xlabel("Peptdie:Nanodisc")
plt.ylabel("# Lipids")
#plt.yticks(np.arange(0,1.01, step=0.25), size=20)
#fig.tight_layout()

plt.show()
