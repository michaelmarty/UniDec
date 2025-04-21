import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import unidec.tools as ud

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 18
rcParams['font.weight'] = 700
rcParams['font.sans-serif'] = "Arial"



File1 = "Z:\\Group Share\\Dash\\Text Files\\Ratios\\FDAP 18to1\\FDAP 18C NQ.txt"
File2 = "Z:\\Group Share\\Dash\\Text Files\\Ratios\\FDAP 18to1\\FDAP 18D NQ.txt"
File3 = "Z:\\Group Share\\Dash\\Text Files\\Ratios\\FPOP FDAP 18TO1\\FPOP FDAP 18A.txt"
File4 = "Z:\\Group Share\\Dash\\Text Files\\Ratios\\FPOP FDAP 18TO1\\FPOP FDAP 18B.txt"

#convert to numpy array
iso_1 = np.loadtxt(File1)
iso_2 = np.loadtxt(File2)
iso_3 = np.loadtxt(File3)
iso_4 = np.loadtxt(File4)

#chopdata
iso1A = ud.datachop(iso_1, 810, 835)
iso2A = ud.datachop(iso_2, 810, 835)
iso3A = ud.datachop(iso_3, 810, 835)
iso4A = ud.datachop(iso_4, 810, 835)

iso1B = ud.datachop(iso_1, 540, 540)
iso2B = ud.datachop(iso_2, 540, 570)
iso3B = ud.datachop(iso_3, 540, 570)
iso4B = ud.datachop(iso_4, 540, 570)


fig, (ax2) = plt.subplots(1)
#ax1.plot(iso1B[:,0], iso1B[:,1]/np.amax(iso1B[:,1])*100, '#491d8b', label= 'PLIMB 40s Rep 1')
#ax1.plot(iso2B[:,0], iso2B[:,1]/np.amax(iso2B[:,1])*100, '#004144', label= 'PLIMB 40s Rep 2')
#ax1.plot(iso3B[:,0], iso3B[:,1]/np.amax(iso3B[:,1])*((-1)*100), 'magenta',  label= 'FPOP Rep 1')
#ax1.plot(iso4B[:,0], iso4B[:,1]/np.amax(iso4B[:,1])*((-1)*100), 'blue', label= 'FPOP Rep 2')
#ax1.legend()
#ax1.set_box_aspect(0.75)

ax2.plot(iso1A[:,0], iso1A[:,1]/np.amax(iso1A[:,1])*100, '#491d8b', label= 'PLIMB Rep 1')
ax2.plot(iso2A[:,0], iso2A[:,1]/np.amax(iso2A[:,1])*100, '#004144', label= 'PLIMB Rep 2')
ax2.plot(iso3A[:,0], iso3A[:,1]/np.amax(iso3A[:,1])*((-1)*100), 'magenta', label= 'FPOP Rep 1')
ax2.plot(iso4A[:,0], iso4A[:,1]/np.amax(iso4A[:,1])*((-1)*100), 'blue', alpha=0.9, label= 'FPOP Rep 2')
ax2.legend()
#ax2.set_box_aspect(1)

plt.show()
plt.savefig("Z:\Group Share\Dash\Results.pdf")
