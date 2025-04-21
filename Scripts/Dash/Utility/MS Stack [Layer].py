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
rcParams['font.size'] = 12
rcParams['font.weight'] = 700
rcParams['font.sans-serif'] = "Arial"
rcParams['axes.titleweight'] = 700
rcParams['axes.titlesize'] = 16

file1 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 3to1\FDAP 3C NQ.txt"
file2 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 3to1\FDAP 3D NQ.txt"
file3 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 9to1\FDAP 9C NQ.txt"
file4 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 9to1\FDAP 9D NQ.txt"
file5 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 18to1\FDAP 18C NQ.txt"
file6 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 18to1\FDAP 18D NQ.txt"



iso1 = np.loadtxt(file1)
iso2 = np.loadtxt(file2)
iso3 = np.loadtxt(file3)
iso4 = np.loadtxt(file4)
iso5 = np.loadtxt(file5)
iso6 = np.loadtxt(file6)

iso_1 = ud.datachop(iso1, 810, 840)
iso_2 = ud.datachop(iso2, 810, 840)
iso_3 = ud.datachop(iso3, 810, 840)
iso_4 = ud.datachop(iso4, 810, 840)
iso_5 = ud.datachop(iso5, 810, 840)
iso_6 = ud.datachop(iso6, 810, 840)

fig, (ax1, ax2, ax3) = plt.subplots(3)
ax1.plot(iso_1[:,0], iso_1[:,1]/np.amax(iso_1[:,1])*100, '#FF7043', alpha=1, label='rep1')
ax1.plot(iso_2[:,0], iso_2[:,1]/np.amax(iso_2[:,1])*100, '#26C6DA', alpha=1, label='rep2' )
ax1.legend()
#ax1.set_title("3to1")
ax1.set_box_aspect(0.5)
#ax1.grid(axis= 'x')


ax2.plot(iso_3[:,0], iso_3[:,1]/np.amax(iso_3[:,1])*100, '#C25432', alpha=1, label='rep1')
ax2.plot(iso_4[:,0], iso_4[:,1]/np.amax(iso_4[:,1])*100, '#006C7A', alpha=1, label='rep2')
ax2.legend()
#ax2.set_title("9to1")
ax2.set_box_aspect(0.5)
#ax2.grid(axis= 'x')

ax3.plot(iso_5[:,0], iso_5[:,1]/np.amax(iso_5[:,1])*100, '#853A22', alpha=1, label='rep1')
ax3.plot(iso_6[:,0], iso_6[:,1]/np.amax(iso_6[:,1])*100, '#004851', alpha=1, label='rep2')
ax3.legend()
#ax3.set_title("18to1")
ax3.set_box_aspect(0.5)
#ax3.grid(axis= 'x')

fig.supxlabel('m/z')
#plt.legend()
plt.xticks()
plt.show()
