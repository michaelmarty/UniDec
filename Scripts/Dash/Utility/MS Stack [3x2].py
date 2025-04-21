import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import unidec.tools as ud
from matplotlib import gridspec

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

file1 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 18to1\FDAP 18A.txt"
file2 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 18to1\FDAP 18B.txt"
file3 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 9to1\FDAP 9A.txt"
file4 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 9to1\FDAP 9B.txt"
file5 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 18to1\FDAP 18A.txt"
file6 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 18to1\FDAP 18B.txt"

file7 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 3to1\FDAP 3C.txt"
file8 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 3to1\FDAP 3D.txt"
file9 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 9to1\FDAP 9C.txt"
file10 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 9to1\FDAP 9D.txt"
file11 = "Z:\Group Share\Dash\Text Files\Ratios\FDAP 18to1\FDAP 18C.txt"
file12 = "Z:\Group Share\Dash\E Protein Raw.txt"


iso1 = np.loadtxt(file12)
iso2 = np.loadtxt(file12)
iso3 = np.loadtxt(file12)
iso4 = np.loadtxt(file12)
iso5 = np.loadtxt(file12)
iso6 = np.loadtxt(file12)

iso7 = np.loadtxt(file12)
iso8 = np.loadtxt(file12)
iso9 = np.loadtxt(file12)
iso10 = np.loadtxt(file12)
iso11 = np.loadtxt(file12)
iso12 = np.loadtxt(file12)


iso_1 = ud.datachop(iso1, 2000, 10000)
iso_2 = ud.datachop(iso2, 2000, 10000)
iso_3 = ud.datachop(iso3, 2000, 10000)
iso_4 = ud.datachop(iso4, 2000, 10000)
iso_5 = ud.datachop(iso5, 2000, 10000)
iso_6 = ud.datachop(iso6, 2000, 10000)
iso_7 = ud.datachop(iso7, 2000, 10000)
iso_8 = ud.datachop(iso8, 2000, 10000)
iso_9 = ud.datachop(iso9, 2000, 10000)
iso_10 = ud.datachop(iso10, 2000, 10000)
iso_11 = ud.datachop(iso11, 2000, 10000)
iso_12 = ud.datachop(iso12, 2000, 10000)

fig, ((ax1, ax4), (ax2, ax5), (ax3, ax6)) = plt.subplots(3, 2, sharex=False, constrained_layout=True)



ax1.plot(iso_1[:,0], iso_1[:,1]/np.amax(iso_1[:,1])*100, 'magenta', alpha=0.75, label='rep1')
ax1.plot(iso_2[:,0], iso_2[:,1]/np.amax(iso_2[:,1])*100, 'blue', alpha=0.5, label='rep2' )
ax1.legend()
ax1.set_xlabel("m/z")
ax1.set_title("3to1_20s")
ax1.set_box_aspect(0.5)
#ax1.grid()


ax2.plot(iso_3[:,0], iso_3[:,1]/np.amax(iso_3[:,1])*100, 'magenta', alpha=0.8, label='rep1')
ax2.plot(iso_4[:,0], iso_4[:,1]/np.amax(iso_4[:,1])*100, 'blue', alpha=0.75, label='rep2')
ax2.legend()
ax2.set_xlabel("m/z")
ax2.set_title("9to1_20s")
ax2.set_box_aspect(0.5)
#ax2.grid()


ax3.plot(iso_5[:,0], iso_5[:,1]/np.amax(iso_5[:,1])*100, 'magenta', alpha=1, label='rep1')
ax3.plot(iso_6[:,0], iso_6[:,1]/np.amax(iso_6[:,1])*100, 'blue', alpha=0.9, label='rep2')
ax3.legend()
ax3.set_xlabel("m/z")
ax3.set_title("18to1_20s")
ax3.set_box_aspect(0.5)
#ax3.grid()



ax4.plot(iso_7[:,0], iso_7[:,1]/np.amax(iso_7[:,1])*100, '#8a3ffc', alpha=1, label='rep1')
ax4.plot(iso_8[:,0], iso_8[:,1]/np.amax(iso_8[:,1])*100, '#007d79', alpha=1, label='rep2' )
ax4.legend()
ax4.set_title("3to1_40s")
ax4.set_box_aspect(0.5)
#ax4.grid()


ax5.plot(iso_9[:,0], iso_9[:,1]/np.amax(iso_9[:,1])*100, '#6929c4', alpha=1, label='rep1')
ax5.plot(iso_10[:,0], iso_10[:,1]/np.amax(iso_10[:,1])*100, '#005d5d', alpha=1, label='rep2')
ax5.legend()
ax5.set_title("9to1_40s")
ax5.set_box_aspect(0.5)
#ax5.grid()


ax6.plot(iso_11[:,0], iso_11[:,1]/np.amax(iso_11[:,1])*100, '#491d8b', alpha=1, label='rep1')
ax6.plot(iso_12[:,0], iso_12[:,1]/np.amax(iso_12[:,1])*100, '#004144', alpha=1, label='rep2')
ax6.legend()
ax6.set_title("18to1_40s")
ax6.set_box_aspect(0.5)
#ax6.grid()

plt.show()
