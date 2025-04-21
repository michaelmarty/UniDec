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

<<<<<<< .mine
#file1 = "Z:\Group Share\Dash\DaptoMSMS\811 isolation.txt"
#file2 = "Z:\Group Share\Dash\DaptoMSMS\819 isolation.txt"
#file3 = "Z:\Group Share\Dash\DaptoMSMS\827 Isolation.txt"
file4 = "Z:\Group Share\Dash\Dapto 35.txt"
||||||| .r2685
file1 = "Z:\Group Share\Dash\DaptoMSMS\811 isolation.txt"
file2 = "Z:\Group Share\Dash\DaptoMSMS\819 isolation.txt"
file3 = "Z:\Group Share\Dash\DaptoMSMS\827 Isolation.txt"
file4 = "Z:\Group Share\Dash\DaptoMSMS\MS1.txt"
=======
file1 = "Z:\\Group Share\\Seun\\VP4 raw data\\mass list\\230217 VP4 C8E4 0.2 M AA.txt"
file2 = "Z:\\Group Share\\Seun\\VP4 raw data\\mass list\\230217 VP4 LDAO 0.2 M AA.txt"
file3 = "Z:\\Group Share\\Seun\\VP4 raw data\\mass list\\230221 VP4 C8E4 pH5.txt"
file4 = "Z:\Group Share\Dash\E Protein Raw.txt"
>>>>>>> .r2895

#iso1 = np.loadtxt(file1)
#iso2 = np.loadtxt(file2)
#iso3 = np.loadtxt(file3)
#iso4 = np.loadtxt(file3)

<<<<<<< .mine
#iso_1 = ud.datachop(iso1, 100, 2000)
#iso_2 = ud.datachop(iso2, 100, 2000)
#iso_3 = ud.datachop(iso3, 100, 2000)
iso_4 = ud.datachop(iso4, 100, 2000)
||||||| .r2685
iso_1 = ud.datachop(iso1, 1275, 1325)
iso_2 = ud.datachop(iso2, 1275, 1325)
iso_3 = ud.datachop(iso3, 1275, 1325)
iso_4 = ud.datachop(iso4, 800, 900)
=======
iso_1 = ud.datachop(iso1, 1000, 2100)
iso_2 = ud.datachop(iso2, 1000, 2100)
iso_3 = ud.datachop(iso3, 1000, 2100)
#iso_4 = ud.datachop(iso3, 1000, 3000)
>>>>>>> .r2895

<<<<<<< .mine
fig,(ax4) = plt.subplots()
#ax1.plot(iso_1[:,0], iso_1[:,1]/np.amax(iso_1[:,1])*100, 'green')
#ax2.plot(iso_2[:,0], iso_2[:,1]/np.amax(iso_2[:,1])*100, 'blue')
#ax3.plot(iso_3[:,0], iso_3[:,1]/np.amax(iso_3[:,1])*100, 'darkorchid')
ax4.plot(iso_4[:,1], iso_4[:,0]/np.amax(iso_4[:,0])*100, 'red')
||||||| .r2685
fig, (ax1, ax2, ax3, ax4) = plt.subplots(4)
ax1.plot(iso_1[:,0], iso_1[:,1]/np.amax(iso_1[:,1])*100, 'green')
ax2.plot(iso_2[:,0], iso_2[:,1]/np.amax(iso_2[:,1])*100, 'blue')
ax3.plot(iso_3[:,0], iso_3[:,1]/np.amax(iso_3[:,1])*100, 'darkorchid')
ax4.plot(iso_4[:,0], iso_4[:,1]/np.amax(iso_4[:,1])*100, 'red')
=======
fig, (ax1, ax2, ax3) = plt.subplots(3)
ax1.plot(iso_1[:,0], iso_1[:,1]/np.amax(iso_1[:,1])*100, 'green')
ax1.set_box_aspect(0.5)
ax2.plot(iso_2[:,0], iso_2[:,1]/np.amax(iso_2[:,1])*100, 'blue')
ax2.set_box_aspect(0.5)
ax3.plot(iso_3[:,0], iso_3[:,1]/np.amax(iso_3[:,1])*100, 'darkorchid')
ax3.set_box_aspect(0.5)
#ax4.plot(iso_4[:,0], iso_4[:,1]/np.amax(iso_4[:,1])*100, 'red')
#ax4.set_box_aspect(0.5)
>>>>>>> .r2895

<<<<<<< .mine
#plt.savefig("Z:\Group Share\Dash\DaptoMSMS\HighMass Output.pdf")
||||||| .r2685
plt.savefig("Z:\Group Share\Dash\DaptoMSMS\HighMass Output.pdf")
=======
#plt.savefig("Z:\Group Share\Dash\DaptoMSMS\HighMass Output.pdf")
plt.rcParams['figure.figsize'] = [4, 4]
>>>>>>> .r2895
plt.show()

