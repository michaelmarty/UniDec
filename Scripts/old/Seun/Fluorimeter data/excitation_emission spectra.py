import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['ps.useafm'] = True
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['font.size'] = 12
mpl.rcParams['font.sans-serif'] = "Arial"



data=np.loadtxt("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20241010 Fluoresecence Data\\mruby3 1000x diluted (blank + background corrected).txt")

plt.plot(data[:,0],data[:,1], label="mRuby 3", color="red")
plt.gca().set_xlim(500,800)
plt.ylim(0.4,)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("Wavelength (nm)")
plt.ylabel("Intensity")
plt.title('Emission spectrum')
plt.legend()


# Adjust layout
plt.tight_layout()

# Save figure
#plt.savefig("Z:\\Group Share\\Seun\\LX-MS\\LC\\Chromatogram\\20240818_mRuby 3 Em blnk_bg_crt.png", dpi=1000)

plt.show()