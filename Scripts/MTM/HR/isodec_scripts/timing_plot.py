import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os

mpl.use("WxAgg")
# create bar charts for this data

timing = {"FlashDeconv": 1295.0, "TopFD": 637., "IsoDec Python": 298., "IsoDec C": 118.0}

ydata = list(timing.values())
ydata = np.array(ydata)

plt.figure(figsize = (5, 3))
plt.bar(timing.keys(), ydata/60, color=['blue', 'orange', 'green', 'red'])
plt.ylabel("Time (min)")
#plt.legend()
# Turn off top and right axis
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.tight_layout()

os.chdir("C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures")
plt.savefig("timing.png", dpi=300, transparent=True)
plt.savefig("timing.pdf")

plt.show()
