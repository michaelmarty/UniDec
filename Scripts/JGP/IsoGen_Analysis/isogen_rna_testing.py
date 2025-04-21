from unidec.IsoDec.IsoGen.isogenrna import *
from unidec.IsoDec.IsoGen.isogen_tools import *

# Set backend to Agg
mpl.use('WxAgg')

eng = IsoGenRNAEngine()

seq = "AUGCUAGC"
dist = eng.predict(seq)

truedist = rnaseq_to_dist(seq)
plt.plot(dist * 100, label="AI", color="b")
plt.plot(truedist * 100, color="k", label="True")
plt.legend()
plt.title(seq)
plt.xlabel("Isotope Number")
plt.ylabel("%")
plt.xlim(0, 10)
plt.show()