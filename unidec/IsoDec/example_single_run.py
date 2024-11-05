import os
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import *
import time

# Set backend to Agg
mpl.use('WxAgg')
# Load Main Engine
eng = IsoDecEngine(phaseres=4, use_wrapper=True)
# Set File Name and Dir
file2 = "C:\\Data\\IsoNN\\test2.txt"
os.chdir(os.path.dirname(file2))

# Load spectrum and process
spectrum = np.loadtxt(file2, skiprows=0)
print("Running")
starttime = time.perf_counter()
eng.batch_process_spectrum(spectrum, centroided=True) # Need to set this false if data isn't centroided

# Load Thermo data and process
# eng.process_file(file2)

print("Number of peaks:", len(eng.pks.peaks), "Time:", time.perf_counter() - starttime)

# Export pks
print("Writing Outputs")
outfile = os.path.basename(file2).replace(".txt", "_pks.txt")
eng.pks.export_tsv(outfile)
# eng.pks.export_prosightlite(outfile) # Your call with formats

# Make Plot
plot_pks(eng.pks, centroids=spectrum, zcolor=True)
plt.show()
