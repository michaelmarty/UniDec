import matplotlib
from LipiDecEng import *
import time

matplotlib.use('WxAgg')

ignore_list = ["Cyclo-PE", "Ala-PE", "Gly-PE", "Succ-PE", "Ala-PG", "Gly-PG", "PI-O", "PI-P"]
#ignore_list = []

# Read the lib file
nonderlibfile = "C:\\Data\\Lipidomics\\Libraries\\nonder1.npz"
derlibfile = "C:\\Data\\Lipidomics\\Libraries\\der1.npz"

# Read the data file
# dfile = "C:\Data\Lipidomics\Infusion\Raw data\"
# dfile = "C:\Data\Lipidomics\Infusion\Raw data\B01-CLMut__BJ-5ta--B1_10in4_MS_MN_211110_NonDer_Neg_211110.raw"

st = time.perf_counter()

topdir = "C:\Data\Lipidomics\Infusion\Raw data"
files = ud.match_files(topdir, "*.raw")
#files = ["B2SW620_3_10in40_RF10_NonDer_Pos.raw", "B2SW620_3_10in40_RF10_NonDer_Neg.raw"]
files = ["C01AB1T1-Human__KANJ--B1_4X_MS_MN_211203_Der_Pos_1.raw"]
files = ["B3_-5ta_CLMut_2in1_FAIMSCV62(5).raw"]
plot=True
for dfile in files:
    drange = [700,750]
    #drange = None
    if "NonDer" in dfile:
        libfile = nonderlibfile
    elif "Der" in dfile:
        libfile = derlibfile
    else:
        libfile = nonderlibfile
    runner = LipiDecRunner(dfile, libfile, dir=topdir, datarange=drange)
    runner.eng.tolerance = 5
    #runner.eng.minisomatches = 3
    runner.run(ignore_list=ignore_list)
    if plot:
        runner.make_plot(show=False)
if plot:
    plt.show()

print("All Complete...Time:", time.perf_counter()-st)

