from unidec.IsoDec.engine import IsoDecEngine
import os
import matplotlib.pyplot as plt
from unidec.IsoDec.encoding import small_data_dirs, data_dirs

outdir = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\IsoDec\\Figures"
os.chdir(outdir)

if __name__ == "__main__":

    eng = IsoDecEngine(phaseres=4)
    topdirectory = "C:\\Data\\IsoNN\\training"

    dirs = [os.path.join(topdirectory, d) for d in small_data_dirs]
    eng.create_merged_dataloader(dirs, "phase4", noise_percent=0.0, batchsize=32, double_percent=0.0)

    zdata = [x[1] for x in eng.training_data]

    # Set y-axis to log scale
    plt.yscale('log')
    plt.hist(zdata, bins=range(0, 50))

    #turn off top and right axis
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)

    plt.xlabel("Charge State")
    plt.ylabel("Number of Peaks in Training Data")
    plt.tight_layout()
    os.chdir(outdir)
    print(os.getcwd())
    #plt.savefig("zdist.png", dpi=300)
    #plt.savefig("zdist.pdf")

    plt.show()