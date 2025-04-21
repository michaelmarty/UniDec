import unidec.tools as ud

from unidec.UniDecImporter.ImporterFactory import ImporterFactory
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


def calc_percent(path, m1, is1, m2, is2, start1, end1, start2, end2, start3, end3, isratio=1, plot=True, f1=0, f2=0):
    importer = ImporterFactory.create_importer(path, function=f1)
    data1 = importer.get_eic(m1, tolerance=0.01)
    data1s = importer.get_eic(is1, tolerance=0.01)
    data1 = ud.datachop(data1, start1, end1)
    data1s = ud.datachop(data1s, start1, end1)

    if f2 != f1:
        importer = ImporterFactory.create_importer(path, function=f2)

    data2 = importer.get_eic(m2, tolerance=0.05)
    data2 = ud.datachop(data2, start2, end2)
    i1 = ud.integrate(data1, start1, end1)[0]
    i2 = ud.integrate(data2, start2, end2)[0]


    data2s = importer.get_eic(is2, tolerance=0.05)
    data2s = ud.datachop(data2s, start3, end3)
    i1s = ud.integrate(data1s, start1, end1)[0]
    i2s = ud.integrate(data2s, start3, end3)[0]
    i1norm = i1 / i1s
    i2norm = (i2 / i2s) * isratio

    percent = (i1norm / (i1norm + i2norm)) * 100


    print("Peak Areas: ", i1, i2, i1s, i2s)
    print("Peak Area Ratios (Norm): ", i1norm, i2norm)
    print("Percent Cholesterol: ", percent)

    if plot:

        plt.subplot(211)
        plt.plot(data1[:, 0], data1[:, 1], color="royalblue", label="Cholesterol", linewidth=3)
        plt.plot(data1s[:, 0], data1s[:, 1], color="aqua", label="D7 Cholesterol", linewidth=3)
        plt.xlabel("Retention Time (mins)", fontdict={"size": 14})
        plt.ylabel("Intensity", fontdict={"size": 14})
        plt.xlim(5,10)
        #plt.ylim(0, 11000)
        plt.legend()

        plt.subplot(212)
        plt.plot(data2[:, 0], data2[:, 1], color="darkgreen", label="POPC", linewidth=3)
        plt.plot(data2s[:, 0], data2s[:, 1], color="lime", label="D7 PC", linewidth=3)
        plt.xlabel("Retention Time (mins)", fontdict={"size": 14})
        plt.ylabel("Intensity", fontdict={"size": 14})

        plt.xlim(6.5, 8.5)
        #plt.ylim(0, 40000000)
        plt.legend()

        plt.tight_layout()
        #plt.savefig(r"Z:\Group Share\Annika\PythonFigs\new(-)controlfig.pdf")
        plt.savefig(r"Z:\Group Share\Annika\PythonFigs\BCP_talk_08-6-24_t3_(+).pdf")
        #plt.savefig(r"Z:\Group Share\Annika\PythonFigs\BrainDisc.pdf")
        plt.show()


    return percent
