from unidec.IsoDec.runtime import  *
from unidec.IsoDec.engine import *
from UniDecImporter.ImporterFactory import *
from unidec.IsoDec.plots import *
import numpy as np
from unidec.modules.isotopetools import *
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import *

os.chdir("Z:\\Group Share\\JGP\\IsoGen\\TurnoverData")

def extract_dist(matched_peak, isolength=32):
    dist = np.zeros(isolength)
    for i,m in enumerate(matched_peak.matchedcentroids):
        dist[i] = m[1]
    dist = dist / np.max(dist)
    return dist

if __name__ == "__main__":
    if False:
        eng = IsoDecEngine(phaseres=8, verbose=False, use_wrapper=False)
        eng.config.knockdown_rounds = 1

        file = "20221010_dSILAC_iMN_Day6_Fraction16.raw"
        reader = ImporterFactory.create_importer(file)
        for i in range(1, len(reader.scans)+1):
            order = reader.get_ms_order(i)
            if order != 1:
                continue
            spectrum = reader.get_single_scan(i)
            eng.batch_process_spectrum(spectrum, centroided=True)

        monoisos = []
        dists = []

        for p in eng.pks.peaks:
            monoisos.append(p.monoiso)
            dists.append(np.array(extract_dist(p)))

        np.savez_compressed("silac_dists_day6.npz", dists=dists, masses=monoisos)

    if False:
        data = np.load("silac_dists_day1.npz")

        dists = data['dists']
        masses = data['masses']

        chis = np.zeros(len(masses))
        for i in range(len(dists)):
            averagine_dist = calc_averagine_isotope_dist(masses[i], isolen=32, mono=True)
            chi = calculate_pearson_chisquared(averagine_dist[:,1], dists[i])
            chis[i] = chi
            if False:
                plt.plot(averagine_dist[:,1], label='Averagine')
                plt.plot(dists[i], label="Experimental")
                plt.title("Mass:" + str(masses[i]) + "Chi:" + str(chi))
                plt.xlim(0,10)
                plt.legend()
                plt.show()

        print("Mean chi:", np.mean(chis))