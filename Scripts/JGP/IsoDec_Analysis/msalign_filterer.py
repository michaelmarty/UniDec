from unidec.IsoDec.compare_tools import read_msalign_to_matchedcollection_fast
from unidec.IsoDec.match import *
from unidec.IsoDec.msalign_export import *
import fnmatch
from matplotlib import rcParams

from unidec.UniDecImporter import ImporterFactory

# Set backend to Agg
mpl.use('WxAgg')
rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['font.family'] = 'Arial'

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)

config = IsoDecConfig()

msalign_dir = "C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247\\msalign_for_toppic\\8bin\\"
rawfile_dir = "C:\\Users\\MartyLabPC\\Desktop\\MartyLab\\IsoDec\\PXD019247\\Raw\\"
filt_fract = 0.5

for file in match_files(msalign_dir, "*ms2.msalign"):
    base_filename = file.replace(".msalign_ms2.msalign", "")
    print("Processing file:", base_filename)

    rawfile = rawfile_dir + base_filename + ".raw"
    reader = ImporterFactory.create_importer(rawfile)

    ms2_mc = read_msalign_to_matchedcollection_fast(msalign_dir + file, decon_engine="IsoDec", load_precursors=True, msorder=2)
    ms1_filename = file.replace("ms2.msalign", "ms1.msalign")
    ms1_mc = read_msalign_to_matchedcollection_fast(msalign_dir + ms1_filename, decon_engine=None, load_precursors=False, msorder=1)
    ms1_scan_dict = sort_by_scan_order(ms1_mc, 1)
    ms2_scan_dict = sort_by_scan_order(ms2_mc, 2)
    ms2_scan_intensities = {}
    for scan in ms2_scan_dict.keys():
        sum_intensity = 0
        peak_count = 0
        for p in ms2_scan_dict[scan]:
            sum_intensity += p.matchedintensity
        ms2_scan_intensities[scan] = np.log(sum_intensity)

    num_to_retain = int(len(ms2_scan_intensities)*filt_fract)
    scans_to_retain = []
    ms2_scan_intensities = dict(sorted(ms2_scan_intensities.items(), key=lambda item: item[1], reverse=True))
    for i in range(num_to_retain+1):
        scans_to_retain.append(list(ms2_scan_intensities.keys())[i])

    ms2_scan_dict_filt = {}
    for scan in scans_to_retain:
        ms2_scan_dict_filt[scan] = ms2_scan_dict[scan]
    write_ms2_msalign(ms2_scan_dict_filt, ms1_scan_dict, reader, msalign_dir + base_filename + "_filt", config,
                      act_type="ETD", max_precursors=1, existing_precursors=ms2_mc.precursors)






