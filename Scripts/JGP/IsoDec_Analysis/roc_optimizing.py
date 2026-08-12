from unidec.IsoDec.compare_tools import read_manual_annotations, compare_annotated, read_msalign_to_matchedcollection
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import *
from unidec.IsoDec.datatools import *
from unidec.IsoDec.match import *

def get_estimated_most_abundant_mass(monoiso):
    most_intense_iso = (int)(0.0006*monoiso + 0.4074)
    return monoiso + (most_intense_iso * 1.0033)

def combination_to_string(combination):
    return "CSS: " + str(combination[0]) + " AC: " + str(combination[1]) + " KD: " + str(combination[2])

def combination_to_string_no_kd(combination):
    return "CSS: " + str(combination[0]) + " AC: " + str(combination[1])

def get_possible_combinations(current_css, current_ac, current_kd,
                              css_min, css_max, css_step, ac_min, ac_max, ac_step, kd_min, kd_max, kd_step,
                              best_combinations):

    #Should also allow for the addition of a KD round along with INCREASING the CSS threshold or AC threshold
    combinations = []
    if current_css - css_step >= css_min:
        combination = [current_css - css_step, current_ac, current_kd]
        if combination not in best_combinations:
            combinations.append(combination)


    # if current_css + css_step <= css_max and current_kd + kd_step <= kd_max:
    #     combination = [current_css + css_step, current_ac, current_kd + kd_step]
    #     if combination not in best_combinations:
    #         combinations.append(combination)

    if current_ac - ac_step >= ac_min:
        combination = [current_css, current_ac - ac_step, current_kd]
        if combination not in best_combinations:
            combinations.append(combination)


    # if current_ac + ac_step <= ac_max and current_kd + kd_step <= kd_max:
    #     combination = [current_css, current_ac + ac_step, current_kd + kd_step]
    #     if combination not in best_combinations:
    #         combinations.append(combination)

    if current_kd + kd_step <= kd_max:
        combination = [current_css, current_ac, current_kd + kd_step]
        if combination not in best_combinations:
            combinations.append(combination)

    return combinations

def get_possible_combinations_no_kd(current_css, current_ac,
                                    css_min, css_max, css_step, ac_min, ac_max, ac_step,
                                    best_combinations):
    # Should also allow for the addition of a KD round along with INCREASING the CSS threshold or AC threshold
    combinations = []
    if current_css - css_step >= css_min:
        combination = [current_css - css_step, current_ac]
        if combination not in best_combinations:
            combinations.append(combination)

    # if current_css + css_step <= css_max and current_kd + kd_step <= kd_max:
    #     combination = [current_css + css_step, current_ac, current_kd + kd_step]
    #     if combination not in best_combinations:
    #         combinations.append(combination)

    if current_ac - ac_step >= ac_min:
        combination = [current_css, current_ac - ac_step]
        if combination not in best_combinations:
            combinations.append(combination)


    return combinations



def get_best_combination(combinations, engine, etd_annotated, etd_spec,
                         uvpd_annotated, uvpd_spec, tof_annotated, tof_spec, combinations_scores):

    assess_etd = True
    assess_uvpd = True
    assess_tof = False

    current_f1_scores = []
    current_precisions = []
    current_recalls = []

    for combination in combinations:
        if combination_to_string(combination) in combinations_scores.keys():
            precision = combinations_scores[combination_to_string(combination)][0]
            recall = combinations_scores[combination_to_string(combination)][1]
            current_f1_scores.append(2 * (precision * recall) / (precision + recall))
            current_precisions.append(precision)
            current_recalls.append(recall)
            continue

        truepos = 0
        falseneg = 0
        falsepos = 0

        engine.config.css_thresh = combination[0]
        engine.config.minareacovered = combination[1]
        engine.config.knockdown_rounds = combination[2]

        if assess_etd:
            engine.batch_process_spectrum(etd_spec, centroided=True)

            shared, u1, u2 = compare_annotated(etd_annotated, engine.pks, ppmtol=50, maxshift=3)

            engine.pks = MatchedCollection()

            truepos += len(shared)
            falseneg += len(u1)
            falsepos += len(u2)

        if assess_uvpd:
            engine.batch_process_spectrum(uvpd_spec, centroided=True)

            shared, u1, u2 = compare_annotated(uvpd_annotated, engine.pks, ppmtol=50, maxshift=3)

            engine.pks = MatchedCollection()

            truepos += len(shared)
            falseneg += len(u1)
            falsepos += len(u2)

        if assess_tof:
            engine.batch_process_spectrum(tof_spec, centroided=True)

            shared, u1, u2 = compare_annotated(tof_annotated, engine.pks, ppmtol=50, maxshift=3)

            engine.pks = MatchedCollection()

            truepos += len(shared)
            falseneg += len(u1)
            falsepos += len(u2)

        recall = truepos / (truepos + falseneg)
        precision = truepos / (truepos + falsepos)
        current_f1_scores.append(2 * (precision * recall) / (precision + recall))
        current_precisions.append(precision)
        current_recalls.append(recall)
        combinations_scores[combination_to_string(combination)] = [precision, recall]



    # Now find the combination that gives the highest F1 score
    best_f1_score = max(current_f1_scores)
    best_f1_score_index = current_f1_scores.index(best_f1_score)

    best_precision = current_precisions[best_f1_score_index]
    best_recall = current_recalls[best_f1_score_index]
    best_combination = combinations[best_f1_score_index]

    return best_combination, best_precision, best_recall


def get_best_combination_no_kd(combinations, engine, etd_annotated, etd_spec,
                         uvpd_annotated, uvpd_spec, tof_annotated, tof_spec, combinations_scores):

    assess_etd = True
    assess_uvpd = True
    assess_tof = False

    current_f1_scores = []
    current_precisions = []
    current_recalls = []

    for combination in combinations:
        if combination_to_string_no_kd(combination) in combinations_scores.keys():
            precision = combinations_scores[combination_to_string_no_kd(combination)][0]
            recall = combinations_scores[combination_to_string_no_kd(combination)][1]
            current_f1_scores.append(2 * (precision * recall) / (precision + recall))
            current_precisions.append(precision)
            current_recalls.append(recall)
            continue

        truepos = 0
        falseneg = 0
        falsepos = 0

        engine.config.css_thresh = combination[0]
        engine.config.minareacovered = combination[1]
        engine.config.knockdown_rounds = 10
        if assess_etd:
            engine.batch_process_spectrum(etd_spec, centroided=True)

            shared, u1, u2 = compare_annotated(etd_annotated, engine.pks, ppmtol=50, maxshift=3)

            engine.pks = MatchedCollection()

            truepos += len(shared)
            falseneg += len(u1)
            falsepos += len(u2)

        if assess_uvpd:
            engine.batch_process_spectrum(uvpd_spec, centroided=True)

            shared, u1, u2 = compare_annotated(uvpd_annotated, engine.pks, ppmtol=50, maxshift=3)

            engine.pks = MatchedCollection()

            truepos += len(shared)
            falseneg += len(u1)
            falsepos += len(u2)

        if assess_tof:
            engine.batch_process_spectrum(tof_spec, centroided=True)

            shared, u1, u2 = compare_annotated(tof_annotated, engine.pks, ppmtol=50, maxshift=3)

            engine.pks = MatchedCollection()

            truepos += len(shared)
            falseneg += len(u1)
            falsepos += len(u2)

        recall = truepos / (truepos + falseneg)
        precision = truepos / (truepos + falsepos)
        current_f1_scores.append(2 * (precision * recall) / (precision + recall))
        current_precisions.append(precision)
        current_recalls.append(recall)
        combinations_scores[combination_to_string_no_kd(combination)] = [precision, recall]



    # Now find the combination that gives the highest F1 score
    best_f1_score = max(current_f1_scores)
    best_f1_score_index = current_f1_scores.index(best_f1_score)

    best_precision = current_precisions[best_f1_score_index]
    best_recall = current_recalls[best_f1_score_index]
    best_combination = combinations[best_f1_score_index]

    return best_combination, best_precision, best_recall


colors = ['b', 'r', 'c', 'm', 'y', 'k', 'olive']

# Set backend to Agg
mpl.use('WxAgg')

etd_spectrum = np.loadtxt(r"Z:\Group Share\JGP\js8b05641_si_001\etd_spectrum.txt")
uvpd_spectrum = np.loadtxt(r"Z:\Group Share\JGP\JPST001885\UVPD_Test\uvpd_spectrum.txt")
tof_spectrum = np.loadtxt("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_spectrum.txt")
#tof_spectrum = spectrum = np.loadtxt("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_spectrum.txt")

etd_annotated = read_manual_annotations("Z:\\Group Share\\JGP\\js8b05641_si_001\\ETD Manual Annotations.txt", delimiter=' ')
uvpd_annotated = read_manual_annotations("Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\manual_annotations.txt", delimiter='\t')
tof_annotated = read_manual_annotations("Z:\\Group Share\\JGP\\MSV000094542\\raw\\RAW\\apomyoglobin_tof_manual_annotations.txt", delimiter='\t')

#Now add a point for default TopFD and eTHRASH
truepos = 0
falseneg = 0
falsepos = 0


#region TopFD Initial Processing
uvpd_topfd = read_msalign_to_matchedcollection(r"Z:\\Group Share\\JGP\\JPST001885\\UVPD_Test\\MASH_Workflows\\10uMCA25+UVPDLPT1,0mJ1pFullP_3_20240813_142813_TopFD_1_1\\input_topfd_ms2.msalign", uvpd_spectrum)
etd_topfd = read_msalign_to_matchedcollection(r"Z:\\Group Share\\JGP\\js8b05641_si_001\\MASH_Workflows\\1500_scans_200K_16 fills-qb1_20240815_120332_TopFD_1_1\\1500_scans_200K_16 fills-qb1_ms2.msalign", etd_spectrum)
tof_topfd = read_msalign_to_matchedcollection(r"Z:\\Group Share\\JGP\\MSV000094542\\raw\RAW\\MASH_Workflows\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_20250102_120613_TopFD_1_1\\ECD Apomyoglobin (16+) 6560 (2.5A, 25v) 0002_peaks_ms2.msalign", tof_spectrum)

shared, u1, u2 = compare_annotated(uvpd_annotated, uvpd_topfd, ppmtol=50, maxshift=3)

truepos += len(shared)
falseneg += len(u1)
falsepos += len(u2)

shared, u1, u2 = compare_annotated(etd_annotated, etd_topfd, ppmtol=50, maxshift=3)

truepos += len(shared)
falseneg += len(u1)
falsepos += len(u2)

shared, u1, u2 = compare_annotated(tof_annotated, tof_topfd, ppmtol=50, maxshift=3)

truepos += len(shared)
falseneg += len(u1)
falsepos += len(u2)

topFD_recall = truepos / (truepos + falseneg)
print("TopFD Recall: " + str(topFD_recall))
topFD_precision = truepos / (truepos + falsepos)
print("TopFD Precision: " + str(topFD_precision))
topFD_F1_score = 2 * (topFD_precision * topFD_recall) / (topFD_precision + topFD_recall)
print("TopFD F1 Score: " + str(topFD_F1_score))
#endregion TopFD Initial Processing

truepos = 0
falseneg = 0
falsepos = 0

#region Set initial parameters
steps_to_take = 15

css_min = 0.5
css_max = 0.999
css_step = 0.05

ac_min = 0.01
ac_max = 0.70
ac_step = 0.05

kd_min = 1
kd_max = 20
kd_step = 1

css_initial = 0.75
ac_initial = 0.2
kd_initial = 4

current_step = 0

#endregion Set initial parameters

engines = []
engine_names = []
precisions = []
recalls = []
best_combinations = []
combinations_scores = []

#Create the engines of interest
engines.append(IsoDecEngine(phaseres=8, use_wrapper=False, verbose=False))
engine_names.append("IsoDec-8bit")

engines.append(IsoDecEngine(phaseres=4, use_wrapper=False, verbose=False))
engine_names.append("IsoDec-4bit")

eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=False)
eng.predmode = 2
engines.append(eng)
engine_names.append("IsoDec-THRASH")

# engines.append(IsoDecEngine(phaseres=1, use_wrapper=False, verbose=False))
# engine_names.append("IsoDec-1bit")

# eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=False)
# eng.predmode = 3
# engines.append(eng)
# engine_names.append("IsoDec-THRASH3")


precisions = [[] for i in range(len(engines))]
recalls = [[] for i in range(len(engines))]
best_combinations = [[] for i in range(len(engines))]
combinations_scores = [{} for i in range(len(engines))]

for i in range(len(engines)):
    #Process all 3 spectra using the initial parameters
    truepos = 0
    falseneg = 0
    falsepos = 0

    eng = engines[i]

    eng.config.css_thresh = css_initial
    eng.config.minareacovered = ac_initial
    eng.config.knockdown_rounds = kd_initial

    eng.batch_process_spectrum(etd_spectrum, centroided=True)
    shared, u1, u2 = compare_annotated(etd_annotated, eng.pks, ppmtol=50, maxshift=3)
    truepos += len(shared)
    falseneg += len(u1)
    falsepos += len(u2)

    eng.pks = MatchedCollection()

    eng.batch_process_spectrum(uvpd_spectrum, centroided=True)
    shared, u1, u2 = compare_annotated(uvpd_annotated, eng.pks, ppmtol=50, maxshift=3)
    truepos += len(shared)
    falseneg += len(u1)
    falsepos += len(u2)

    eng.pks = MatchedCollection()

    eng.batch_process_spectrum(tof_spectrum, centroided=True)
    shared, u1, u2 = compare_annotated(tof_annotated, eng.pks, ppmtol=50, maxshift=3)
    truepos += len(shared)
    falseneg += len(u1)
    falsepos += len(u2)

    eng.pks = MatchedCollection()

    recall = truepos / (truepos + falseneg)
    precision = truepos / (truepos + falsepos)
    precisions[i].append(precision)
    recalls[i].append(recall)
    print("Initial precision: " + str(precision))
    print("Initial recall: " + str(recall))
    print("Initial F1 Score: " + str(2 * (precision * recall) / (precision + recall)))
    best_combinations[i].append([css_initial, ac_initial, kd_initial])


timestart = time.time()
while current_step < steps_to_take:
    for i in range(len(engines)):
        combinations = get_possible_combinations(best_combinations[i][-1][0], best_combinations[i][-1][1], best_combinations[i][-1][2],
                                                    css_min, css_max, css_step, ac_min, ac_max, ac_step, kd_min, kd_max, kd_step,
                                                    best_combinations[i])

        best_combination, best_precision, best_recall = (
            get_best_combination(combinations, engines[i],
                                 etd_annotated, etd_spectrum, uvpd_annotated, uvpd_spectrum, tof_annotated, tof_spectrum, combinations_scores[i]))


        precisions[i].append(best_precision)
        recalls[i].append(best_recall)
        best_combinations[i].append(best_combination)
        print("Best combination for " + engine_names[i] + " at step " + str(current_step) + " is: " + combination_to_string(best_combination))
        print("Precision: " + str(best_precision))
        print("Recall: " + str(best_recall))
        print("F1 Score: " + str(2 * (best_precision * best_recall) / (best_precision + best_recall)))

    print("Current step: ", current_step)
    current_step += 1

print("Time:", time.time() - timestart)

roc_df = pd.DataFrame()
for i in range(len(engines)):
    roc_df[engine_names[i] + "-Recall"] = recalls[i]
    roc_df[engine_names[i] + "-Precision"] = precisions[i]
    roc_df[engine_names[i] + "-CSS"] = [x[0] for x in best_combinations[i]]
    roc_df[engine_names[i] + "-AC"] = [x[1] for x in best_combinations[i]]
    roc_df[engine_names[i] + "-KD"] = [x[2] for x in best_combinations[i]]
try:
    roc_df.to_csv("Z:\\Group Share\\JGP\\IsoDec_figs\\uvpd_etd_p_r_curve_data3.csv")
except:
    pass



plt.figure()
# Create ROC plot

for i in range(len(engines)):
    plt.plot(recalls[i], precisions[i], linestyle='-', marker='', color=colors[i], label=engine_names[i])
plt.scatter([topFD_recall], [topFD_precision], color='purple', label="TopFD")
plt.legend()
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.xlabel("Recall", fontsize=24, fontweight='bold')
plt.ylabel("Precision", fontsize=24, fontweight='bold')
plt.xlim(0, 1.05)
plt.ylim(0, 1.05)
plt.xticks(np.arange(0,1.1,0.5), fontsize=16)
plt.yticks(np.arange(0,1.1,0.5), fontsize=16)
plt.tight_layout()
plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\uvpd_etd_roc3.pdf", dpi=300, transparent=True)
plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\uvpd_etd_roc3.png", dpi=300, transparent=True)
plt.show()
plt.figure()

exit()