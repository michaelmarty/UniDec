from unidec.IsoDec.IsoGen.isogenpep import *
import os
import numpy as np
from Scripts.JGP.IsoGen_Analysis.distribution_assessment import calculate_cosinesimilarity
import fnmatch

def match_files(directory, string, exclude=None):
    files = []
    for file in os.listdir(directory):
        if fnmatch.fnmatch(file, string):
            if exclude is None or exclude not in file:
                files.append(file)
    return np.array(files)


if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    topdir = r"C:\Users\Admin\Documents\martylab\Protein\IsotopeScrambling\Turnover"
    os.chdir(topdir)

    train_folder = os.path.join(topdir, "Training")

    trainfiles = match_files(r"C:\Users\Admin\Documents\martylab\Protein\IsotopeScrambling\Turnover\Training", "*.npz")
    trainfiles = [os.path.join(train_folder, str(f)) for f in trainfiles]
    assessment_files = match_files(r"C:\Users\Admin\Documents\martylab\Protein\IsotopeScrambling\Turnover\Assessment", "*.npz")

    assessment_dists = []
    assessment_seqs = []

    for file in assessment_files:
        data = np.load("Assessment\\" + file)
        assessment_seqs.extend(data['seqs'])
        assessment_dists.extend((data['dists']))

    eng = IsoGenPepEngine(isolen=16)


    untrained_css_vals = []
    for i in range(len(assessment_seqs)):
        dist = eng.predict(assessment_seqs[i])
        css = calculate_cosinesimilarity(dist, assessment_dists[i])
        untrained_css_vals.append(css)

        # plt.plot(dist, label="NN")
        # plt.plot(assessment_dists[i], label="Exp")
        # plt.legend()
        # plt.xlim(0, 15)
        # plt.ylim(bottom=0)
        # plt.show()

    avg_untrained = np.average(np.array(untrained_css_vals))
    print("Untrained Average CSS: " + str(avg_untrained))

    eng.train_multiple(trainfiles, epochs=10, forcenew=False)


    transfer_trained_css_vals = []

    for i in range(len(assessment_seqs)):
        dist = eng.predict(assessment_seqs[i])
        css = calculate_cosinesimilarity(dist, assessment_dists[i])
        transfer_trained_css_vals.append(css)

        # plt.plot(dist, label="NN")
        # plt.plot(dists[i], label="Exp")
        # plt.legend()
        # plt.xlim(0,15)
        # plt.ylim(bottom=0)
        # plt.show()

    avg_trained = np.average(np.array(transfer_trained_css_vals))
    print("Transfer-Trained Average CSS: ", str(avg_trained))

    plt.hist(untrained_css_vals, alpha=0.7)
    plt.hist(transfer_trained_css_vals, alpha=0.7)
    plt.show()

