import os
import matplotlib.pyplot as plt
import numpy as np
import unidec.tools as ud

dir = "Z:\\Group Share\\Skippy\\Projects\\HT\\Good data sets\\Bgal bit len"

lens = [1, 3, 5]
colors = ["g", "b", "m"]
noiseRanges = [(0, 6), (12.5, 15), (12.5, 36)]

for ind, len in enumerate(lens):
    plt.figure()
    print (len, "-bit")
    top_folder = os.path.join(dir, str(len) + " bit")

    chromfiles = ud.match_files_recursive(top_folder, "chrom.txt")
    print(chromfiles)
    cdata = [np.loadtxt(file) for file in chromfiles]
    #print(len(cdata))
    #plt.plot(cdata[0][:, 0], cdata[0][:, 1])

    files = ud.match_files_recursive(top_folder, "TIC_decon.txt")
    #print(files)
    datas = [np.loadtxt(file) for file in files]

    if len == 1:
        datas = cdata
    for i in range(3):
        plt.plot(datas[i][:, 0], datas[i][:, 1], label="Chromatogram", color=colors[i])
        #plt.plot(datas[i][:, 0], datas[i][:, 1], label="Deconvoluted")

    # Noise level of 1-bit data
    for i in range(3):
        # Get noise level as three times the standard deviation of x values in the first 6 minutes
        #noise = np.std(datas[i][:, 1][datas[i][:, 0] < 6])*3

        lbound = np.argmax(np.array(datas) > noiseRanges[ind][0])
        print('lbound', lbound)
        hbound = np.argmax(np.array(datas) > noiseRanges[ind][1])
        print('hbound', hbound)

        noise = np.std(datas[i][lbound:hbound]) * 3

        # plot noise level as a horizontal line
        plt.axhline(noise, color=colors[i], linestyle='--', label='Noise level')

        # calculate SNR
        snr = np.max(cdata[i][:, 1])/noise

        print('SNR', snr)

plt.legend()
plt.show()