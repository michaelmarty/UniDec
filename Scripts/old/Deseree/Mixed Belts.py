import sys

sys.path.append("C:\\Python\\UniDec")

dir = "Z:\Group Share\Bremen\POPC_D1T0-1m_HCD\\20170207_P1C_POPC_ND_D1T0-1m_HCD_RAMP_0_240_20_1_100.0_unidecfiles"
file = "20170207_P1C_POPC_ND_D1T0-1m_HCD_RAMP_0_240_20_1_100.0_rawdata.txt"

path = os.path.join(dir, file)

dir = "Z:\Group Share\James\FIGURES\\20170727_JEK_ND_AmtB_mixedbelt_T0T1_HCD_Ramp_0_240_20\\20170727_JEK_ND_AmtB_mixedbelt_T0T1_HCD_Ramp_0_240_20_120_unidecfiles"
file = "20170727_JEK_ND_AmtB_mixedbelt_T0T1_HCD_Ramp_0_240_20_120_rawdata.txt"

path2 = os.path.join(dir, file)

dir = "Z:\Group Share\James\FIGURES\\20170801_JEK_ND_AmtB_T0T2_HCD_Ramp_0_240_20\\20170801_JEK_ND_AmtB_T0T2_HCD_Ramp_0_240_20_120_unidecfiles"
file = "20170801_JEK_ND_AmtB_T0T2_HCD_Ramp_0_240_20_120_rawdata.txt"

path3 = os.path.join(dir, file)

dir = "Z:\Group Share\James\FIGURES\\20170804_JEK_ND_AmtB_T1T2_HCD_Ramp_0_240_20\\20170804_JEK_ND_AmtB_T1T2_HCD_Ramp_0_240_20_120_unidecfiles"
file = "20170804_JEK_ND_AmtB_T1T2_HCD_Ramp_0_240_20_120_rawdata.txt"

path4 = os.path.join(dir, file)

paths = [path, path2, path3, path4, ]
colors = ["purple", "orange", "blue", "red"]

for i, path in enumerate(paths):
    print(path)

    data = np.loadtxt(path)
    data = datachop(data, 11200, 11600)
    xvalues = data[:,0]
    yvalues = data[:,1]
    yvalues/=np.amax(yvalues)
    plt.plot(data[:, 0] / 1000.0, data[:, 1] / np.amax(data[:, 1]) - i * 0.70, color=colors[i])
    #plt.plot (xvalues, yvalues,color=colors[i])
plt.show()

# plt.savefig('Z:\\Group Share\\Deseree\\BeltsShift.pdf')
