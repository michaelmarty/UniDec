import sys

sys.path.append("C:\\Python\\UniDec")

dir = "Z:\Group Share\Bremen\POPC_D1T0-1m_HCD\\20170207_P1C_POPC_ND_D1T0-1m_HCD_RAMP_0_240_20_1_120.0_unidecfiles"
file = "20170207_P1C_POPC_ND_D1T0-1m_HCD_RAMP_0_240_20_1_120.0_rawdata.txt"

path = os.path.join(dir, file)

dir = "Z:\Group Share\Bremen\POPC_D1T0-2m_HCD\\20170207_P1D_POPC_ND_D1T0-2m_HCD_RAMP_0_240_20_1_120.0_unidecfiles"
file = "20170207_P1D_POPC_ND_D1T0-2m_HCD_RAMP_0_240_20_1_120.0_rawdata.txt"

path2 = os.path.join(dir, file)

dir = "Z:\Group Share\Bremen\POPC_D1T1-2m_HCD\\20170207_PB_POPC_ND_D1T1-2m_HCD_RAMP_0_240_20_1_120.0_unidecfiles"
file = "20170207_PB_POPC_ND_D1T1-2m_HCD_RAMP_0_240_20_1_120.0_rawdata.txt"

path3 = os.path.join(dir, file)

dir = "Z:\Group Share\Deseree\Last Two Mixed Belts\\20170629_DJR_POPC_T0_R1_Ramp_0_200_20\\20170629_DJR_POPC_T0_R1_Ramp_0_200_20_140_unidecfiles"
file = "20170629_DJR_POPC_T0_R1_Ramp_0_200_20_140_rawdata.txt"

path4 = os.path.join(dir, file)

paths = [path, path2, path3, path4, ]
colors = ["purple", "red", "blue", "green"]

for i, path in enumerate(paths):
    print(path)

    data = np.loadtxt(path)
    data = datachop(data, 12010, 12080)
    xvalues = data[:,0]
    yvalues = data[:,1]
    yvalues/=np.amax(yvalues)
    plt.plot(data[:, 0] / 1000.0, data[:, 1] / np.amax(data[:, 1]) - i * 1.00, color=colors[i])
    #plt.plot (xvalues, yvalues,color=colors[i])
plt.show()

# plt.savefig('Z:\\Group Share\\Deseree\\BeltsShift.pdf')
