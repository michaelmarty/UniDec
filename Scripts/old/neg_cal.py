import numpy as np
import matplotlib.pyplot as plt
import os

refions = [386.714387, 646.524297, 906.334207, 1166.144117, 1425.954027, 1685.763937, 1945.573847, 2205.383757,
           2465.193667, 2725.003577, 2984.813487, 3244.623397, 3504.433307, 3764.243217, 4024.053127, 4283.863037,
           4543.672947, 4803.482857, 5063.292767, 5323.102677, 5582.912587, 5842.722497, 6102.532407, 6362.342317,
           6622.152227, 6881.962137, 7141.772047, 7401.581957, 7661.391867, 7921.201777, 8181.011687, 8440.821597,
           8700.631507, 8960.441417, 9220.251327, 9480.061237, 9739.871147, 9999.681057, 10779.110787, 11558.540517]

ces = 132.905433
iod = 126.904477


def csi(n):
    return iod + (ces + iod) * n


print(csi(41))


def localmaxpos(data, start, end):
    try:
        boo1 = data[:, 0] < end
        boo2 = data[:, 0] > start
        intdat = data[np.all([boo1, boo2], axis=0)]
        pos = np.argmax(intdat[:, 1])
        return intdat[pos, 0]
    except (ValueError, TypeError):
        return 0


def get_diff(d1, d2, value):
    ppmwindow = 2000
    window = ppmwindow / 1000000. * value
    print(window)
    p1 = localmaxpos(d1, value - window, value + window)
    p2 = localmaxpos(d2, value - window, value + window)
    print(p1, p2, value)
    return [p2 - p1, value - p2, value - p1, p1, p2]


os.chdir("C:\\Data\\cal")

d1 = np.loadtxt("CSI_neg_before.txt")
d2 = np.loadtxt("CSI_neg_after.txt")

diffs = np.array([get_diff(d1, d2, r) for r in refions])

order = 2
cal1 = np.polyfit(diffs[:, 3], refions, order)
cdat1 = np.polyval(cal1, diffs[:, 3])
odat1 = cdat1 - diffs[:, 3]
print(cal1)
print(np.std(cdat1 - refions))

"""
Best values: [  1.54101412e-08   1.00077531e+00   1.86397570e-02]
"""

cal2 = np.polyfit(diffs[:, 3], diffs[:, 4], order)
cdat2 = np.polyval(cal2, diffs[:, 3])
print(cal2)
print(np.std(cdat2 - refions))

newx = np.polyval(cal1, d1[:, 0])

if True:
    plt.plot(refions, diffs[:, 0] ,label="Before vs. After")
    plt.plot(refions, diffs[:, 1])
    plt.plot(refions, diffs[:, 2])
    #plt.plot(refions, cdat1 - refions)
    #plt.plot(refions, cdat2 - refions)
    #plt.plot(refions, odat1)
    plt.xlabel("m/z")
    plt.ylabel("Delta m/z")
    plt.show()

d1[:, 1] /= np.amax(d1[:, 1])
d2[:, 1] /= np.amax(d2[:, 1])

plt.plot(d1[:, 0], d1[:, 1])
plt.plot(newx, d1[:, 1])
plt.plot(d2[:, 0], d2[:, 1])
plt.show()
