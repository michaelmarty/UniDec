#
#
#
# path = "C:\\Users\\Michael Marty\\Box Sync\\Grants\\2016_R21\\Data"
#
# files = ["peptide 1to118C.asc"]
# row = [1, 1, 1, 1]
# corr = [0, 0, 0, 1.52]
# fig=plt.figure(figsize=(4, 3))
#
# for i, f in enumerate(files):
#     fname = os.path.join(path, f)
#     data = np.loadtxt(fname, skiprows=3)
#     ydata = data[:, row[i]]
#     ydata = ydata - np.amin(ydata)
#     ydata = ydata / np.amax(ydata)
#     xdata = data[:, 0] - corr[i]
#
#     b1 = xdata > 5
#     b2 = xdata < 18
#     b3 = np.logical_and(b1, b2)
#
#     plt.plot(xdata[b3], ydata[b3])
#     pex = peakdetect(data, window=3, threshold=0.1)
#
#     print(pex)
#
# calfile = "calibration.asc"
# cdata = np.loadtxt(os.path.join(path, calfile), skiprows=3)
# cdata[:, 1] = cdata[:, 1] / np.amax(cdata[:, 1])
# # plt.plot(cdata[:, 0], cdata[:, 1])
#
# peaks = peakdetect(cdata, window=3, threshold=0.5)
# diams = [17, 12.2, 10.4, 7.1]
#
# pdata = np.transpose([peaks[:, 0], diams])
# for p in pdata:
#     plt.vlines(p[0], 0.9, 1)
#     plt.text(p[0]-0.5, 1.15, str(p[1]),rotation=45,ha="left")
# for p in pex:
#     try:
#         print(linterpolate(pdata, p[0]))
#     except:
#         pass
#
# b1 = xdata > 10
# b2 = xdata < 12
# b4 = np.logical_and(b1, b2)
# popt, junk, fitdat= fit_peak(xdata[b4],ydata[b4],0,10.7,0.44,1,0)
#
# fitdat = psfit(xdata[b3], popt[0], popt[1], popt[2], 0, 0)
#
# plt.plot(xdata[b3],fitdat)
# plt.text(10.5,0.4,"12.5 nm")
# plt.vlines(10.73, 0.25,0.35)
#
# fitdata=np.transpose([xdata[b3],fitdat])
# outdata=np.transpose([xdata[b3], ydata[b3]])
# print(integrate(fitdata,5,18)[0]/integrate(outdata,5,18)[0])
#
#
# plt.xlabel("Retention Time (min)\nStokes Diameter (nm)")
# plt.ylabel("A280")
#
# fig.tight_layout()
# plt.subplots_adjust(left=0.2, right=0.9, top=0.8, bottom=0.2)
# plt.savefig(os.path.join(path,"figure_1.pdf"))
# plt.show()