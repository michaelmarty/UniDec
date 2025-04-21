path = "C:\\Data\\imstest2.mzML"
path = "C:\Data\IMS Example Data\mAbs5uL_test.mzML.gz"

im = mzMLimporter(path)

data = im.get_im_data(mzbins=1)

# np.savetxt("C:\\Data\\IMS Example Data\\test.txt", data)
#plt.hist2d(data[:,0], data[:,1], weights=data[:,2], bins=100)
#plt.show()
#exit()

x = np.unique(data[:, 0])
y = np.unique(data[:, 1])
Z = np.reshape(data[:, 2], (len(x), len(y))).transpose()
print(Z.shape)
plt.contourf(x, y, np.sqrt(Z), 100, cmap="nipy_spectral")
plt.show()

exit()

'''
path = "D:\Data\ShortCourse\strep.mzML"
msrun = pymzml.run.Reader(path)

spec1 = msrun[1]
spec2 = msrun[2]

peaks1 = spec1.peaks('centroided')
peaks2 = spec2.peaks('centroided')

average_spec = spec1 + spec2

spec3 = np.array(average_spec.peaks('reprofiled'))

plt.plot(spec3[:,0], spec3[:,1]/2, color="r")
plt.plot(spec1.mz, spec1.i, color="k")
plt.plot(spec2.mz, spec2.i, color="b")
plt.xlim(2500,5000)
plt.show()
'''
path = "D:\Data\ShortCourse\strep.RAW"
from unidec import modules as dr

im = dr.DataImporter(path)

spec1 = im.get_data(scan_range=[0, 0])
spec2 = im.get_data(scan_range=[1, 1])
spec3 = im.get_data(scan_range=[2, 2])
avg = im.get_data(scan_range=[0, 2])

start = 0
end = 5
simpavg = np.average([spec1[start:end, 0], spec2[start:end, 0]], axis=0)
out = np.transpose([spec1[start:end, 0], spec2[start:end, 0], spec3[start:end, 0], avg[start:end, 0], simpavg])
print(out)
res = spec2[0, 0] / np.diff(spec2[start:end, 0])[0]
res2 = spec1[0, 0] / np.diff(spec1[start:end, 0])[0]
res = np.average([res, res2])
print(res)
startpoint = spec2[0, 0]
print(startpoint + spec2[0, 0] / res)
