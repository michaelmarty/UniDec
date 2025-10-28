from unidec.IsoDec.compare_tools import compare_matchedcollections
from unidec.IsoDec.engine import IsoDecEngine
from unidec.IsoDec.plots import *
import pyperclip
from unidec.IsoDec.encoding import *

# Set backend to Agg
mpl.use('WxAgg')

file2 = "Z:\\Group Share\\JGP\\js8b05641_si_001\\1500_scans_200K_16 fills-qb1.raw"
os.chdir(os.path.dirname(file2))

#spectrum = np.loadtxt("test893.txt")
spectrum = np.loadtxt("test.txt")
spectrum = ud.datachop(spectrum, 891.18, 894.285)
emat = encode_phase(spectrum)

flat = emat.flatten()


modelpath = "C:\\Python\\UniDec3\\unidec\\IsoDec\\phase_model_8.pth"
model = torch.load(modelpath, weights_only=True)

w1 = model['linear_relu_stack.0.weight'].cpu().numpy()
b1 = model['linear_relu_stack.0.bias'].cpu().numpy()
w2 = model['linear_relu_stack.2.weight'].cpu().numpy()
b2 = model['linear_relu_stack.2.bias'].cpu().numpy()

l1 = np.dot(flat, np.transpose(w1)) + b1

l1 = np.maximum(l1, 0)

l2 = np.dot(l1, np.transpose(w2)) + b2


#plt.imshow(np.transpose([l2]), aspect=0.25, cmap="viridis", interpolation="none")
#plt.show()

'''
reader = ImporterFactory.create_importer(file2)
spectrum = reader.grab_centroid_data(1)

spectrum = remove_noise_cdata(spectrum, 100, factor=1.5, mode="median")

mzrange= [893, 896]
if mzrange is not None:
    spectrum = ud.datachop(spectrum, mzrange[0], mzrange[1])
np.savetxt("test893.txt", spectrum)
exit()'''

verbose = True

eng = IsoDecEngine(phaseres=8, use_wrapper=False, verbose=verbose)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2

print("Simple:", eng.phase_predictor(spectrum))

exit()
timestart = time.time()
eng.batch_process_spectrum(spectrum, centroided=True)

mc1 = eng.pks

print("Time Python:", time.time() - timestart)
#exit()
eng = IsoDecEngine(phaseres=8, use_wrapper=True, verbose=verbose)
eng.config.activescan = 1
eng.config.activescanrt = -1
eng.config.activescanorder = 2

timestart = time.time()
eng.batch_process_spectrum(spectrum, centroided=True)

mc2: List[MatchedPeak] = eng.pks

print("Time C:", time.time() - timestart)

shared, u1, u2 = compare_matchedcollections(mc1, mc2, ppmtol=50, objecttocompare="monoiso", maxshift=3)

plt.figure()
plot_pks(u1,forcecolor="r", centroids=spectrum, ccolor="k", show=False, plotmass=False, nocentroids=True)
plot_pks(u2, forcecolor="b", show=False, plotmass=False, nocentroids=True)
plot_pks(shared, forcecolor="g", show=False, plotmass=False, nocentroids=True)
plt.show()
exit()

