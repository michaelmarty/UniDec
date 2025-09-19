from matplotlib import rcParams

calpath = "Z:\Group Share\Marius Kostelic\CD-MS\Replicates\Calibration.csv"
calpath = "Z:\Group Share\Marius Kostelic\Baker Lab AAVs\\02092022 recalibration EFT n & off\\02092022_calibration_STORI_240K_.csv"
calpath = "Z:\Group Share\Marius Kostelic\Baker Lab AAVs\\02092022 recalibration EFT n & off\\02092022_calibration_STORI_240K_withoutBSA_.csv"

caldat = np.genfromtxt(calpath, delimiter=",", dtype=np.str)[1:]
print(caldat)

eng = UniDecCD()
extracts=[]
snextracts = []
for i, l in enumerate(caldat):
    f = l[0]
    eng.open_file(f)
    m = float(l[1])
    minz = float(l[2])
    maxz = float(l[3])
    ext, snext = eng.extract_intensities(m, minz, maxz, window=5, noise_mult=0, sdmult=1)

    if i == 0:
        extracts = ext
    else:
        extracts = np.concatenate((extracts, ext), axis=0)

    if i == 0:
        snextracts = snext
    else:
        snextracts = np.concatenate((snextracts, snext), axis=0)

    print("Length", len(ext))

fits, fitdat = eng.get_fit(extracts)
print(fits)
snfits, snfitdat = eng.get_fit(snextracts)
print(snfits)

fits2, fitdat2 = eng.get_fit_linear(extracts)
print(fits2)

fits3, fitdat3 = eng.get_fit_quadratic(extracts)
print(fits3)

print(slopefunc(170, 112597))
print(slopefunc(170, *fits))
print(linearmodel(170, *fits2))
print(quadraticmodel(170, *fits3))

rcParams['font.size']=16
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype']=42

bins = [int(np.amax(fitdat[:,0])-np.amin(fitdat[:,0])), 50]
plt.subplot(121)
plt.hist2d(extracts[:, 0], extracts[:, 1], bins=bins, cmap="viridis")#, cmax=100)
plt.plot(fitdat[:,0], fitdat[:,1], color="red")
plt.plot(fitdat2[:,0], fitdat2[:,1], color="orange")
plt.plot(fitdat3[:,0], fitdat3[:,1], color="y")
plt.plot(fitdat[:,0], fitdat[:,0]*112597 , color="green")
plt.subplot(122)
plt.hist2d(snextracts[:, 0], snextracts[:, 1], bins=bins, cmap="viridis")#, cmax=100)
plt.plot(snfitdat[:,0], snfitdat[:,1], color="red")
plt.show()
