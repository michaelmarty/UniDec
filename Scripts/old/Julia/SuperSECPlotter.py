path = "C:\\Users\\marty\\University of Arizona\\martylab - Documents\\Papers\\M2\\SEC data for AM2 compiled.csv"

data = np.genfromtxt(path, delimiter=",", skip_header=2)

samples = ["LDAO 5", "C8E4 5", "OG 9", "C8E4 9", "C8E4 9 AMT"]
colors = ["r", "b", "g", "y", "k"]

# Loop through each sample
for i, s in enumerate(samples):
    #ax = plt.subplot(2, 3, i + 1)
    # Loop through each replicate
    fwhms = []
    peaks = []
    for j in [2, 4]:
        # Extract the required columns
        start = i * 6 + j
        end = i * 6 + j + 2
        d1 = data[:, start:end]
        # Remove the nan data points
        d1 = d1[~np.isnan(d1[:, 1])]
        # Narrow time range
        d1 = datachop(d1, 10, 16)
        # Subtract the minimum
        d1[:, 1] -= np.amin(d1[:, 1])
        # Normalize to the maximum
        d1[:, 1] /= np.amax(d1[:, 1])
        # Plot
        plt.plot(d1[:, 0], d1[:, 1]-i, color=colors[i], label=s)

        # Get Peaks and FWHM
        peak = d1[np.argmax(d1[:,1]), 0]
        peaks.append(peak)
        fwhm, range = calc_FWHM(peak, d1)
        fwhms.append(fwhm)
        #ax.set_title(s)
    print(s)
    print(np.mean(peaks), np.std(peaks))
    print(np.mean(fwhms), np.std(fwhms))
#plt.vlines(12.5, -4,1)
plt.legend()
plt.grid(which="both")
plt.show()
