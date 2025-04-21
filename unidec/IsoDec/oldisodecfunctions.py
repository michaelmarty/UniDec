

"""OLD ISODEC FUNCTIONS NO LONGER RELEVANT TO USE"""


"""ENCODING.PY"""
# @njit(fastmath=True)
# def get_digit(number, n):
#     Get the nth digit of a number in a given base
#     :param number: Float number
#     :param n: Digit to get
#     :return: Digit normalized to 0 to 1
#     """
#     a = number // 10. ** n % 10.
#     return a / 10.


# @njit(fastmath=True)
# def mz_to_bit(mz, maxlen=16):
#     """
#     Convert m/z to a bit array of digits in different bases.
#
#     :param mz: m/z value float
#     :param maxlen: Max length of the bit array
#     :return: List of digits scaled to 0 to 1 in different bases
#     """
#     if mz == 0:
#         return np.zeros(maxlen)
#
#     # Separate integer part and decimal part
#     intpart = int(mz)
#     decpart = mz - intpart
#     floatpart = decpart
#     decpart = round(decpart, 4)
#
#     # Get digits for integer part and decimal part in base 10
#     intpartdigits = [get_digit(intpart, i) for i in prange(4)]
#     decimaldigits = [get_digit(decpart * 1e4, i) for i in prange(4)]
#
#     # combine digits together into a list
#     if maxlen > 16:
#         l = maxlen
#     else:
#         l = 16
#
#     digits = np.zeros(l, dtype=float)
#     digits[2:6] = intpartdigits[::-1]
#     digits[6:10] = decimaldigits[::-1]
#     # Sacrifice the 100k digit for the decimal part
#     digits[0] = encode_mz_to_single(mz)
#     digits[1] = floatpart
#
#     # Get digits for decimal part in other bases
#     dd12 = decimal_to_other(decpart, base=12)
#     digits[10:13] = dd12
#
#     if maxlen >= 16:
#         dd7 = decimal_to_other(decpart, base=7)
#         digits[13:16] = dd7
#     if maxlen >= 19:
#         dd11 = decimal_to_other(decpart, base=11)
#         digits[16:19] = dd11
#     else:
#         return digits
#     if maxlen >= 22:
#         dd13 = decimal_to_other(decpart, base=13)
#         digits[19:22] = dd13
#     if maxlen >= 25:
#         dd15 = decimal_to_other(decpart, base=15)
#         digits[22:25] = dd15
#     if maxlen >= 28:
#         dd17 = decimal_to_other(decpart, base=17)
#         digits[25:28] = dd17
#     if maxlen >= 31:
#         dd19 = decimal_to_other(decpart, base=19)
#         digits[28:31] = dd19
#     # if maxlen >= 34:
#     #    dd23 = decimal_to_other(decpart, base=23)
#     #    digits[31:34] = dd23
#
#     return digits

# @njit(fastmath=True)
# def bit_matrix(x, maxlen=16):
#     """
#     Create a matrix of bit arrays for a list of m/z values
#     :param x: list of m/z values
#     :param maxlen: Max lengths of the bit arrays
#     :return: Matrix of bit arrays (maxlen x maxlen)
#     """
#     z = np.empty((maxlen, maxlen))
#     for i in prange(maxlen):
#         z[i] = mz_to_bit(x[i], maxlen=maxlen)
#     return z

# @njit(fastmath=True)
# def diff_matrix(x):
#     """
#     Create a matrix of differences between m/z values. Faster than numpy approach of subtract.outer
#     Set all above 1 to 0 and all below 0 to 1/abs(diff)/50.
#     This is a simple way to encode the differences and use the full range of the matrix.
#     :param x: m/z values
#     :return: matrix of differences (maxlen x maxlen)
#     """
#     lx = len(x)
#     z = np.empty((lx, lx))
#     for i in prange(lx):
#         for j in prange(lx):
#             val = (x[i] - x[j])
#             if val < 0:
#                 val = (1 / abs(val)) / 50.
#             if val > 1:
#                 val = 0
#             z[i, j] = val
#     return z

#@njit(fastmath=True)
# def encode_isodist(isodist, maxlen=16):
#     """
#     Encode an isotope distribution into a 3 channel matrix.
#     :param isodist: Isotope distribution (m/z, intensity)
#     :param maxlen: Max length of the matrix
#     :return: Encoded matrix (3 x maxlen x maxlen), indexes of the original isodist after sorting
#     """
#     indexes = np.arange(len(isodist))
#     # sort isodist by intensity and remove lower values
#     sortindex = np.argsort(isodist[:, 1])[::-1]
#     isodist = isodist[sortindex]
#     indexes = indexes[sortindex]
#
#     # Fill to correct length, pad with 0s or -1s
#     isodist2 = np.zeros((maxlen, 2))
#     indexes2 = np.zeros(maxlen) - 1
#     i1 = isodist[:maxlen]
#     isodist2[:len(i1)] = i1
#     indexes2[:len(i1)] = indexes[:maxlen]
#     isodist = isodist2
#
#     # Set up distaince matrix
#     # distmatrix = np.subtract.outer(isodist[:, 0], isodist[:, 0])
#     distmatrix = diff_matrix(isodist[:, 0])
#
#     # Create the weights matrix
#     weightsmatrix = np.outer(isodist[:, 1], isodist[:, 1])
#     weightsmatrix /= np.amax(weightsmatrix)
#
#     # Create the digit matrix
#     digitmatrix = bit_matrix(isodist[:, 0], maxlen=maxlen)
#
#     # Put it all together
#     emat = np.empty((3, maxlen, maxlen), dtype=float)
#     emat[0] = digitmatrix
#     emat[1] = distmatrix
#     emat[2] = weightsmatrix
#
#     return emat, indexes2


# def decode_emat(emat):
#     """
#     Decode an encoded matrix back into an isotope distribution
#     :param emat: Encoded matrix (3 x maxlen x maxlen)
#     :return: Isotope distribution (m/z, intensity)
#     """
#     digitmatrix = emat[0]
#     distmatrix = emat[1]
#     weightsmatrix = emat[2]
#
#     # Decode the digit matrix
#     mzs = np.zeros(len(digitmatrix))
#     for i in range(len(digitmatrix)):
#         mzs[i] = digitmatrix[i, 0] * 100000
#
#     # Decode the weights matrix
#     ints = np.zeros(len(weightsmatrix))
#     for i in range(len(weightsmatrix)):
#         ints[i] = np.sqrt(weightsmatrix[i, i])
#
#     # b1 = ints == 0
#     # ints[b1] = np.amax(ints) * 0.05
#     b1 = ints > 0
#     mzs = mzs[b1]
#     ints = ints[b1]
#
#     return np.vstack([mzs, ints]).T


# def plot_emat(emat):
#     """
#     Simple plot to view the encoded matrix
#     :param emat: Encoded matrix
#     :return: None
#     """
#     plt.subplot(221)
#     plt.imshow(emat[0], cmap='gray', aspect='auto')
#     plt.subplot(222)
#     plt.imshow(emat[1], cmap='gray', aspect='auto')
#     plt.subplot(223)
#     plt.imshow(emat[2], cmap='gray', aspect='auto')
#     plt.subplot(224)
#     # Color plot for rgb channels
#     rgb = emat.transpose(1, 2, 0)
#     rgb *= 255 / np.amax(rgb)
#     rgb = rgb.astype(np.uint8)
#     plt.imshow(rgb, aspect='auto')
#     plt.show()

# def save_encoding(data, outfile):
#     emat = [d[0] for d in data]
#     centroids = np.array([d[1] for d in data], dtype=object)
#     z = [d[2] for d in data]
#     print("Saving to:", outfile)
#     np.savez_compressed(outfile, emat=emat, centroids=centroids, z=z)
#
#
# def encode_dir(pkldir, outdir=None, name="medium", maxfiles=None, plot=False, **kwargs):
#     startime = time.perf_counter()
#     training = []
#     test = []
#     zdist = []

    # files = ud.match_files_recursive(pkldir, ".pkl")
    #
    # if maxfiles is not None:
    #     files = files[:maxfiles]
    #
    # print("Files:", files)
    #
    # for file in files:
    #     if "bad_data" in file:
    #         continue
    #
    #     tr, te, zd = encode_phase_file(file, **kwargs)
    #
    #     training.extend(tr)
    #     test.extend(te)
    #     zdist.extend(zd)
    #
    # # Write out everything
    # print(len(training), len(test))
    # if outdir is not None:
    #     if not os.path.isdir(outdir):
    #         os.mkdir(outdir)
    #     os.chdir(outdir)
    #
    # save_encoding(training, "training_data_" + name + ".npz")
    # save_encoding(test, "test_data_" + name + ".npz")
    #
    # # torch.save(training, "training_data_" + name + ".pth")
    # # torch.save(test, "test_data_" + name + ".pth")
    #
    # endtime = time.perf_counter()
    # print("Time:", endtime - starttime, len(zdist))
    # if plot:
    #     plt.hist(zdist, bins=np.arange(0.5, 50.5, 1))
    #     plt.show()

#@njit(fastmath=True)
# def charge_phase_calculator(centroids, maxz=50, phaseres=16, remove_harmonics=True):
#     """
#     Calculate the charge phases for a set of centroids
#     :param centroids: Centroids (m/z, intensity)
#     :param maxz: Maximum charge state to calculate
#     :return: Charge phases (maxz x maxlen)
#     """
#     phases = np.zeros((maxz, phaseres))
#     indexes = np.zeros((maxz, len(centroids)))
#     for i in range(1, maxz + 1):
#         rescale = centroids[:, 0] * 2 * np.pi * i / mass_diff_c
#         y = np.sin(rescale)
#         z = np.cos(rescale)
#         phase = (np.arctan2(z, y) / (2 * np.pi)) % 1
#         phaseindexes = np.floor(phase * phaseres)
#         indexes[i - 1] = phaseindexes
#         for j in range(len(centroids)):
#             phases[i - 1, int(phaseindexes[j])] += centroids[j, 1]
#
#     maxes = np.array([np.max(p) for p in phases])
#
#     if remove_harmonics:
#         peakbool = maxes > np.amax(maxes) * .95
#         if len(peakbool) > 1:
#             # Find first peak that is true
#             top_z = np.argmax(peakbool) + 1
#         else:
#             top_z = np.argmax(maxes) + 1
#     else:
#         top_z = np.argmax(maxes) + 1
#
#     best_phase = np.argmax(phases[top_z - 1])
#     mask = indexes[top_z - 1] == best_phase
#     return top_z, mask




"""DATATOOLS.PY"""
# @njit(fastmath=True)
# def calculate_cosinesimilarity(dist1, dist2):
#     ss = 10
#     ab = 0
#     a2 = 0
#     b2 = 0
#     for i in range(len(dist1)):
#         ab += dist1[i] * dist2[i]
#         a2 += dist1[i] ** 2
#         b2 += dist2[i] ** 2
#     return ab / (a2 ** 0.5 * b2 ** 0.5)


# @njit(fastmath=True)
# def fastwithinppmtol(array, target, ppmtol):
#     result = np.empty(0, dtype=np.int64)
#
#     if array.size == 0:
#         return result
#
#     nearest_idx = int(fastnearest(array, target))
#
#     if nearest_idx < 0 or nearest_idx >= array.size:
#         return result
#
#
#     if ud.within_ppm(array[nearest_idx], target, ppmtol):
#         result = np.append(result, nearest_idx)
#
#     adding_upper = True
#     adding_lower = True
#
#     current_upper = nearest_idx + 1
#     current_lower = nearest_idx - 1
#
#     while adding_upper or adding_lower:
#         if adding_upper:
#             if current_upper >= array.size:
#                 adding_upper = False
#             elif ud.within_ppm(array[current_upper], target, ppmtol):
#                 result = np.append(result, current_upper)
#                 current_upper += 1
#             else:
#                 adding_upper = False
#
#         if adding_lower:
#             if current_lower < 0:
#                 adding_lower = False
#             elif ud.within_ppm(array[current_lower], target, ppmtol):
#                 result = np.append(result, current_lower)
#                 current_lower -= 1
#             else:
#                 adding_lower = False
#
#     return result


# @njit(fastmath=True)
# def bisect_left(a, x):
#     """Similar to bisect.bisect_left(), from the built-in library."""
#     M = len(a)
#     for i in range(M):
#         if a[i] >= x:
#             return i
#     return M
#
#
# @njit(fastmath=True)
# def bisect_left(a, x):
#     """Return the index where to insert item x in list a, assuming a is sorted.
#
#     The return value i is such that all e in a[:i] have e < x, and all e in
#     a[i:] have e >= x.  So if x already appears in the list, a.insert(i, x) will
#     insert just before the leftmost x already there.
#     """
#
#     lo = 0
#     hi = len(a)
#
#     while lo < hi:
#         mid = (lo + hi) // 2
#         if a[mid] < x:
#             lo = mid + 1
#         else:
#             hi = mid
#
#     return lo'''
#
# '''
# @njit(fastmath=True)
# def fastnearest2(array, target):
#     """
#     In a sorted array, quickly find the position of the element closest to the target.
#     :param array: Array
#     :param target: Value
#     :return: np.argmin(np.abs(array - target))
#     """
#     i = int(bisect_left(array, target))
#     print(i)
#     if i <= 0:
#         return 0
#     elif i >= len(array) - 1:
#         return len(array) - 1
#
#     if np.abs(array[i] - target) > np.abs(array[i + 1] - target):
#         return i+1
#     elif np.abs(array[i] - target) > np.abs(array[i - 1] - target):
#         return i-1
#     return int(i)

# def remove_noise_cdata2(data, min_snr=2):
#     ydata = data[:, 1]
#     noise = np.median(ydata)
#     # plt.hist(ydata, bins=1000)
#     # plt.show()
#     max_signal = data[np.argmax(ydata), 1]
#     print("Max SNR:", max_signal / noise)
#     min_signal = noise * min_snr
#     data = data[ydata > min_signal]
#     return data

# def get_top_peak_mz(data):
#     """
#     Get the m/z value of the top peak in the data.
#     :param data: 2D numpy array of data
#     :return: float, m/z value of the top peak
#     """
#     # get the index of the maximum value in the data
#     maxindex = np.argmax(data[:, 1])
#     # return the m/z value of the peak
#     return data[maxindex, 0]


# def calc_match_simp(centroids, isodist):
#     spectrum1 = matchms.Spectrum(mz=centroids[:, 0], intensities=centroids[:, 1], metadata={"precursor_mz": 1})
#     spectrum2 = matchms.Spectrum(mz=isodist[:, 0], intensities=isodist[:, 1], metadata={"precursor_mz": 1})
#     cosine_greedy = matchms.similarity.CosineGreedy(tolerance=0.01)
#     score = cosine_greedy.pair(spectrum1, spectrum2)
#     print(score)
#     return score




"""PLOTS.PY"""

# def plot_zdist(eng):
#     zdata = [x[1] for x in eng.training_data]
#     plt.hist(zdata, bins=range(0, 50))
#     plt.show()


"""MZML.PY"""
#Slower get data memory safe

# Old Slow Method
# for i in range(int(scan_range[0]) + 1, scan_range[1] + 1):
#     try:
#         data = get_data_from_spectrum(self.msrun[self.ids[i]])
#         newdat = ud.mergedata(template, data)
#         template[:, 1] += newdat[:, 1]
#     except Exception as e:
#         print("Error", e, "With scan number:", i)