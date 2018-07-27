# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 16:36:48 2018

@author: Sean
"""
from PyQt5 import QtGui
import numpy as np
from scipy import interpolate as interpolate
import os
from unidec_modules import unidectools as ud


def plot(self):
    name = QtGui.QFileDialog.getOpenFileName()
    namestr = str(name[0])
    try:
        skip = ud.header_test(namestr)
        #print("Header test on ", namestr, "found", skip)
    except:
        skip = 0

    namebase = os.path.splitext(namestr)[0]
    rawdata = open(name[0], "r")

    if namestr.endswith('.csv'):
        x, y = np.loadtxt(rawdata,
                          unpack=True,
                          delimiter=',')
    elif namestr.endswith('.txt') or namestr.endswith('.dat'):
        x, y = np.loadtxt(rawdata,
                          unpack=True,
                          skiprows=skip)
    else:
        warningMS = QtGui.QWidget()
        msg2 = QtGui.QMessageBox.warning(warningMS, "Message", "Not a Txt or CSV file")
        msg2.show()

    return x, y, namebase


#########################################################################################
# look for duplicate m/z values and replace them with averages between themselves and their closest different m/z value
#########################################################################################
def plot_function(x, y):
    numbadx = 0

    for i in range(1, (len(x) - 1)):
        if (x[i] == x[i + 1]):
            badx = x[i]
            x[i] = (x[i - 1] + badx) / 2
            x[i + 1] = (x[i] + badx) / 2
            numbadx += 1

    #########################################################################################
    #########################################################################################

    #########################################################################################
    ###########################Interpolates the data#########################################
    #########################################################################################
    datainterpolation = interpolate.InterpolatedUnivariateSpline(x, y, k=3)
    xnew = np.linspace(np.min(x), np.max(x), np.size(x), endpoint=True)  # number of equally-spaced m/z values is equal
    # to number of m/z values in original data
    ynew = datainterpolation(xnew)  # evalues the interpolation at the equally-spaced m/z values

    negvalues = 0

    for i in range(0, len(ynew)):  # set any negative interpolated abundances to 0

        if ynew[i] < 0:
            ynew[i] = 0
            negvalues += 1

        #########################################################################################
    #########################################################################################

    #########################################################################################
    ###########
    #########################################################################################
    next2power = np.ceil(np.log2(len(ynew)))
    newpadding = pow(2, next2power) - len(ynew)
    paddedynew = np.pad(ynew, (0, np.int_(newpadding)), 'constant', constant_values=(0, 0))
    paddedxnew = np.pad(xnew, (0, np.int_(newpadding)), 'constant', constant_values=(0, 0))

    span = np.max(x) - np.min(x)  # m/z span of mass spec
    expandedspan = span + newpadding / len(ynew) * span
    temp = np.array(paddedynew)
    yflip = np.fliplr([temp])[0]
    yfull = np.append(yflip, paddedynew)  # mirror the mass spectrum about its highest m/z value
    maxfreq = np.size(yfull) / expandedspan / 2
    ftspacing = maxfreq / (np.size(yfull) - 1)

    return yfull, expandedspan, xnew, ftspacing, paddedxnew, ynew


def Fourier(self, maxfreq, yfull):
    ftx = np.linspace(0, maxfreq, np.size(yfull), endpoint=True)
    FT = np.fft.fft(yfull)
    ABFT = np.absolute(FT)
    # ABFT = np.real(FT)
    return ftx, ABFT, FT


def maxfreq(self, yfull, expandedspan):
    maxfreq = np.size(yfull) / expandedspan / 2
    return maxfreq


def inputmax(chargestaterCalc, submassCalc, ABFT, ftx):
    refmaxtab = []
    for i in range(0, len(chargestaterCalc)):
        start = chargestaterCalc[i] / submassCalc
        for i in range(0, len(ftx)):
            if ftx[i] < start:
                continue
            else:
                refmaxtab.append([ftx[i], ABFT[i]])
                break
    return refmaxtab


def findmax(expandedspan, ABFT, ftx, lowend, delta, pctpkht):  # find all the local maxima in the FFT data beyond a
    # lowest frequency, with user-tunable tolerances
    lowend = lowend  # the first bin in the FFT data that gets used for finding FFT peaks and subunit mass, etc.
    pctpkht = pctpkht / 100
    delta = delta  # the half-width of the span of FFT data around a given point used to determine whether a peak is a
    # local maximum
    ABFTRange = []
    HalfLen = int(len(ftx) / 2)
    for i in range(0, HalfLen):
        if ftx[i] < lowend:
            continue
        else:
            ABFTRange.append(ABFT[i])
    minNum = pctpkht * max(ABFTRange)
    refmaxtab = []

    for i in range(0 + delta, HalfLen - delta):
        if ftx[i] < lowend:
            continue
        else:
            large = True
            baseline = True
            if ABFT[i] < minNum:
                baseline = False
            for j in range(i - delta, i + delta):
                if ABFT[i] < ABFT[j]:
                    large = False
            if large == True and baseline == True:
                refmaxtab.append([ftx[i], ABFT[i]])
    return refmaxtab


def spacing(self, refmaxtabCalc):  # calculates an initial estimate for the spacing of the group of evenly-spaced
    # FFT spectrum peaks that probably aren't overtone peaks

    numcharCalc = 0
    for i in range(1, len(refmaxtabCalc) - 1):
        if refmaxtabCalc[i + 1][0] - refmaxtabCalc[i][0] > (1.5 * (refmaxtabCalc[i][0] - refmaxtabCalc[i - 1][0])):
            numcharCalc = i
            break
        else:
            if (i < len(refmaxtabCalc) - 2):
                continue
            else:
                numcharCalc = len(refmaxtabCalc) - 1

    return numcharCalc


def omega(refmaxtabCalc, numcharCalc):
    omegaCalc = (refmaxtabCalc[numcharCalc][0] - refmaxtabCalc[0][0]) / numcharCalc
    return omegaCalc


def charge(refmaxtabCalc, numcharCalc, omegaCalc):
    chargestatesCalc = []
    chargestatesrCalc = []

    for i in range(0, numcharCalc + 1):
        ch = refmaxtabCalc[i][0]
        ch1 = (ch / omegaCalc)
        chround = round(ch1)
        chargestatesCalc.append(ch1)
        chargestatesrCalc.append(chround)

    return chargestatesCalc, chargestatesrCalc


def subunit(refmaxtabCalc, numcharCalc, omegaCalc, chargestatesCalc, chargestatesrCalc, ftxCalc, ABFTCalc):
    # uses the initial estimate of the FFT peak spacing to find the centroids of each FFT peak, then uses the
    # spacing between them to find the subunit mass
    centroids = []

    for i in range(0, numcharCalc + 1):
        numerator = 0
        denominator = 0
        submasstot = 0
        stdevtot = 0
        rawpeakcenter = (omegaCalc * chargestatesCalc[i])

        for j in range(np.int_(np.ceil((rawpeakcenter - omegaCalc / 2) / ftxCalc[1] - 1)),
                       np.int_(np.ceil((rawpeakcenter + omegaCalc / 2) / ftxCalc[1]))):
            numerator += ABFTCalc[j] * ftxCalc[j]
            denominator += ABFTCalc[j]

        centroids.append(numerator / denominator)

    for i in range(0, numcharCalc + 1):
        print(len(centroids))
        print(len(chargestatesrCalc))
        deltasmtot = (1 / centroids[i]) * chargestatesrCalc[i]
        submasstot += deltasmtot

    submass = submasstot / (numcharCalc + 1)

    for i in range(0, numcharCalc + 1):
        deltastdevtot = (1 / centroids[i] * chargestatesrCalc[i] - submass) * (
                    1 / centroids[i] * chargestatesrCalc[i] - submass)
        stdevtot += deltastdevtot

    stdevmass = np.sqrt(stdevtot / (numcharCalc))
    return submass, stdevmass


def FTRealNoiAvg(minnoi, maxnoi, FT):
    realnoiseavg = np.average(np.real(FT[minnoi:maxnoi]))
    return realnoiseavg


def FTImagNoiAvg(minnoi, maxnoi, FT):
    imagnoiseavg = np.average(np.imag(FT[minnoi:maxnoi]))
    return imagnoiseavg


def FTRealNoiStdev(minnoi, maxnoi, FT):
    realnoisestdev = np.std(np.real(FT[minnoi:maxnoi]))
    return realnoisestdev


def FTImagNoiStdev(minnoi, maxnoi, FT):
    imagnoisestdev = np.std(np.imag(FT[minnoi:maxnoi]))
    return imagnoisestdev


def IFTEnvelope(chargestatesrCalc, expandedspanCalc, submassCalc, ftxCalc, ftspacingCalc, FTCalc, paddedxnewCalc,
                xnewCalc, ynewCalc, ydataCalc, realnoisestdevCalc, imagnoisestdevCalc):
    omegafinal = expandedspanCalc / submassCalc * 2
    ABIFT = []
    ABIFTnoise = []
    ABIFTmax = []
    ABIFTmaxfinal = 0
    ABIFTintegral = 0
    msintegral = sum(ydataCalc)
    print(msintegral)
    for i in range(0, len(chargestatesrCalc)):
        freqmax = (chargestatesrCalc[i] + 1 / 2) * omegafinal
        freqmin = (chargestatesrCalc[i] - 1 / 2) * omegafinal
        condition = np.logical_and((ftxCalc / ftspacingCalc) > freqmin, (ftxCalc / ftspacingCalc) < freqmax)
        csdata = np.extract(condition, FTCalc)  # extracts the FFT data from the FFT spectrum that are within 1/2
        # the peak spacing of each maximum
        extlen = np.size(csdata)
        leftzeroes = np.ceil(freqmin)
        rightzeroes = np.size(FTCalc) - extlen - leftzeroes
        paddedcsdata = np.lib.pad(csdata, (np.int_(leftzeroes), np.int_(rightzeroes)), 'constant',
                                  constant_values=(0, 0))
        IFT = np.fft.ifft(paddedcsdata)
        ABIFT.append(abs(IFT[int((len(IFT)) / 2):]))

        csdatawnoise = []
        paddedcsdatawnoise = []
        IFTwnoise = []
        ABIFTwrannoise = []
        ABIFTwrannoiseT = []

        for k in range(0, 20):
            csdatawnoise.append(csdata + np.random.normal(0, realnoisestdevCalc, len(csdata)) +
                                np.random.normal(0, imagnoisestdevCalc, len(csdata)) * 1j)
            paddedcsdatawnoise.append(np.lib.pad(csdatawnoise[k], (np.int_(leftzeroes), np.int_(rightzeroes)),
                                                 'constant', constant_values=(0, 0)))
            IFTwnoise.append(np.fft.ifft(paddedcsdatawnoise[k]))
            ABIFTwrannoise.append(abs(IFTwnoise[k][int(len(IFTwnoise[k]) / 2):]))

        tempnoise = []

        ABIFTwrannoiseT = np.transpose(ABIFTwrannoise)

        for m in range(0, len(ABIFTwrannoise[0])):

            if m % 50 == 0:
                tempnoise.append(np.std(ABIFTwrannoiseT[m][:]))

        ABIFTnoise.append(tempnoise)

    ############### Normalization of the IFFT Data ##########################
    for i in range(0, len(chargestatesrCalc)):
        ABIFTmax.append(max(ABIFT[i]))

    ABIFTmaxfinal = max(ABIFTmax)

    for i in range(0, len(chargestatesrCalc)):
        for j in range(0, len(ABIFT[i])):
            ABIFT[i][j] = ABIFT[i][j] / ABIFTmaxfinal

    for i in range(0, len(chargestatesrCalc)):
        ABIFTintegral += sum(ABIFT[i])
    print(ABIFTintegral)
    for i in range(0, len(chargestatesrCalc)):
        for j in range(0, len(ABIFT[i])):
            ABIFT[i][j] = ABIFT[i][j] / ABIFTintegral * msintegral
    ############### Normalization of the IFFT Data ##########################

    return ABIFT, ABIFTnoise, ABIFTintegral, msintegral


def analyzerCalc(self, inputcharge, inputbasemass, inputsubmass):
    name2 = QtGui.QFileDialog.getOpenFileName()
    namestr2 = str(name2[0])
    rawdata = open(name2[0], "r")
    x, y = np.loadtxt(name2[0],
                      unpack=True,
                      delimiter=',')
    avgmz = 0
    varmz = 0
    totarea = np.sum(y)
    for i in range(0, len(y)):
        avgmz += x[i] * y[i]
    avgmz = avgmz / totarea
    avgtotmass = avgmz * inputcharge
    avgtotsubunits = (avgtotmass - inputbasemass) / inputsubmass
    for i in range(0, len(y)):
        varmz += (x[i] ** 2) * y[i]
    varmz = varmz / totarea - avgmz ** 2
    stdevmz = np.sqrt(varmz)
    stdevlip = (inputcharge * stdevmz / inputsubmass)
    avgtotsubunitsdis = str(avgtotsubunits)
    stdevlipdis = str(stdevlip)
    return avgtotsubunitsdis, stdevlipdis


def FFTFilter(expandedspanCalc, submassCalc, ZFreqData, ftxCalc, ftspacingCalc, FTCalc, chargestatesrCalc, OTnum):
    omegafinal = expandedspanCalc / submassCalc * 2
    ABIFT = []
    for i in range(0, len(chargestatesrCalc)):  ###here we include all the harmonics, including in calculating noise
        # bands, so we'll have to add a new loop just to calculate single-harmonic envelopes with and without noise
        freqmax = (chargestatesrCalc[i] + 1 / 2) * omegafinal
        freqmin = (chargestatesrCalc[i] - 1 / 2) * omegafinal
        paddedcsdata = np.zeros(np.size(FTCalc))

        for j in range(1, OTnum):
            condition = np.logical_and((ftxCalc / ftspacingCalc) > j * freqmin, (ftxCalc / ftspacingCalc) < j * freqmax)
            csdata = 2 * np.extract(condition, FTCalc)  # extracts the FFT data from the FFT spectrum that are
            #  within 1/2 the peak spacing of each maximum, and the 2 is because the amplitude was halved due to mirroring the data
            extlen = np.size(csdata)
            print(j * freqmax)
            leftzeroes = np.ceil(j * freqmin)
            rightzeroes = np.size(FTCalc) - extlen - leftzeroes
            paddedcsdata = np.add(paddedcsdata, np.lib.pad(csdata, (np.int_(leftzeroes), np.int_(rightzeroes)),
                                                           'constant', constant_values=(0, 0)))

        IFT = np.fft.ifft(paddedcsdata)
        ABIFT.append(IFT[int((len(IFT)) / 2):])
    reconst = np.zeros(len(ABIFT[1]))
    baseline = np.zeros(len(ABIFT[1]))
    zerofreqwindow = ZFreqData  # governs how far out, in multiples of the fundamental spacing, data for the
    # zero-frequency (FT baseline) peak should be used
    maxfreq = np.max(ftxCalc) / ftspacingCalc
    conditionleft = ftxCalc / ftspacingCalc < zerofreqwindow * omegafinal  # extract Fourier data near zero
    # frequency to add to reconstructed spectrum
    csdataleft = np.extract(conditionleft, FTCalc)  # extracts the FFT data from the FFT spectrum that are within 1/2
    # the peak spacing of each maximum, and the 2 is because the amplitude was halved due to mirroring the data
    csdataleft[0] = csdataleft[0]
    extlenleft = np.size(csdataleft)
    conditionright = ftxCalc / ftspacingCalc > (maxfreq - zerofreqwindow * omegafinal)
    csdataright = np.extract(conditionright, FTCalc)
    extlenright = np.size(csdataright)
    middlezeroes = np.size(FTCalc) - extlenleft - extlenright
    zerofreqdata = np.append(np.lib.pad(csdataleft, (0, np.int_(middlezeroes)),
                                        'constant', constant_values=(0, 0)), csdataright)

    IFTzero = np.fft.ifft(zerofreqdata)
    baseline = np.add(baseline, IFTzero[int((len(IFTzero)) / 2):])

    for i in range(0, len(chargestatesrCalc)):
        reconst = np.add(reconst, ABIFT[i])

    # self.reconstenv = np.abs(self.reconst + self.baseline)
    reconstspec = np.real(reconst + baseline)
    reconstbaseline = np.real(baseline)
    return reconstspec, reconstbaseline


def SinGauAnalyzerCalc(charge, lipidH, lipidL, subunitM, peakstdev, basemass):
    #######inputs########
    lipidnumber = np.arange(lipidL, lipidH)
    lipidstd = np.std(lipidnumber) * subunitM
    #######inputs########

    ###########uses inputs to assemble spetra graphs################
    datamax = int((lipidH * subunitM + basemass) / charge) + 1000  # makes that max of your data set.  Takes the highest
    # possible amount of subunits divded by the lowest charge, and adds 1000 data points
    data = np.zeros(datamax)  # makes the initial data set
    datax = np.arange(1, datamax + 1)  # makes the x axis in the figure
    centroidE = ((((
                lipidH + lipidL)) / 2) * subunitM) + basemass  # makes the centroid of each charge state distribution, using
    # the middle amunt of lipids as the centroid
    ###########uses inputs to assemble spetra graphs################

    return data, datax, lipidnumber, centroidE, lipidstd


def MulGauAnalyzerCalc(chargeH, chargeL, lipidH, lipidL, subunitM,
                       basemass):  # The Analysis functions for Multiple Gaussian
    # Spectra Assembler

    #######inputs########
    charge = np.arange(chargeL, chargeH + 1)
    lipidnumber = np.arange(lipidL, lipidH)
    lipidstd = np.std(lipidnumber) * subunitM

    #######inputs########

    ###########uses inputs to assemble spetra graphs################
    datamax = int(
        (lipidH * subunitM + basemass) / chargeL) + 1000  # makes that max of your data set.  Takes the highest possible
    # amount of subunits divded by the lowest charge, and adds 1000 data points
    data = np.zeros(datamax)  # makes the initial data set
    datax = np.arange(1, datamax + 1)  # makes the x axis in the figure
    centroidT = ((((lipidH + lipidL) / 2) * subunitM) + basemass) / (
                (chargeH + chargeL) / 2)  # makes the centroid of the overall
    # Gaussian.  Takes the middle point of the middle charge state, and calls it that.
    centroidE = ((((lipidH + lipidL)) / 2) * subunitM) + basemass
    return charge, lipidnumber, lipidstd, centroidT, data, datax, centroidE


def analyzer2(y1, y2, y3):
    ############### normalize FT data ##############################################
    AvNumLil = np.zeros(len(y3))
    AvNumLilfinal = np.zeros(len(y3))
    y1new = np.zeros(len(y3))
    y2new = np.zeros(len(y3))
    DevFinal = []
    for i in range(1, 101):
        for j in range(0, len(y3)):
            AvNumLil[j] = ((i * y1[j] + (100 - i) * y2[j]) / 100)
        # AvNumLilmax = np.max(AvNumLil)
        # for j in range (0, len(AvNumLil)):
        # AvNumLil[j] = AvNumLil[j]/AvNumLilmax
        Dev = 0
        for j in range(0, len(y3)):
            Dev += (AvNumLil[j] - y3[j]) ** 2
        DevFinal.append(Dev)

    Min = np.min(DevFinal)
    for i in range(0, len(DevFinal)):
        if DevFinal[i] == Min:
            print(i)
            for j in range(0, len(y1)):
                y1new[j] = i * y1[j]
                y2new[j] = (100 - i) * y2[j]
                AvNumLilfinal[j] = (y1new[j] + y2new[j]) / 100
            AvNumLilfinalMax = np.max(AvNumLilfinal)
            for j in range(0, len(AvNumLilfinal)):
                AvNumLilfinal[j] = AvNumLilfinal[j] / AvNumLilfinalMax
        else:
            continue
    return AvNumLilfinal, DevFinal


def ManCompare(y1, y2, in1, in2):
    AvNumLilfinal = np.zeros(len(y2))
    for i in range(0, len(y1)):
        AvNumLilfinal[i] = (in1 * y1[i] + in2 * y2[i]) / (in1 + in2)
    return AvNumLilfinal


def AverageHarmFun(chargestatesrCalc, expandedspanCalc, submassCalc, ftxCalc, ftspacingCalc, FTCalc, paddedxnewCalc,
                   xnewCalc,
                   ynewCalc, ydataCalc, ov):
    ABIFT = []
    ABIFTmax = []
    ABIFTintegral = 0
    msintegral = sum(ydataCalc)
    print(msintegral)
    for k in range(1, ov + 1):
        ABIFTtemp = []
        for i in range(0, len(chargestatesrCalc)):
            omegafinal = expandedspanCalc / (submassCalc / k) * 2
            freqmax = (chargestatesrCalc[i] + 1 / 2) * omegafinal
            freqmin = (chargestatesrCalc[i] - 1 / 2) * omegafinal
            condition = np.logical_and((ftxCalc / ftspacingCalc) > freqmin, (ftxCalc / ftspacingCalc) < freqmax)
            csdata = np.extract(condition, FTCalc)  # extracts the FFT data from the FFT spectrum that are within 1/2
            # the peak spacing of each maximum
            extlen = np.size(csdata)
            leftzeroes = np.ceil(freqmin)
            rightzeroes = np.size(FTCalc) - extlen - leftzeroes
            paddedcsdata = np.lib.pad(csdata, (np.int_(leftzeroes), np.int_(rightzeroes)), 'constant',
                                      constant_values=(0, 0))
            IFT = np.fft.ifft(paddedcsdata)
            ABIFTtemp.append(abs(IFT[int((len(IFT)) / 2):]))
        ############### Normalization of the IFFT Data ##########################
        for i in range(0, len(chargestatesrCalc)):
            ABIFTmax.append(max(ABIFTtemp[i]))

        ABIFTmaxfinal = max(ABIFTmax)

        for i in range(0, len(chargestatesrCalc)):
            for j in range(0, len(ABIFTtemp[i])):
                ABIFTtemp[i][j] = ABIFTtemp[i][j] / ABIFTmaxfinal

        for i in range(0, len(chargestatesrCalc)):
            ABIFTintegral += sum(ABIFTtemp[i])
        print(ABIFTintegral)
        for i in range(0, len(chargestatesrCalc)):
            for j in range(0, len(ABIFTtemp[i])):
                ABIFTtemp[i][j] = ABIFTtemp[i][j] / ABIFTintegral * msintegral
            ABIFT.append(ABIFTtemp[i])
    ############### Normalization of the IFFT Data ##########################
    print(len(ABIFT))
    Average = FinalAvg(ABIFT, len(chargestatesrCalc), msintegral)
    return Average


def FinalAvg(Average, cs, msintegral):
    if len(Average) > cs:
        print('1')
        AverageNew = []
        for i in range(0, len(Average) - cs):
            Averagetemp = []
            for j in range(0, (len(Average[0]))):
                Averagetemp.append((Average[i][j] + Average[i + cs][j]) / 2)
            AverageNew.append(Averagetemp)
        print(AverageNew)
        return FinalAvg(AverageNew, cs, msintegral)
    else:
        print('2')
        ABIFTmax = []
        ABIFTintegral = 0
        for i in range(0, cs):
            ABIFTmax.append(max(Average[i]))

        ABIFTmaxfinal = max(ABIFTmax)

        for i in range(0, cs):
            for j in range(0, len(Average[i])):
                Average[i][j] = Average[i][j] / ABIFTmaxfinal

        for i in range(0, cs):
            ABIFTintegral += sum(Average[i])
        for i in range(0, cs):
            for j in range(0, len(Average[i])):
                Average[i][j] = Average[i][j] / ABIFTintegral * msintegral
        return Average


def zerocharge(ABIFT, xnew, chargestatesr):
    yfull = []
    for i in range(0, len(ABIFT)):
        yfulltemp = ABIFT[i][0:int(len(xnew))]
        yfull.append(yfulltemp)
    xfull = []
    for i in range(0, len(chargestatesr)):
        if i != len(chargestatesr) + 1:
            xfull.append(xnew)
    xfull = np.array(xfull)

    ######### multiplies each x component by its charge state #############################
    for i in range(0, len(xfull)):
        for j in range(0, len(xfull[0])):
            xfull[i][j] = chargestatesr[i] * xfull[i][j]
    ######### multiplies each x component by its charge state #############################

    ######### interpolates that data to make them equally spaced ##########################
    xnewfinal = []
    ynewfinal = []
    for i in range(0, len(xfull)):
        datainterpolation = interpolate.InterpolatedUnivariateSpline(xfull[i], yfull[i], k=3)
        xnew = np.linspace(np.min(xfull[i]), np.max(xfull[i]), np.size(xfull),
                           endpoint=True)  # number of equally-spaced m/z
        # values is equal to number of
        # m/z values in original data
        xnewfinal.append(xnew)
        ynew = datainterpolation(xnew)
        ynewfinal.append(ynew)
    ######### interpolates that data to make them equally spaced ##########################

    ######### creates range in 10's for the entire zero charge spectrum ######
    xmax = []
    xmin = []
    for i in range(0, len(xnewfinal)):
        xmax.append(xnewfinal[i][-1] / 10)
        xmin.append(xnewfinal[i][0] / 10)
    xmaxfinal = (np.ceil(xmax)) * 10
    xminfinal = (np.floor(xmin)) * 10
    xrange = np.arange((min(xminfinal)), (max(xmaxfinal) + 10), 10)
    ######### creates range in 10's for the entire zero charge spectrum ######

    ########### creates charge specific x axis in 10's for each charge state###########
    xrangespec = []
    xmin = np.floor(xmin) * 10
    xmax = np.ceil(xmax) * 10
    yrangespec = []
    for i in range(0, len(xmin)):
        xrangespectemp = np.arange(xmin[i], xmax[i] + 10, 10)
        xrangespec.append(xrangespectemp)
    for i in range(0, len(xrangespec)):
        datainterpolation = interpolate.InterpolatedUnivariateSpline(xfull[i], yfull[i], k=3)
        yrangespectemp = datainterpolation(xrangespec[i])
        yrangespec.append(yrangespectemp)
    ########### creates charge specific x axis in 10's fir each charge state###########

    ####### places a zero if the data is outside the range of the original charge state specific spectrum #######
    for i in range(0, len(xrangespec)):
        Ytempnum = (max(xrange) - max(xrangespec[i])) / 10
        yzerosright = np.zeros(int(Ytempnum))
        yrangespec[i] = np.append(yrangespec[i], yzerosright)
        Ytempnum2 = (min(xrangespec[i]) - min(xrange)) / 10
        yzerosleft = np.zeros(int(Ytempnum2))
        yrangespec[i] = np.append(yzerosleft, yrangespec[i])
    yfinal = np.zeros(len(yrangespec[0]))
    for i in range(len(yrangespec)):
        yfinal += yrangespec[i]
    ####### places a zero if the data is outside the range of the original charge state specific spectrum #######

    return xrange, yfinal, yrangespec
