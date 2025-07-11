//
// Created by mm96978 on 7/10/2025.
//

#include "UniDecIM_Main.h"

int run_unidec_IM(int argc, char *argv[], Config config) {
    time_t starttime, endtime;
    starttime = time(NULL);

    printf("Opening File: %s\n", config.infile);
    int lines = getfilelengthbin(config.infile, sizeof(float), 3);

    printf("config.length of data: %d\n", lines);

    int *ztab = NULL, *closemind = NULL,
            *closezind = NULL,
            *closecind = NULL,
            *barr = NULL, *closetab = NULL;
    int size[4] = {0, 0, 0, 0};

    float *mzdat = NULL,
            *dtdat = NULL,
            *dataInt = NULL,
            *mzext = NULL,
            *dtext = NULL,
            *masstab = NULL,
            *ccstab = NULL,
            *blur = NULL,
            *newblur = NULL,
            *peakshape = NULL,
            *denom = NULL,
            *deltas = NULL,
            *testmasses = NULL,
            *testmassvals = NULL,
            *testmassCCSavg = NULL,
            *testmassCCSstd = NULL,
            *testmassZavg = NULL,
            *testmassZstd = NULL;

    //Mass Limit File
    int tmlen = 0;
    if (config.mflag == 1) {
        tmlen = getfilelength(config.mfile);
        testmasses = calloc(tmlen, sizeof(float));
        testmassvals = calloc(tmlen, sizeof(float));
        testmassCCSavg = calloc(tmlen, sizeof(float));
        testmassCCSstd = calloc(tmlen, sizeof(float));
        testmassZavg = calloc(tmlen, sizeof(float));
        testmassZstd = calloc(tmlen, sizeof(float));
        readmfile(config.mfile, tmlen, testmasses);
        printf("Read mass file of length: %d\n", tmlen);
    }

    //CCS Constants
    float tempo = 273.15f;
    float ccsconst = 0;
    int ccsconsttype = 0;
    config.temp = config.temp + tempo;
    if (config.press == 0 || config.volt == 0 || config.temp == 0) {
        printf("Either pressure, voltage, or temp is 0. Assuming length is beta (Agilent style)\n");
        if (config.len != 0) {
            ccsconst = 1.0f / config.len;
            ccsconsttype = 1;
        } else {
            printf("Need length parameter to be set to beta");
            exit(1);
        }
    } else {
        double pi = 3.14159265359f;
        float po = 760.f;
        double n = 2.6867774E25;
        double kb = 1.3806488E-23;
        double e = 1.60217657E-19;
		// I know this looks weird, but I had issues with single point floads on the b term.
    	// Switched all to double for consistency and precision.
        float a = (po / config.press) * (config.temp / tempo) * (config.volt / powf(config.len, 2));
        double b = e / sqrt(kb * (double) config.temp) / n;
        double c = sqrt(18 * pi) / 16;
        double d = c * b * (double) a * 1E20;
        ccsconst = (float) d; // Convert back to float for consistency
        printf("CCS Const: %.10e\n", ccsconst);
    }

    if (config.twaveflag > 0) { printf("Ridin' the T-Wave!\n"); }

    //Reading In Data
    mzdat = calloc(lines, sizeof(float));
    dtdat = calloc(lines, sizeof(float));
    dataInt = calloc(lines, sizeof(float));
    readfile3bin(config.infile, lines, mzdat, dtdat, dataInt);

    //Charge States
    int numz = config.endz - config.startz + 1;
    ztab = calloc(numz, sizeof(int));
    for (int i = 0; i < numz; i++) {
        ztab[i] = config.startz + i;
    }

    //Getting Dimension Sizes of Data
    size[0] = GetSize0(mzdat, lines);
    size[1] = GetSize1(dtdat, lines);
    size[2] = numz;
    int totlen = size[0] * size[1] * size[2];
    printf("Dimensions of data: %d mz by %d dt by %d z: %d lines: %d total\n", size[0], size[1], size[2],
           size[0] * size[1], totlen);
    if (totlen > 155E6) { printf("Warning: May exceed system memory capacity\n"); }
    //Extracting mz and dt ranges
    mzext = calloc(size[0], sizeof(float));
    dtext = calloc(size[1], sizeof(float));
    PullXY(mzext, dtext, mzdat, dtdat, size);
    float mzranges[4] = {0, 0, 0, 0};
    mzranges[0] = mzext[0];
    mzranges[1] = mzext[size[0] - 1];
    mzranges[2] = dtext[0];
    mzranges[3] = dtext[size[1] - 1];
    printf("MZ Range: %f to %f\n", mzranges[0], mzranges[1]);
    printf("DT Range: %f to %f\n", mzranges[2], mzranges[3]);

    peakshape = calloc(lines, sizeof(float));
    MakeKernel2D(peakshape, size, mzext, dtext, config.mzsig, config.dtsig, config.psfun, 0);
    printf("Peak Shape Set\n");

    //Filling the mass table
    ccstab = calloc(totlen, sizeof(float));
    int l = size[0] * size[2];
    masstab = calloc(l, sizeof(float));
    barr = calloc(totlen, sizeof(int));
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[2]; j++) {
            masstab[index2D(size[2], i, j)] = calcmass(mzext[i], ztab[j], config.adductmass);
        }
    }
    //Fill CCS table and simultaneously set limit array


#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++) {
            for (int k = 0; k < size[2]; k++) {
                float tempmass = masstab[index2D(size[2], i, k)];
                float tempccs = 0;
                if (config.twaveflag == 0) {
                    tempccs = calcCCS(masstab[index2D(size[2], i, k)], ztab[k], dtext[j], ccsconst, config.hmass,
                                      config.to, ccsconsttype);
                } else if (config.twaveflag == 1) {
                    tempccs = calcCCSTwaveLog(masstab[index2D(size[2], i, k)], ztab[k], dtext[j], config.tcal1,
                                              config.tcal2, config.hmass, config.edc);
                } else if (config.twaveflag == 2) {
                    tempccs = calcCCSTwaveLinear(masstab[index2D(size[2], i, k)], ztab[k], dtext[j], config.tcal1,
                                                 config.tcal2, config.hmass, config.edc);
                } else if (config.twaveflag == 3) {
                    tempccs = calcCCSTwavePower(masstab[index2D(size[2], i, k)], ztab[k], dtext[j], config.tcal1,
                                                config.tcal2, config.hmass, config.edc);
                } else if (config.twaveflag == 4) {
                    tempccs = calcCCSSLIMpoly2(masstab[index2D(size[2], i, k)], ztab[k], dtext[j], config.tcal1,
                                               config.tcal2, config.tcal3, config.hmass, config.edc);
                } else if (config.twaveflag == 5) {
                    tempccs = calcCCSSLIMpoly3(masstab[index2D(size[2], i, k)], ztab[k], dtext[j], config.tcal1,
                                               config.tcal2, config.tcal3, config.tcal4, config.hmass, config.edc);
                } else {
                    TwaveError(config.twaveflag);
                }
                ccstab[index3D(size[1], size[2], i, j, k)] = tempccs;
                float zlimit = nativecharge(tempmass, 0);
                float climit = nativeCCS(tempmass, 0, config.hmass);
                if (tempmass < config.massub && tempmass > config.masslb && tempccs < config.ccsub && tempccs > config.
                    ccslb && (float) ztab[k] < zlimit + config.nativezub && (float) ztab[k] > zlimit + config.nativezlb
                    && tempccs < climit + config.nativeccsub && tempccs > climit + config.nativeccslb) {
                    if (config.mflag == 1) {
                        float testmassclose = nearfastval(testmasses, tempmass, tmlen);
                        if (fabsf(testmassclose - tempmass) <= config.mtabsig) {
                            barr[index3D(size[1], size[2], i, j, k)] = 1;
                        } else { barr[index3D(size[1], size[2], i, j, k)] = 0; }
                    } else {
                        barr[index3D(size[1], size[2], i, j, k)] = 1;
                    }
                } else {
                    barr[index3D(size[1], size[2], i, j, k)] = 0;
                }
            }
        }
    }

    if (config.manualflag == 1) {
        ManualAssign_IM(config.manualfile, size, mzdat, dtdat, ztab, barr);
    }


    //Working out the Blur
    int mlength = 2 * (int) config.msig + 1;
    int zlength = 2 * (int) config.zsig + 1;
    int clength = 2 * (int) config.csig + 1;
    if (config.csig < 0) { clength = 1; }
    int numclose = mlength * zlength * clength;
    size[3] = numclose;
    closemind = calloc(numclose, sizeof(int));
    closezind = calloc(numclose, sizeof(int));
    closecind = calloc(numclose, sizeof(int));
    int n2 = numclose * totlen;
    closetab = calloc(n2, sizeof(int));

#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < mlength; i++) {
        for (int j = 0; j < zlength; j++) {
            for (int k = 0; k < clength; k++) {
                int index = index3D(zlength, clength, i, j, k);
                closemind[index] = (int) config.msig - i;
                closezind[index] = (int) config.zsig - j;
                closecind[index] = (int) config.csig - k;
                //printf("m z c: %d %d %d %d\n",index,config.msig-i,config.zsig-j,config.csig-k);
            }
        }
    }


    printf("Number of Blurs: %d\n", numclose);


#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++) {
            for (int k = 0; k < size[2]; k++) {
                for (int n = 0; n < size[3]; n++) {
                    float ccspt = ccstab[index3D(size[1], size[2], i, j, k)];
                    ccspt = ccspt + (float) closecind[n] * config.ccsbins;
                    int indz = k + closezind[n];
                    int indm = 0;
                    int indc = 0;
                    int badflag = 0;
                    if (indz < 0 || indz >= numz) { badflag = 1; } else {
                        float point = calcmz(masstab[index2D(size[2], i, k)] + (float) closemind[n] * config.molig,
                                             ztab[k] + closezind[n], config.adductmass);
                        if (point < mzranges[0] || point > mzranges[1]) { badflag = 1; } else {
                            float point2 = 0;
                            indm = nearfast(mzext, point, size[0]);
                            if (config.twaveflag ==
                                0) {
                                point2 =
                                        calcDt(point, ztab[k] + closezind[n], ccspt, ccsconst, config.hmass, config.to,
                                               ccsconsttype);
                            } else if (config.twaveflag ==
                                       1) {
                                point2 =
                                        calcDtTwaveLog(point, ztab[k] + closezind[n], ccspt, config.tcal1, config.tcal2,
                                                       config.hmass, config.edc);
                            } else if (config.twaveflag ==
                                       2) {
                                point2 =
                                        calcDtTwaveLinear(point, ztab[k] + closezind[n], ccspt, config.tcal1,
                                                          config.tcal2, config.hmass, config.edc);
                            } else if (config.twaveflag ==
                                       3) {
                                point2 =
                                        calcDtTwavePower(point, ztab[k] + closezind[n], ccspt, config.tcal1,
                                                         config.tcal2, config.hmass, config.edc);
                            } else if (config.twaveflag == 4) {
                                point2 = calcDtSLIMpoly2(point, ztab[k] + closezind[n], ccspt, config.tcal1,
                                                         config.tcal2, config.tcal3, config.hmass, config.edc);
                            } else if (config.twaveflag == 5) {
                                point2 = calcDtSLIMpoly3(point, ztab[k] + closezind[n], ccspt, config.tcal1,
                                                         config.tcal2, config.tcal3, config.tcal4, config.hmass,
                                                         config.edc);
                            } else { TwaveError(config.twaveflag); }
                            if ((point2 < mzranges[2] || point2 > mzranges[3]) && config.csig >=
                                0) { badflag = 1; } else {
                                // printf("Point2: %f\n", point2);
                                if (point2 > 0) { indc = nearfast(dtext, point2, size[1]); } else { badflag = 1; }
                            }
                        }
                    }
                    int index;
                    if (badflag == 1) { index = -1; } else { index = index3D(size[1], size[2], indm, indc, indz); }
                    int index2 = index4D(size, i, j, k, n);
                    closetab[index2] = index;
                }
            }
        }
    }
    printf("Finished Blur\n");


    //Setting Up the Iteration
    blur = calloc(totlen, sizeof(float));
    newblur = calloc(totlen, sizeof(float));
    //memset(barr,1,totlen);
#pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < size[0]; i++) {
        for (int j = 0; j < size[1]; j++) {
            for (int k = 0; k < size[2]; k++) {
                blur[index3D(size[1], size[2], i, j, k)] = dataInt[index2D(size[1], i, j)];
            }
        }
    }

    if (config.intthresh > 0) { KillB_IM(dataInt, barr, size, config.intthresh); }

    denom = calloc(lines, sizeof(float));
    deltas = calloc(lines, sizeof(float));
    printf("Iterating: \n");

    //Iterating
    for (int m = 0; m < config.numit; m++) {
        blur_it_IM(size, blur, newblur, closetab, barr, config.csig);

        sumdeltas(size, deltas, newblur);
        fftconvolve2D(denom, deltas, peakshape, size);
        bayes(size, denom, blur, newblur, dataInt, barr);
        if (config.numit < 10) { printf("Iteration: %d\n", m); } else {
            if (m % (config.numit / 10) == 0) { printf("Iteration: %d\n", m); }
        }
        //printf("Iteration: %d\n",m);
    }

    //Writing outputs

    //Write the Fit config.to the Experimental Data
    sumdeltas(size, deltas, blur);
    fftconvolve2D(denom, deltas, peakshape, size);

    float denommax = Max(denom, lines);
    float datamax = Max(dataInt, lines);
    Normalize(lines, dataInt, datamax);
    Normalize(lines, denom, denommax);
    ApplyCutoff(denom, 0, lines);
    //printf("maxes: %f %f\n",denommax,datamax);
    char *suffixfit = "fitdat";
    write1D(config.outfile, suffixfit, denom, lines);
    float error = errfun(lines, dataInt, denom);

    //Convolve the deltas with the peak shape and put that inconfig.to newblur
    printf("Convolving...");
    convolve3D(newblur, blur, peakshape, size);
    printf("Done\n");

    //Get Max and Min mass and ccs values
    float ranges[4] = {0, 0, 0, 0};
    ranges[0] = config.masslb;
    ranges[1] = config.massub;
    ranges[2] = config.ccslb;
    ranges[3] = config.ccsub;
    if (config.fixedmassaxis == 0) { getranges(size, newblur, masstab, ccstab, ranges, barr); }
    printf("Mass Range: %f to %f  Mass Bins: %f Da\n", ranges[0], ranges[1], config.massbins);
    printf("CCS Range: %f to %f CCS Bins: %f\n", ranges[2], ranges[3], config.ccsbins);

    int maaxle = 1 + (int) ((ranges[1] - ranges[0]) / config.massbins);
    int ccaxle = 1 + (int) ((ranges[3] - ranges[2]) / config.ccsbins);
    int newsize[3] = {0, 0, 0};
    newsize[0] = maaxle;
    newsize[1] = ccaxle;
    newsize[2] = numz;
    printf("Dimensions of Output: %d m by %d ccs by %d z: %d total\n", newsize[0], newsize[1], newsize[2],
           newsize[0] * newsize[1] * newsize[2]);
    float *massaxis = NULL, *massvals = NULL, *ccsaxis = NULL, *ccsvals = NULL, *newgrid = NULL;
    if (newsize[0] * newsize[1] * newsize[2] * sizeof(float) > 2000000) {
        printf("Warning: May exceed system memory capacity");
    }

    massaxis = calloc(maaxle, sizeof(float));
    massvals = calloc(maaxle, sizeof(float));
    ccsaxis = calloc(ccaxle, sizeof(float));
    ccsvals = calloc(ccaxle, sizeof(float));
    int lnew = newsize[0] * newsize[1] * newsize[2];
    newgrid = calloc(lnew, sizeof(float));
    memset(newgrid, 0, lnew);
    makeaxis(massaxis, maaxle, ranges[0], config.massbins);
    makeaxis(ccsaxis, ccaxle, ranges[2], config.ccsbins);

    //Transform the Grid config.to Mass CCS space
    if (config.poolflag == 0) {
        //#pragma omp parallel for private (i,j,k,tempmass,tempccs,indm,indc), schedule(dynamic)
        for (int i = 0; i < size[0]; i++) {
            for (int j = 0; j < size[1]; j++) {
                for (int k = 0; k < size[2]; k++) {
                    float tempmass = masstab[index2D(size[2], i, k)];
                    float tempccs = ccstab[index3D(size[1], size[2], i, j, k)];
                    if (tempmass > massaxis[0] && tempmass < massaxis[newsize[0] - 1] && tempccs > ccsaxis[0] && tempccs
                        < ccsaxis[newsize[1] - 1]) {
                        int indm = nearfast(massaxis, tempmass, newsize[0]);
                        int indc = nearfast(ccsaxis, tempccs, newsize[1]);

                        float val;
                        if (config.rawflag == 0) {
                            val = newblur[index3D(size[1], size[2], i, j, k)];
                        } else { val = blur[index3D(size[1], size[2], i, j, k)]; }

                        int indm2;
                        int indc2;
                        if (massaxis[indm] < tempmass) { indm2 = indm + 1; } else { indm2 = indm - 1; }
                        if (ccsaxis[indc] < tempccs) { indc2 = indc + 1; } else { indc2 = indc - 1; }
                        if (indm2 >= 0 && indm2 < maaxle && indc2 >= 0 && indc2 < ccaxle) {
                            float interposm = LinearInterpolatePosition(massaxis[indm], massaxis[indm2], tempmass);
                            float interposc = LinearInterpolatePosition(ccsaxis[indc], ccsaxis[indc2], tempccs);
                            float val1 = (1 - interposm) * (1 - interposc) * val;
                            float val2 = (1 - interposm) * (interposc) * val;
                            float val3 = (interposm) * (1 - interposc) * val;
                            float val4 = (interposm) * (interposc) * val;
                            newgrid[index3D(newsize[1], newsize[2], indm, indc, k)] += val1;
                            newgrid[index3D(newsize[1], newsize[2], indm, indc2, k)] += val2;
                            newgrid[index3D(newsize[1], newsize[2], indm2, indc, k)] += val3;
                            newgrid[index3D(newsize[1], newsize[2], indm2, indc2, k)] += val4;
                        }
                    }
                }
            }
        }
        printf("Tranformed m/z and dt to mass and ccs by Integration\n");
    } else {
        // #pragma omp parallel for private (i,j,k,tempmass,tempccs,indm,indc), schedule(dynamic)
        for (int i = 0; i < newsize[0]; i++) {
            for (int j = 0; j < newsize[1]; j++) {
                for (int k = 0; k < newsize[2]; k++) {
                    float tempmass = massaxis[i];
                    float tempmz = calcmz(tempmass, ztab[k], config.adductmass);
                    float tempccs = ccsaxis[j];
                    float tempdt, endval = 0;

                    if (config.twaveflag ==
                        0) {
                        tempdt = calcDt(tempmass, ztab[k], tempccs, ccsconst, config.hmass, config.to, ccsconsttype);
                    } else if (config.twaveflag ==
                               1) {
                        tempdt =
                                calcDtTwaveLog(tempmass, ztab[k], tempccs, config.tcal1, config.tcal2, config.hmass,
                                               config.edc);
                    } else if (config.twaveflag ==
                               2) {
                        tempdt =
                                calcDtTwaveLinear(tempmass, ztab[k], tempccs, config.tcal1, config.tcal2, config.hmass,
                                                  config.edc);
                    } else if (config.twaveflag ==
                               3) {
                        tempdt = calcDtTwavePower(tempmass, ztab[k], tempccs, config.tcal1, config.tcal2, config.hmass,
                                                  config.edc);
                    } else if (config.twaveflag == 4) {
                        tempdt = calcDtSLIMpoly2(tempmass, ztab[k], tempccs, config.tcal1, config.tcal2, config.tcal3,
                                                 config.hmass, config.edc);
                    } else if (config.twaveflag == 5) {
                        tempdt = calcDtSLIMpoly3(tempmass, ztab[k], tempccs, config.tcal1, config.tcal2, config.tcal3,
                                                 config.tcal4, config.hmass, config.edc);
                    } else { TwaveError(config.twaveflag); }

                    if (tempmz > mzext[1] && tempmz < mzext[size[0] - 2] && tempdt > dtext[1] && tempdt < dtext[
                            size[1] - 2]) {
                        int indm = nearfast(mzext, tempmz, size[0]);
                        int indc = nearfast(dtext, tempdt, size[1]);
                        // printf("indm: %d indc: %d k: %d tempmass: %f tempmz: %f\n", indm, indc, k, tempmass, tempmz);
                        endval = bicubicinterpolation(size, indm, indc, k, mzext, dtext, tempmz, tempdt, config.rawflag,
                                                      newblur, blur); //bicubic interpolation
                        //endval=bilinearinterpolation(size,indm,indc,k,mzext,dtext,tempmz,tempdt,config.rawflag,newblur,blur);//bilinear interpolation
                        //endval=newblur[index3D(size[1],size[2],indm,indc,k)]; //(Nearest Neighbor)
                    }
                    newgrid[index3D(newsize[1], newsize[2], i, j, k)] = endval;
                }
            }
        }
        printf("Transformed m/z and dt to mass and ccs by Interpolation\n");
    }
    //Sum all config.to mass axis
    Sum1D(newsize, massvals, newgrid, 0);
    Sum1D(newsize, ccsvals, newgrid, 1);
    char *suffixmass = "mass";
    write2D(config.outfile, suffixmass, massaxis, massvals, maaxle);
    char *suffixccs = "ccs";
    write2D(config.outfile, suffixccs, ccsaxis, ccsvals, ccaxle);

    //Sum across charge and CCS
    float *massccsgrid = NULL, *masszgrid = NULL, *ccszgrid = NULL; //,*ccsaxis=NULL,*ccsvals=NULL,*newgrid=NULL;
    int l1 = maaxle * ccaxle;
    int l2 = maaxle * numz;
    int l3 = ccaxle * numz;
    massccsgrid = calloc(l1, sizeof(float));
    masszgrid = calloc(l2, sizeof(float));
    ccszgrid = calloc(l3, sizeof(float));
    Sum2D(newsize, massccsgrid, newgrid, 2);
    Sum2D(newsize, masszgrid, newgrid, 1);
    Sum2D(newsize, ccszgrid, newgrid, 0);
    char *suffix3 = "massccs";
    write1D(config.outfile, suffix3, massccsgrid, maaxle * ccaxle);
    char *suffix4 = "massgrid";
    write1D(config.outfile, suffix4, masszgrid, maaxle * numz);
    char *suffix6 = "ccsz";
    write1D(config.outfile, suffix6, ccszgrid, ccaxle * numz);

    if (config.mflag == 1) {
        MFileInt(tmlen, massaxis, massvals, testmasses, testmassvals, maaxle);
        MFileCCS(tmlen, massaxis, ccsaxis, massccsgrid, testmasses, testmassCCSavg, testmassCCSstd, maaxle, ccaxle);
        MFileZ(tmlen, massaxis, ztab, masszgrid, testmasses, testmassZavg, testmassZstd, maaxle, numz);
        char *suffixmfileint = "mfileresults";
        writemfileres(config.outfile, suffixmfileint, testmasses, testmassvals, testmassCCSavg, testmassCCSstd,
                      testmassZavg, testmassZstd, tmlen);
    }

    //Write out specific charge slices
    char *suffix5 = "zout";
    if (config.zout > 0) {
        writezslice(newsize, config.outfile, suffix5, ztab, newgrid, nearint(ztab, config.zout, numz));
    }
    if (config.zout < 0) {
        for (int k = 0; k < numz; k++) {
            writezslice(newsize, config.outfile, suffix5, ztab, newgrid, k);
        }
    }
    char *suffixgrid = "mzgrid";
    writemzgrid(config.outfile, suffixgrid, newblur, size);

    char outstring[500];
    char *suffixerr = "error";
    FILE *out_ptrIM = NULL;
    sprintf(outstring, "%s_%s.txt", config.outfile, suffixerr);
    out_ptrIM = fopen(outstring, "w");
    if (out_ptrIM == 0) {
        printf("Error Opening %s\n", outstring);
        exit(1);
    }
    fprintf(out_ptrIM, "error %f\n", error);
    fclose(out_ptrIM);
    printf("File written to: %s\n", outstring);


    printf("\nR Squared: %f\n", error);

    free(massaxis);
    free(massvals);
    free(ccsaxis);
    free(ccsvals);
    free(newgrid);
    free(barr);
    free(mzdat);
    free(dtdat);
    free(dataInt);
    free(ccstab);
    free(ztab);
    free(masstab);
    free(blur);
    free(newblur);
    free(peakshape);
    free(denom);
    free(deltas);
    free(mzext);
    free(dtext);
    free(closemind);
    free(closezind);
    free(closecind);
    free(closetab);
    free(testmasses);
    free(testmassvals);
    free(testmassCCSavg);
    free(testmassCCSstd);
    free(testmassZavg);
    free(testmassZstd);
    free(massccsgrid);
    free(masszgrid);
    free(ccszgrid);
    endtime = time(NULL);

    printf("\nDone in %ds!\n", (int) difftime(endtime, starttime));
    return 0;
}
