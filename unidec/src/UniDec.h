/*
 * UniDec.h
 *
 *  Created on: 1 Jul 2014
 *      Author: Michael.Marty
 */

/*
#ifndef UNIDEC_H_
#define UNIDEC_H_
#endif*/ /* UNIDEC_H_ */

#ifndef UNIDEC_HEADER
#define UNIDEC_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
// #include <time.h>

#include "udtools.h"
#include "udstruct.h"
#include "udio.h"
#include "udcore.h"
// #include "UD_sg.h"




void TestMassListWindowed(int lengthmz, int numz, char *barr, float *mtab, float nativezub, float nativezlb,
                          float massub, float masslb, int *nztab, float *testmasses, int mfilelen, float mtabsig) {
    unsigned int i, j;
    float testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(auto)
    for (i = 0; i < lengthmz; i++) {
        for (j = 0; j < numz; j++) {
            barr[index2D(numz, i, j)] = 0;
            testmass = mtab[index2D(numz, i, j)];
            nativelimit = nativecharge(testmass, 0);
            if (testmass < massub && testmass > masslb && nztab[j] < nativelimit + nativezub && nztab[j] > nativelimit +
                nativezlb) {
                if (neartest(testmasses, testmass, mfilelen, mtabsig) == 1) {
                    barr[index2D(numz, i, j)] = 1;
                }
            }
        }
    }
}

void TestMassListLimit(int lengthmz, int numz, char *barr, float *mtab, float nativezub, float nativezlb, float massub,
                       float masslb, int *nztab, int *testmasspos, int mfilelen) {
    unsigned int i, j, k;
    float testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(auto)
    for (i = 0; i < lengthmz; i++) {
        for (j = 0; j < numz; j++) {
            barr[index2D(numz, i, j)] = 0;
            testmass = mtab[index2D(numz, i, j)];
            nativelimit = nativecharge(testmass, 0);
            if (testmass < massub && testmass > masslb && nztab[j] < nativelimit + nativezub && nztab[j] > nativelimit +
                nativezlb) {
                for (k = 0; k < mfilelen; k++) {
                    if (testmasspos[index2D(numz, k, j)] == i) {
                        barr[index2D(numz, i, j)] = 1;
                    }
                }
            }
        }
    }
}

void TestMass(int lengthmz, int numz, char *barr, float *mtab, float nativezub, float nativezlb, float massub,
              float masslb, int *nztab) {
    unsigned int i, j;
    float testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(auto)
    for (i = 0; i < lengthmz; i++) {
        for (j = 0; j < numz; j++) {
            testmass = mtab[index2D(numz, i, j)];
            nativelimit = nativecharge(testmass, 0);
            if (testmass < massub && testmass > masslb && nztab[j] < nativelimit + nativezub && nztab[j] > nativelimit +
                nativezlb) {
                barr[index2D(numz, i, j)] = 1;
            } else {
                barr[index2D(numz, i, j)] = 0;
            }
        }
    }
}


void MakeSparseBlur(const int numclose, char *barr, const int *closezind,
                    const int *closemind, const float *mtab, const int *nztab,
                    const float *dataMZ, int *closeind, float *closeval, float *closearray, const Config config) {
    int lengthmz = config.lengthmz;
    int numz = config.numz;
    float molig = config.molig;
    float adductmass = config.adductmass;

#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++) {
        for (int j = 0; j < numz; j++) {
            if (barr[index2D(numz, i, j)] == 1) {
                int num = 0;

                //Reset the threshold if it is zero
                float mzsig = config.mzsig;
                if (mzsig == 0) {
                    int i1 = i - 1;
                    int i2 = i + 1;
                    if (i >= lengthmz - 1) { i2 = i; }
                    if (i == 0) { i1 = i; }
                    mzsig = 2 * fabs(dataMZ[i2] - dataMZ[i1]);
                    if (mzsig > config.massbins || mzsig == 0) { mzsig = config.massbins * 2; }
                }
                float newthreshold = mzsig * 2;

                for (int k = 0; k < numclose; k++) {
                    //Find the z index and test if it is within the charge range
                    int indz = (int) (j + closezind[k]);
                    if (indz < 0 || indz >= numz || (nztab[j] + closezind[k]) == 0) {
                        closeind[index3D(numz, numclose, i, j, k)] = -1;
                        closearray[index3D(numz, numclose, i, j, k)] = 0;
                    } else {
                        //Find the nearest m/z value and test if it is close enough to the predicted one and within appropriate ranges
                        float point = (float) ((mtab[index2D(numz, i, j)] + closemind[k] * molig + adductmass * (float)
                                                (nztab[j] + closezind[k])) / (float) (nztab[j] + closezind[k]));
                        if (point < dataMZ[0] - newthreshold || point > dataMZ[lengthmz - 1] + newthreshold) {
                            closeind[index3D(numz, numclose, i, j, k)] = -1;
                            closearray[index3D(numz, numclose, i, j, k)] = 0;
                        } else {
                            int ind = nearfast(dataMZ, point, lengthmz);
                            float closepoint = dataMZ[ind];
                            int newind = index2D(numz, ind, indz);
                            if (barr[newind] == 1 && fabs(point - closepoint) < newthreshold) {
                                closeind[index3D(numz, numclose, i, j, k)] = newind;
                                closearray[index3D(numz, numclose, i, j, k)] = closeval[k] * mzpeakshape(
                                                                                   point, closepoint, mzsig,
                                                                                   config.psfun);
                                num += 1;
                            } else {
                                closeind[index3D(numz, numclose, i, j, k)] = -1;
                                closearray[index3D(numz, numclose, i, j, k)] = 0;
                            }
                            //printf("%d %d %d %f %f %d\n", i, j, k, point, closepoint, closeind[index3D(numz, numclose, i, j, k)]);
                        }
                    }
                }
                if (num < 2 && config.manualflag == 0) { barr[index2D(numz, i, j)] = 0; }
                // printf("%d %d \n", i, j);}
            } else {
                for (int k = 0; k < numclose; k++) {
                    closeind[index3D(numz, numclose, i, j, k)] = -1;
                    closearray[index3D(numz, numclose, i, j, k)] = 0;
                }
            }
        }
    }
}


void MakePeakShape2D(int lengthmz, int maxlength, int *starttab, int *endtab, float *dataMZ, float mzsig, int psfun,
                     int speedyflag, float *mzdist, float *rmzdist, int makereverse) {
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++) {
        int start = starttab[i];
        int end = endtab[i];
        for (int j = start; j <= end; j++) {
            int j2 = fixk(j, lengthmz);
            mzdist[index2D(maxlength, i, j2 - start)] = mzpeakshape(dataMZ[i], dataMZ[j2], mzsig, psfun);
            if (makereverse == 1) {
                rmzdist[index2D(maxlength, i, j2 - start)] = mzpeakshape(dataMZ[j2], dataMZ[i], mzsig, psfun);
            }
        }
    }
}

void MakePeakShape1D(float *dataMZ, float threshold, int lengthmz, int speedyflag, float mzsig, int psfun,
                     float *mzdist, float *rmzdist, int makereverse) {
    float binsize = dataMZ[1] - dataMZ[0];
    float newrange = threshold / binsize;
    int n;
    for (n = (int) -newrange; n < (int) newrange; n++) {
        mzdist[indexmod(lengthmz, 0, n)] = mzpeakshape(0, n * binsize, mzsig, psfun);
        if (makereverse == 1) { rmzdist[indexmod(lengthmz, 0, n)] = mzpeakshape(n * binsize, 0, mzsig, psfun); }
    }
    printf("\nNotice: Assuming linearized data. \n\n");
}


float Reconvolve(const int lengthmz, const int numz, const int maxlength, const int *starttab, const int *endtab,
                 const float *mzdist, const float *blur, float *newblur, const int speedyflag, const char *barr) {
    float newblurmax = 0;
    if (speedyflag == 0) {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            for (int j = 0; j < numz; j++) {
                float cv = 0;
                if (barr[index2D(numz, i, j)] == 1) {
                    for (int k = starttab[i]; k <= endtab[i]; k++) {
                        int k2 = fixk(k, lengthmz);
                        if (blur[index2D(numz, k2, j)] != 0) {
                            int start = starttab[k2];
                            cv += blur[index2D(numz, k2, j)] * mzdist[index2D(maxlength, k2, i - start)];
                        }
                    }
                }
                newblur[index2D(numz, i, j)] = cv;
                if (cv > newblurmax) {
                    newblurmax = cv;
                }
            }
        }
    } else {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            for (int j = 0; j < numz; j++) {
                float cv = 0;
                if (barr[index2D(numz, i, j)] == 1) {
                    for (int k = starttab[i]; k <= endtab[i]; k++) {
                        if (blur[index2D(numz, k, j)] != 0) {
                            cv += blur[index2D(numz, k, j)] * mzdist[indexmod(lengthmz, k, i)];
                        }
                    }
                }
                newblur[index2D(numz, i, j)] = cv;
                if (cv > newblurmax) {
                    newblurmax = cv;
                }
            }
        }
    }
    return newblurmax;
}


void IntegrateTransform(const int lengthmz, const int numz, const float *mtab, float massmax, float massmin,
                        const int maaxle, float *massaxis, float *massaxisval, const float *blur, float *massgrid) {
    for (int i = 0; i < lengthmz; i++) {
        for (int j = 0; j < numz; j++) {
            float testmass = mtab[index2D(numz, i, j)];
            if (testmass < massmax && testmass > massmin) {
                int index = nearfast(massaxis, testmass, maaxle);
                float newval = blur[index2D(numz, i, j)];
                if (massaxis[index] == testmass) {
                    massaxisval[index] += newval;
                    massgrid[index2D(numz, index, j)] += newval;
                }

                if (massaxis[index] < testmass && index < maaxle - 2) {
                    int index2 = index + 1;
                    float interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
                    massaxisval[index] += (1.0 - interpos) * newval;
                    massgrid[index2D(numz, index, j)] += (1.0 - interpos) * newval;
                    massaxisval[index2] += (interpos) * newval;
                    massgrid[index2D(numz, index2, j)] += (interpos) * newval;
                }

                if (massaxis[index] > testmass && index > 0) {
                    int index2 = index - 1;
                    float interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
                    massaxisval[index] += (1 - interpos) * newval;
                    massgrid[index2D(numz, index, j)] += (1 - interpos) * newval;
                    massaxisval[index2] += (interpos) * newval;
                    massgrid[index2D(numz, index2, j)] += (interpos) * newval;
                }
            }
        }
    }
}

void InterpolateTransform(const int maaxle, const int numz, const int lengthmz, const int *nztab, float *massaxis,
                          const float adductmass, const float *dataMZ, float *massgrid, float *massaxisval,
                          const float *blur) {
    float startmzval = dataMZ[0];
    float endmzval = dataMZ[lengthmz - 1];
    //#pragma omp parallel for schedule(auto)
    for (int i = 0; i < maaxle; i++) {
        float val = 0;

        for (int j = 0; j < numz; j++) {
            float newval = 0;
            float mztest = (massaxis[i] + nztab[j] * adductmass) / (float) nztab[j];

            if (mztest > startmzval && mztest < endmzval) {
                int index = nearfast(dataMZ, mztest, lengthmz);
                int index2 = index;
                if (dataMZ[index] == mztest) {
                    newval = blur[index2D(numz, index, j)];
                    val += newval;
                    massgrid[index2D(numz, i, j)] = newval;
                } else {
                    if (dataMZ[index] > mztest && index > 1 && index < lengthmz - 1) {
                        index2 = index;
                        index = index - 1;
                    } else if (dataMZ[index] < mztest && index < lengthmz - 2 && index > 0) {
                        index2 = index + 1;
                    }
                    if (index2 > index && (dataMZ[index2] - dataMZ[index]) != 0) {
                        float mu = (mztest - dataMZ[index]) / (dataMZ[index2] - dataMZ[index]);
                        float y0 = blur[index2D(numz, index - 1, j)];
                        float y1 = blur[index2D(numz, index, j)];
                        float y2 = blur[index2D(numz, index2, j)];
                        float y3 = blur[index2D(numz, index2 + 1, j)];
                        newval = clip(CubicInterpolate(y0, y1, y2, y3, mu), 0);
                        //newval=CRSplineInterpolate(y0,y1,y2,y3,mu);
                        //newval=LinearInterpolate(y1,y2,mu);
                        val += newval;
                        massgrid[index2D(numz, i, j)] = newval;
                    }
                }
            }
        }
        massaxisval[i] = val;
    }
}

void SmartTransform(const int maaxle, const int numz, const int lengthmz, const int *nztab, float *massaxis,
                    const float adductmass, const float *dataMZ, float *massgrid, float *massaxisval,
                    const float *blur) {
    float startmzval = dataMZ[0];
    float endmzval = dataMZ[lengthmz - 1];
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < maaxle; i++) {
        float val = 0;
        for (int j = 0; j < numz; j++) {
            int z = nztab[j];
            float mtest = massaxis[i];
            float mztest = calcmz(mtest, z, adductmass);

            float mzlower;
            float mlower;
            if (i > 0) {
                mlower = massaxis[i - 1];
                mzlower = calcmz(mlower, z, adductmass);
            } else {
                mzlower = mztest;
                mlower = mtest;
            }

            float mzupper;
            float mupper;
            if (i < maaxle - 1) {
                mupper = massaxis[i + 1];
                mzupper = calcmz(mupper, z, adductmass);
            } else {
                mzupper = mztest;
                mupper = mtest;
            }

            float newval = 0;
            if (mzupper > startmzval && mzlower < endmzval) {
                int index = nearfast(dataMZ, mztest, lengthmz);
                int index1 = nearfast(dataMZ, mzlower, lengthmz);
                int index2 = nearfast(dataMZ, mzupper, lengthmz);
                float imz = dataMZ[index];

                if (index2 - index1 < 5) {
                    if (imz == mztest) {
                        newval = clip(blur[index2D(numz, index, j)], 0);
                        val += newval;
                        massgrid[index2D(numz, i, j)] = newval;
                    } else {
                        int edge = 0;
                        index2 = index;
                        if (imz > mztest) {
                            index = index - 1;
                        } else if (imz < mztest) {
                            index2 = index + 1;
                        }

                        if (index < 1 || index2 >= lengthmz - 1) {
                            edge = 1;
                        }
                        if (index < 0 || index2 >= lengthmz) {
                            edge = 2;
                        }


                        if (index2 > index && (dataMZ[index2] - dataMZ[index]) != 0 && edge == 0) {
                            float mu = (mztest - dataMZ[index]) / (dataMZ[index2] - dataMZ[index]);
                            float y0 = blur[index2D(numz, index - 1, j)];
                            float y1 = blur[index2D(numz, index, j)];
                            float y2 = blur[index2D(numz, index2, j)];
                            float y3 = blur[index2D(numz, index2 + 1, j)];
                            newval = clip(CubicInterpolate(y0, y1, y2, y3, mu), 0);
                            //newval=clip(CRSplineInterpolate(y0,y1,y2,y3,mu), 0);
                            //newval=clip(LinearInterpolate(y1,y2,mu),0);
                            val += newval;
                            massgrid[index2D(numz, i, j)] = newval;
                            //printf("0\n");
                        } else if (edge == 1 && (dataMZ[index2] - dataMZ[index]) != 0) {
                            float mu = (mztest - dataMZ[index]) / (dataMZ[index2] - dataMZ[index]);
                            float y1 = blur[index2D(numz, index, j)];
                            float y2 = blur[index2D(numz, index2, j)];
                            newval = clip(LinearInterpolate(y1, y2, mu), 0);
                            val += newval;
                            massgrid[index2D(numz, i, j)] = newval;
                            //printf("1 %d %d %f %f %f\n", index, index2, mztest, imz, newval, massaxis[i]);
                        } else if (edge == 2 && (dataMZ[index2] - dataMZ[index]) != 0) {
                            if (index2 == 0) {
                                index = 0;
                                index2 = 1;
                            }
                            if (index == lengthmz - 1) {
                                index = lengthmz - 1;
                                index2 = lengthmz - 2;
                            }
                            float mu = (mztest - dataMZ[index]) / (dataMZ[index] - dataMZ[index2]);
                            float y1 = blur[index2D(numz, index, j)];
                            float y2 = 0;
                            newval = clip(LinearInterpolate(y1, y2, mu), 0);
                            val += newval;
                            massgrid[index2D(numz, i, j)] = newval;
                            //printf("2\n");
                        }
                    }
                } else {
                    //printf("Integrating\n");
                    float num = 0;
                    for (int k = index1; k < index2 + 1; k++) {
                        float kmz = dataMZ[k];
                        float km = (kmz - adductmass) * z;
                        float scale = 1;

                        if (mztest < kmz && km < mupper) {
                            //scale = LinearInterpolatePosition(mzupper, mztest, kmz);
                            scale = LinearInterpolatePosition(mupper, mtest, km);
                        } else if (kmz < mztest && km > mlower) {
                            //scale = LinearInterpolatePosition(mzlower, mztest, kmz);
                            scale = LinearInterpolatePosition(mlower, mtest, km);
                        } else if (kmz == mztest) { scale = 1; } else { scale = 0; }

                        newval += scale * blur[index2D(numz, k, j)];
                        num += scale;
                    }
                    if (num != 0) { newval /= num; }
                    newval = clip(newval, 0);
                    val += newval;
                    massgrid[index2D(numz, i, j)] = newval;
                }
            }
        }
        massaxisval[i] = val;
    }
}


void ignorezeros(char *barr, const float *dataInt, const int lengthmz, const int numz) {
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++) {
        float val = dataInt[i];
        if (val == 0) {
            for (int j = 0; j < numz; j++) {
                barr[index2D(numz, i, j)] = 0;
            }
        }
    }
}


void KillB(float *I, char *B, float intthresh, int lengthmz, int numz) {
    for (int i = 0; i < lengthmz; i++) {
        for (int j = 0; j < numz; j++) {
            if (I[i] <= intthresh) { B[index2D(numz, i, j)] = 0; }
        }
    }
}


void charge_scaling(float *blur, const int *nztab, const int lengthmz, const int numz) {
    for (int i = 0; i < lengthmz; i++) {
        for (int j = 0; j < numz; j++) {
            int charge = nztab[j];
            float z = (float) charge;
            if (z != 0) { blur[index2D(numz, i, j)] /= z; }
        }
    }
}


void point_smoothing(float *blur, const char *barr, const int lengthmz, const int numz, const int width) {
    float *newblur;
    int l = lengthmz * numz;
    newblur = (float *) calloc(l, sizeof(float));
    if (newblur) {
        size_t newl = (size_t) lengthmz * numz * sizeof(float);
        memcpy(newblur, blur, newl);
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            for (int j = 0; j < numz; j++) {
                if (barr[index2D(numz, i, j)] == 1) {
                    int low = i - width;
                    if (low < 0) { low = 0; }
                    int high = i + width + 1;
                    if (high > lengthmz) { high = lengthmz; }

                    float sum = 0;
                    for (int k = low; k < high; k++) {
                        sum += newblur[index2D(numz, k, j)];
                    }

                    blur[index2D(numz, i, j)] = sum / ((float) 1 + 2 * width);
                }
            }
        }
        free(newblur);
    }
    return;
}

void point_smoothing_peak_width(const int lengthmz, const int numz, const int maxlength, const int *starttab,
                                const int *endtab, const float *mzdist, float *blur, const int speedyflag,
                                const char *barr) {
    float *newblur;
    int ln = lengthmz * numz;
    newblur = (float *) calloc(ln, sizeof(float));
    if (newblur) {
        memcpy(newblur, blur, ln * sizeof(float));
        Reconvolve(lengthmz, numz, maxlength, starttab, endtab, mzdist, newblur, blur, speedyflag, barr);
    }
    return;
}


void softargmax_transposed(float *blur, const int lengthmz, const int numz, const float beta, const char *barr,
                           const int maxlength, const int speedyflag
                           , int *starttab, int *endtab, float *mzdist, const float mzsig) {
    float *newblur, *deltas, *deltas2, *denom, *denom2;
    int newl = lengthmz * numz;
    newblur = (float *) calloc(newl, sizeof(float));
    deltas = (float *) calloc(lengthmz, sizeof(float));
    deltas2 = (float *) calloc(lengthmz, sizeof(float));
    denom = (float *) calloc(lengthmz, sizeof(float));
    denom2 = (float *) calloc(lengthmz, sizeof(float));
    size_t newl2 = (size_t) lengthmz * numz * sizeof(float);
    memcpy(newblur, blur, newl2);

    //Sum deltas
    sum_deltas(lengthmz, numz, blur, barr, deltas);

    if (mzsig != 0) //Convolve with peak shape
    {
        convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, denom, speedyflag);
    } else { memcpy(denom, deltas, sizeof(float) * lengthmz); }

#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz * numz; i++) {
        blur[i] = exp(beta * newblur[i]) - 1;
    }

    //Sum deltas
    sum_deltas(lengthmz, numz, blur, barr, deltas2);

    if (mzsig != 0) //Convolve with peak shape
    {
        convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas2, denom2, speedyflag);
    } else { memcpy(denom2, deltas, sizeof(float) * lengthmz); }

#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++) {
        float factor = 0;
        if (denom2[i] != 0) { factor = denom[i] / denom2[i]; };
        for (int j = 0; j < numz; j++) {
            //blur[index2D(numz, i, j)] -= 1;
            blur[index2D(numz, i, j)] *= factor;
        }
        //printf("%f %f\n", sum1, sum3);
    }
    free(newblur);
    free(deltas);
    free(deltas2);
    free(denom);
    free(denom2);
    return;
}

void softargmax_everything(float *blur, const int lengthmz, const int numz, const float beta) {
    float *newblur;
    int l = lengthmz * numz;
    size_t sizel = (size_t) l * sizeof(float);
    newblur = (float *) calloc(l, sizeof(float));
    if (newblur) {
        memcpy(newblur, blur, sizel);
        //float max1 = 0;
        //float max2 = 0;
        float sum2 = 0;
        float sum1 = 0;
        float min2 = 1000000000000000.0;
        //#pragma omp parallel for schedule(auto), reduction(min:min2), reduction(+:sum1), reduction(+:sum2)
        for (int i = 0; i < lengthmz * numz; i++) {
            float d = newblur[i];
            //if (d > max1) { max1 = d; }
            float e = exp(beta * d);
            //float e = pow(d, -beta);
            //if (e > max2) { max2 = e; }
            blur[i] = e;
            if (e < min2) { min2 = e; }
            sum1 += d;
            sum2 += e;
        }
        float factor = 0;
        //float denom = (max2 - min2);
        //if (denom != 0) { factor = max1 / denom; };
        float denom = (sum2 - min2 * numz * lengthmz);
        if (denom != 0) { factor = sum1 / denom; };
        //if (sum2 != 0) { factor = sum1 / sum2; };
        if (factor > 0) {
#pragma omp parallel for schedule(auto)
            for (int i = 0; i < lengthmz * numz; i++) {
                blur[i] -= min2;
                blur[i] *= factor;
            }
        }
        free(newblur);
    }
    return;
}

void softargmax(float *blur, const int lengthmz, const int numz, const float beta) {
    if (beta < 0) {
        //softargmax_transposed(blur, lengthmz, numz, fabs(beta));
        softargmax_everything(blur, lengthmz, numz, fabs(beta));
        return;
    }

    float *newblur;
    int l = lengthmz * numz;
    size_t sizel = (size_t) l * sizeof(float);
    newblur = (float *) calloc(l, sizeof(float));
    if (newblur) {
        memcpy(newblur, blur, sizel);
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            float sum2 = 0;
            float sum1 = 0;
            float factor = 0;
            float min2 = 1000000000000.0;

            for (int j = 0; j < numz; j++) {
                float d = newblur[index2D(numz, i, j)];
                sum1 += d;

                float e = exp(beta * d);

                if (e < min2) { min2 = e; }

                blur[index2D(numz, i, j)] = e;
                sum2 += e;
            }

            float denom = (sum2 - min2 * numz);
            if (denom != 0) { factor = sum1 / denom; };

            if (factor > 0) {
                for (int j = 0; j < numz; j++) {
                    blur[index2D(numz, i, j)] -= min2;
                    blur[index2D(numz, i, j)] *= factor;
                }
            } else {
                for (int j = 0; j < numz; j++) {
                    blur[index2D(numz, i, j)] = 0;
                }
            }
        }
        free(newblur);
    }
    return;
}


void SetLimits(const Config config, Input *inp) {
    //Determines the indexes of each test mass from mfile in m/z space
    int *testmasspos = (int *) malloc(sizeof(float) * config.mfilelen * config.numz);
    if (testmasspos) {
        if (config.mflag == 1 && config.limitflag == 1) {
            for (int i = 0; i < config.mfilelen; i++) {
                for (int j = 0; j < config.numz; j++) {
                    float mztest = (inp->testmasses[i] + config.adductmass * inp->nztab[j]) / (float) inp->nztab[j];
                    testmasspos[index2D(config.numz, i, j)] = nearfast(inp->dataMZ, mztest, config.lengthmz);
                }
            }
        }

        //If there is a mass file read, it will only allow masses close to those masses within some config.mtabsig window.
        if (config.mflag == 1 && config.limitflag == 0) {
            TestMassListWindowed(config.lengthmz, config.numz, inp->barr, inp->mtab, config.nativezub, config.nativezlb,
                                 config.massub, config.masslb, inp->nztab, inp->testmasses, config.mfilelen,
                                 config.mtabsig);
        }
        //If there is a mass file read and the mass table window (config.mtabsig) is 0, it will only write intensities at the m/z values closest to the m/z values read in from the mfile.
        else if (config.mflag == 1 && config.limitflag == 1) {
            TestMassListLimit(config.lengthmz, config.numz, inp->barr, inp->mtab, config.nativezub, config.nativezlb,
                              config.massub, config.masslb, inp->nztab, testmasspos, config.mfilelen);
        }
        //Normally, write the intensity values if the values fall within the mass upperbound and lower bound
        else {
            TestMass(config.lengthmz, config.numz, inp->barr, inp->mtab, config.nativezub, config.nativezlb,
                     config.massub, config.masslb, inp->nztab);
        }
        free(testmasspos);
    }
}

//Sets the maxlength parameter and the start and end values for the m/z peak shape convolution
//Convolution uses a reflection for the edges, so some care needs to be taken when things are over the edge.
int SetStartsEnds(const Config config, const Input *inp, int *starttab, int *endtab) {
    int maxlength = 1;
    for (int i = 0; i < config.lengthmz; i++) {
        float point = inp->dataMZ[i] - config.psmzthresh;
        int start, end;
        if (point < inp->dataMZ[0] && config.speedyflag == 0) {
            //start = (int)((point - inp->dataMZ[0]) / (inp->dataMZ[1] - inp->dataMZ[0]));
            start = 0 - nearfast(inp->dataMZ, 2 * inp->dataMZ[0] - point, config.lengthmz);
        } else {
            start = nearfast(inp->dataMZ, point, config.lengthmz);
        }
        starttab[i] = start;

        point = inp->dataMZ[i] + config.psmzthresh;
        if (point > inp->dataMZ[config.lengthmz - 1] && config.speedyflag == 0) {
            //end = config.lengthmz - 1 + (int)((point - inp->dataMZ[config.lengthmz - 1]) / (inp->dataMZ[config.lengthmz - 1] - inp->dataMZ[config.lengthmz - 2]));
            end = config.lengthmz - 1 + nearfast(inp->dataMZ, 2 * inp->dataMZ[0] - point, config.lengthmz);
        } else {
            end = nearfast(inp->dataMZ, point, config.lengthmz);
        }
        endtab[i] = end;
        if (end - start > maxlength) { maxlength = end - start; }
        //printf("%d %d\n", start, end);
    }
    //printf("Max Length: %d\t", maxlength);
    return maxlength;
}



void CalcMasses(Config *config, Input *inp) {
    //Fills mtab by multiplying each mz by each z
    int ln = config->lengthmz * config->numz;
    inp->mtab = (float *) calloc(ln, sizeof(float));
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < config->lengthmz; i++) {
        for (int j = 0; j < config->numz; j++) {
            inp->mtab[index2D(config->numz, i, j)] = (
                inp->dataMZ[i] * inp->nztab[j] - config->adductmass * inp->nztab[j]);
        }
    }
}


int SetUpPeakShape(Config config, Input inp, Decon *decon, const int silent, const int verbose) {
    //Create a list of start and end values to box in arrays based on the above threshold
    decon->starttab = (int *) calloc(config.lengthmz, sizeof(int));
    decon->endtab = (int *) calloc(config.lengthmz, sizeof(int));
    int maxlength = 1;
    if (config.mzsig != 0) {
        //Gets maxlength and sets start and endtab
        maxlength = SetStartsEnds(config, &inp, decon->starttab, decon->endtab);

        //Changes dimensions of the peak shape function. 1D for speedy and 2D otherwise
        int pslen = config.lengthmz;
        if (config.speedyflag == 0) { pslen = config.lengthmz * maxlength; }
        if (verbose == 1) {
            printf("Maxlength: %d \t Total Size: %zu GB\n", maxlength, pslen * sizeof(float) / 1000000000);
        }
        decon->mzdist = (float *) calloc(pslen, sizeof(float));
        if (pslen * sizeof(float) / 1000000000 > 4) {
            printf("Danger: Your data may crash the memory. Consider setting the Peak FWHM to 0.\n");
        }
        int makereverse = 0;
        if (config.mzsig < 0 || config.beta < 0) {
            makereverse = 1;
            decon->rmzdist = (float *) calloc(pslen, sizeof(float));
        } else { decon->rmzdist = (float *) calloc(0, sizeof(float)); }
        printf("Test: %d\n", config.psfun);
        //Calculates the distance between mz values as a 2D or 3D matrix

        if (config.speedyflag == 0) {
            if (verbose == 1) { printf("Making Peak Shape 2D\n"); }
            MakePeakShape2D(config.lengthmz, maxlength, decon->starttab, decon->endtab, inp.dataMZ,
                            fabs(config.mzsig) * config.peakshapeinflate, config.psfun, config.speedyflag,
                            decon->mzdist, decon->rmzdist, makereverse);
        } else {
            if (verbose == 1) { printf("Making Peak Shape 1D\n"); }
            //Calculates peak shape as a 1D list centered at the first element for circular convolutions
            MakePeakShape1D(inp.dataMZ, config.psmzthresh, config.lengthmz, config.speedyflag,
                            fabs(config.mzsig) * config.peakshapeinflate,
                            config.psfun, decon->mzdist, decon->rmzdist, makereverse);
        }
        if (silent == 0) { printf("mzdist set: %f\t maxlength: %d\n", decon->mzdist[0], maxlength); }
    } else {
        decon->mzdist = (float *) calloc(0, sizeof(float));
        maxlength = 0;
    }
    return maxlength;
}



#endif
