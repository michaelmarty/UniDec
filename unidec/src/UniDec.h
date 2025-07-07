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


void softargmax(float *blur, const int lengthmz, const int numz, const float beta) {
    // if (beta < 0) {
    //     softargmax_everything(blur, lengthmz, numz, fabs(beta));
    //     return;
    // }

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




#endif
