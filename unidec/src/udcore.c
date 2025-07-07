//
// Created by mm96978 on 7/7/2025.
//

#include "udcore.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "udtools.h"


//..........
//
// Blur Functions
//
//............

//Convolution of neighborhood function with gaussian filter.
void blur_it(const int lengthmz,
             const int numz,
             const int numclose,
             const int *__restrict closeind,
             const float *__restrict closearray,
             float *__restrict newblur,
             const float *__restrict blur,
             const char *__restrict barr) {
    if (numclose == 1) {
        const size_t len = (size_t) lengthmz * numz * sizeof(float);
        memcpy(newblur, blur, len);
    } else {
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz * numz; i++) {
            float temp = 0;
            if (barr[i] == 1) {
                for (int k = 0; k < numclose; k++) {
                    if (closeind[index2D(numclose, i, k)] != -1) {
                        temp += closearray[index2D(numclose, i, k)] * blur[closeind[index2D(numclose, i, k)]];
                    }
                }
            }
            newblur[i] = temp;
        }
    }
}

//Charge state smooth using a mean filter of the log
void blur_it_mean(const int lengthmz,
                  const int numz,
                  const int numclose,
                  const int *__restrict closeind,
                  float *__restrict newblur,
                  const float *__restrict blur,
                  const char *__restrict barr,
                  const float *__restrict closearray,
                  const float zerolog) {
    if (numclose == 1) {
        const size_t len = (size_t) lengthmz * numz * sizeof(float);
        memcpy(newblur, blur, len);
    } else {
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz * numz; i++) {
            float temp = 0;
            if (barr[i] == 1) {
                for (int k = 0; k < numclose; k++) {
                    float temp2 = 0;
                    if (closeind[index2D(numclose, i, k)] != -1) {
                        temp2 = blur[closeind[index2D(numclose, i, k)]] * closearray[index2D(numclose, i, k)];
                    }
                    if (temp2 > 0) { temp += logf(temp2); } else { temp += zerolog; }
                }
                temp = expf(temp / (float) numclose);
            }
            newblur[i] = temp;
        }
    }
}


//Convolution of neighborhood function with gaussian filter.
void blur_it_hybrid1(const int lengthmz,
                     const int numz,
                     const int zlength,
                     const int mlength,
                     const int *__restrict closeind,
                     const int *__restrict closemind,
                     const int *__restrict closezind,
                     const float *__restrict mdist,
                     const float *__restrict zdist,
                     float *__restrict newblur,
                     const float *__restrict blur,
                     const char *__restrict barr,
                     const float *__restrict closearray,
                     const float zerolog) {

    const int numclose = zlength * mlength;
    if (numclose == 1) {
        const size_t len = (size_t) lengthmz * numz * sizeof(float);
        memcpy(newblur, blur, len);
    } else {
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            for (int j = 0; j < numz; j++) {
                float temp = 0;
                if (barr[index2D(numz, i, j)] == 1) {
                    for (int n = 0; n < mlength; n++) {
                        float temp2 = 0;
                        for (int k = 0; k < zlength; k++) {
                            const int m = index2D(mlength, k, n);
                            float temp3 = 0;
                            if (closeind[index3D(numz, numclose, i, j, m)] != -1) {
                                temp3 = blur[closeind[index3D(numz, numclose, i, j, m)]] * closearray[index3D(
                                            numz, numclose, i, j, m)];
                            }
                            if (temp3 > 0) { temp2 += logf(temp3); } else { temp2 += zerolog; }
                        }
                        temp += expf(temp2 / (float) zlength) * mdist[n];
                    }
                }
                newblur[index2D(numz, i, j)] = temp;
            }
        }
    }
}


//Convolution of neighborhood function with gaussian filter.
void blur_it_hybrid2(const int lengthmz,
                     const int numz,
                     const int zlength,
                     const int mlength,
                     const int *__restrict closeind,
                     const int *__restrict closemind,
                     const int *__restrict closezind,
                     const float *__restrict mdist,
                     const float *__restrict zdist,
                     float *__restrict newblur,
                     const float *__restrict blur,
                     const char *__restrict barr,
                     const float *__restrict closearray,
                     const float zerolog) {

    const int numclose = zlength * mlength;
    if (numclose == 1) {
        const size_t len = (size_t) lengthmz * numz * sizeof(float);
        memcpy(newblur, blur, len);
    } else {
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            for (int j = 0; j < numz; j++) {
                float temp = 0;
                if (barr[index2D(numz, i, j)] == 1) {
                    for (int n = 0; n < mlength; n++) {
                        float temp2 = 0;
                        for (int k = 0; k < zlength; k++) {
                            const int m = index2D(mlength, k, n);
                            if (closeind[index3D(numz, numclose, i, j, m)] != -1) {
                                temp2 += blur[closeind[index3D(numz, numclose, i, j, m)]] * zdist[k] * closearray[
                                    index3D(numz, numclose, i, j, m)];
                            }
                        }
                        if (temp2 > 0) { temp += logf(temp2); } // / (float)mlength);}
                        else { temp += zerolog; }
                    }
                    temp = expf(temp / (float) mlength);
                }
                newblur[index2D(numz, i, j)] = temp;
            }
        }
    }
}



void blur_baseline(float *baseline, const int lengthmz, const float *dataMZ, const float mzsig, int mult,
                   const int filterwidth) {
    const int mulin = mult;
    float *temp = NULL;
    temp = (float *) calloc(lengthmz, sizeof(float));
    memcpy(temp, baseline, sizeof(float) * lengthmz);
    if (temp == NULL || baseline == NULL || dataMZ == NULL) {
        fprintf(stderr, "Error allocating memory for temp in blur_baseline.\n");
        exit(11);
    }

    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++) {
        float mzdiff = 0;
        if (i > 0) {
            mzdiff = dataMZ[i] - dataMZ[i - 1];
        } else {
            mzdiff = dataMZ[i + 1] - dataMZ[i];
        }
        if (mulin == 0 && mzdiff > 0) {
            //mult = lengthmz / 500;
            mult = (int) (2 * mzsig / mzdiff);
        }
        if (mult < 1) { mult = 1; }

        float val = 0;
        const int window = filterwidth;
        for (int j = -window; j < window; j++) {
            const int k = i + j * (mult);
            float newval = 0;
            if (k >= 0 && k < lengthmz) {
                newval = temp[k];
            }
            if (k < 0) {
                newval = temp[-k];
            }
            if (k >= lengthmz) {
                newval = temp[2 * lengthmz - k];
            }
            //if(newval>0){
            val += newval;
            //}
        }
        baseline[i] = val / ((float) window * 2 + 1);
    }
    free(temp);
}


void midblur_baseline(float *baseline, const int lengthmz, const float *dataMZ, const float mzsig, int mult) {
    if (mult == 0) {
        mult = lengthmz / 400;
    }
    float *temp = NULL;
    temp = (float *) calloc(lengthmz, sizeof(float));
    if (temp) {
        memcpy(temp, baseline, sizeof(float) * lengthmz);
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            const int window = 25;
            const int len = window * 2;
            float *med = NULL;
            med = (float *) calloc(len, sizeof(float));
            if (med) {
                int index = 0;
                for (int j = -window; j < window; j++) {
                    int k = i + j * mult;
                    float newval = 0;
                    if (k >= 0 && k < lengthmz) {
                        newval = temp[k];
                    }
                    if (k < 0) {
                        newval = temp[-k];
                    }
                    if (k >= lengthmz) {
                        newval = temp[2 * lengthmz - k];
                    }
                    med[index] = newval;
                    index++;
                }
                qsort(med, len, sizeof(float), compare_function);
                index = 0;
                float val = 0;
                for (int j = 0; j < window; j++) {
                    val += med[j];
                    index++;
                }
                if (index != 0) {
                    baseline[i] = val / ((float) index);
                } else { baseline[i] = 0; }
                free(med);
            }
        }
        free(temp);
    }
}


//...............
//
// Convolution Functions and Richardson Lucy Functions
//
//...............


void convolve_simp(const int lengthmz, const int maxlength, const int *starttab, const int *endtab, const float *mzdist,
                   const float *deltas, float *denom, const int speedyflag) {
    if (speedyflag == 0) {
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            float cv = 0;
            for (int k = starttab[i]; k <= endtab[i]; k++) {
                const int k2 = fixk(k, lengthmz);
                const int start = starttab[k2];
                cv += deltas[k2] * mzdist[index2D(maxlength, k2, i - start)];
            }
            denom[i] = cv;
        }
    } else {
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            float cv = 0;
            for (int k = starttab[i]; k <= endtab[i]; k++) {
                cv += deltas[k] * mzdist[indexmod(lengthmz, k, i)];
            }
            denom[i] = cv;
        }
    }
}


void sum_deltas(const int lengthmz, const int numz, const float *__restrict blur, const char *__restrict barr,
                float *deltas) {

        //Collapse the grid into a 1D array of delta function values
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            float temp = 0;
            for (int j = 0; j < numz; j++) {
                if (barr[index2D(numz, i, j)] == 1) {
                    temp += blur[index2D(numz, i, j)];
                }
            }
            deltas[i] = temp;
        }
}

void apply_ratios(const int lengthmz, const int numz, const float *__restrict blur, const char *__restrict barr,
                  const float *__restrict denom, float *blur2) {

    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++) {
        for (int j = 0; j < numz; j++) {
            if (barr[index2D(numz, i, j)] == 1) {
                blur2[index2D(numz, i, j)] = denom[i] * blur[index2D(numz, i, j)];
            } else {
                blur2[index2D(numz, i, j)] = 0;
            }
        }
    }
}

//...............
//
// Larger Deconvolution Functions
//
//...............


void deconvolve_baseline(const int lengthmz, const float *dataMZ, const float *dataInt, float *baseline,
                         const float mzsig) {
    float *denom = NULL;
    denom = (float *) calloc(lengthmz, sizeof(float));
    if (denom) {
        midblur_baseline(baseline, lengthmz, dataMZ, mzsig, 0);
        midblur_baseline(baseline, lengthmz, dataMZ, mzsig, 5);

        memcpy(denom, baseline, sizeof(float) * lengthmz);
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            if (denom[i] != 0 && dataInt[i] >= 0) { denom[i] = dataInt[i] / denom[i]; }
        }

        midblur_baseline(denom, lengthmz, dataMZ, mzsig, 0);
        midblur_baseline(denom, lengthmz, dataMZ, mzsig, 5);
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            baseline[i] = baseline[i] * (denom[i]);
        }
        free(denom);
    }
}


float deconvolve_iteration_speedy(const int lengthmz, const int numz, const int maxlength, const float *__restrict blur,
                                  float *__restrict blur2,
                                  const char *__restrict barr, const int aggressiveflag,
                                  const float *__restrict dataInt,
                                  const int *starttab, const int *endtab, const float *mzdist, const float *rmzdist,
                                  const int speedyflag, const int baselineflag, float *baseline,
                                  float *noise, const float mzsig, const float *dataMZ, const int filterwidth,
                                  const float psig) {

    float *deltas = NULL, *denom = NULL;
    deltas = (float *) calloc(lengthmz, sizeof(float));
    denom = (float *) calloc(lengthmz, sizeof(float));
    if (deltas == NULL || denom == NULL) {
        fprintf(stderr, "Error allocating memory for deltas or denom in deconvolve_iteration_speedy.\n");
        exit(11);
    }

    if (aggressiveflag == 1 && mzsig != 0) {
        blur_baseline(baseline, lengthmz, dataMZ, fabsf(mzsig), 0, filterwidth);
    }
    //printf("1\n");
    //Sum deltas
    sum_deltas(lengthmz, numz, blur, barr, deltas);
    //printf("2\n");
    if (mzsig != 0 && psig >= 0) {
        //Convolve with peak shape
        convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, denom, speedyflag);
    } else {
        memcpy(denom, deltas, sizeof(float) * lengthmz);
    }
    //printf("3\n");
    if (aggressiveflag == 1) {
#pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            denom[i] += baseline[i]; // +noise[i]);
        }
    }
    //printf("4\n");
    //Calculate Ratio
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++) {
        if (denom[i] != 0 && dataInt[i] >= 0) { denom[i] = dataInt[i] / denom[i]; }
    }
    //printf("5\n");
    if (mzsig < 0) {
        //Real Richardson-Lucy Second Convolution
        convolve_simp(lengthmz, maxlength, starttab, endtab, rmzdist, denom, deltas, speedyflag);
        memcpy(denom, deltas, sizeof(float) * lengthmz);
        /*
        for (i = 0; i<lengthmz; i++)
        {
            deltas[lengthmz - 1 - i] = denom[i];
        }
        convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, denom, speedyflag);
        for (i = 0; i<lengthmz; i++)
        {
            deltas[lengthmz - 1 - i] = denom[i];
        }
        memcpy(denom, deltas, sizeof(float)*lengthmz);*/
    }

    //printf("6\n");
    //Multiply Ratio by prior
    apply_ratios(lengthmz, numz, blur, barr, denom, blur2);

    //printf("7\n");
    if (aggressiveflag == 1) {
        //memcpy(deltas, denom, sizeof(float)*lengthmz);
        blur_baseline(denom, lengthmz, dataMZ, fabsf(mzsig), 0, filterwidth);
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            baseline[i] = baseline[i] * (denom[i]);
            //noise[i] = noise[i]*(deltas[i]);
        }
    }
    //printf("8\n");
    free(deltas);
    free(denom);
    return 0;
}


float getfitdatspeedy(float *fitdat, const float *blur, const int lengthmz, const int numz,
                      const int maxlength, const float maxint,
                      const int *starttab,
                      const int *endtab, const float *mzdist, const int speedyflag) {

    float *deltas = NULL;
    deltas = (float *) calloc(lengthmz, sizeof(float));

    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < lengthmz; i++) //Collapse the grid into a 1D array of delta function values
    {
        float temp = 0;
        for (int j = 0; j < numz; j++) {
            temp += blur[index2D(numz, i, j)];
        }
        deltas[i] = temp;
    }

    if (maxlength != 0) {
        convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, fitdat, speedyflag);
    } else {
        memcpy(fitdat, deltas, sizeof(float) * lengthmz);
    }

    free(deltas);
    float fitmax = 0;
    //#pragma omp parallel for private (i), schedule(auto)
    for (int i = 0; i < lengthmz; i++) {
        if (fitdat[i] > fitmax) { fitmax = fitdat[i]; }
    }
    //#pragma omp parallel for private (i), schedule(auto)
    if (fitmax != 0) {
        for (int i = 0; i < lengthmz; i++) {
            if (fitdat[i] < 0) { fitdat[i] = 0; } else { fitdat[i] = fitdat[i] * maxint / fitmax; }
        }
    }
    return fitmax;
}

float errfunspeedy(Config config, Decon decon, const char *barr, const float *dataInt, const int maxlength,
                   const int *starttab, const int *endtab,
                   const float *mzdist, float *rsquared) {
    //Get max intensity
    float maxint = 0;
    for (int i = 0; i < config.lengthmz; i++) {
        if (dataInt[i] > maxint) { maxint = dataInt[i]; }
    }

    getfitdatspeedy(decon.fitdat, decon.blur, config.lengthmz, config.numz, maxlength,
                    maxint, starttab, endtab, mzdist, config.speedyflag);

    if (config.baselineflag == 1) {
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < config.lengthmz; i++) {
            decon.fitdat[i] += decon.baseline[i]; // +decon.noise[i];
            //decon.fitdat[i] = decon.noise[i]+0.1;
        }
    }
    ApplyCutoff(decon.fitdat, 0, config.lengthmz);

    float fitmean = Average(config.lengthmz, dataInt);

    float error = 0;
    float sstot = 0;
    for (int i = 0; i < config.lengthmz; i++) {
        error += pow((decon.fitdat[i] - dataInt[i]), 2);
        sstot += pow((dataInt[i] - fitmean), 2);
    }

    //Calculate R-squared
    if (sstot != 0) { *rsquared = 1 - (error / sstot); }

    return error;
}
