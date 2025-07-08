//
// Created by mm96978 on 7/7/2025.
//

#include "udcore.h"


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


void softargmax(float *blur, const int lengthmz, const int numz, const float beta) {
    const int l = lengthmz * numz;
    const size_t sizel = (size_t) l * sizeof(float);
    float *newblur = (float *) calloc(l, sizeof(float));
    if (newblur) {
        memcpy(newblur, blur, sizel);
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < lengthmz; i++) {
            float sum2 = 0;
            float sum1 = 0;
            float factor = 0;
            float min2 = 1000000000000.0f;

            for (int j = 0; j < numz; j++) {
                const float d = newblur[index2D(numz, i, j)];
                sum1 += d;

                const float e = expf(beta * d);

                if (e < min2) { min2 = e; }

                blur[index2D(numz, i, j)] = e;
                sum2 += e;
            }

            const float denom = (sum2 - min2 * (float) numz);
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
}




void point_smoothing(float *blur, const char *barr, const int lengthmz, const int numz, const int width) {
    const int l = lengthmz * numz;
    const float fwidth = (float) width;
    float *newblur = calloc(l, sizeof(float));
    if (newblur) {
        const size_t newl = (size_t) lengthmz * numz * sizeof(float);
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
                    blur[index2D(numz, i, j)] = sum / (1.0f + 2.0f * fwidth);
                }
            }
        }
        free(newblur);
    }
}



float getfitdatspeedy(float *fitdat, const float *blur, const int lengthmz, const int numz,
                      const int maxlength, const float maxint,
                      const int *starttab,
                      const int *endtab, const float *mzdist, const int speedyflag) {

    float *deltas = NULL;
    deltas = (float *) calloc(lengthmz, sizeof(float));
    if (deltas == NULL) {
        fprintf(stderr, "Error allocating memory for deltas in getfitdatspeedy.\n");
        exit(11);
    }

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

    const float fitmean = Average(config.lengthmz, dataInt);

    float error = 0;
    float sstot = 0;
    for (int i = 0; i < config.lengthmz; i++) {
        error += powf((decon.fitdat[i] - dataInt[i]), 2.f);
        sstot += powf((dataInt[i] - fitmean), 2.f);
    }

    //Calculate R-squared
    if (sstot != 0) { *rsquared = 1 - (error / sstot); }

    return error;
}

// ...................
//
// Test Mass Functions
//
// ....................


void CalcMasses(const Config *config, Input *inp) {
    //Fills mtab by multiplying each mz by each z
    const int ln = config->lengthmz * config->numz;
    inp->mtab = (float *) calloc(ln, sizeof(float));
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < config->lengthmz; i++) {
        for (int j = 0; j < config->numz; j++) {
            inp->mtab[index2D(config->numz, i, j)] = calcmass(inp->dataMZ[i], inp->nztab[j], config->adductmass);
        }
    }
}

void TestMassListWindowed(const Config config, Input *inp) {
#pragma omp parallel for schedule(auto)
    for (int i = 0; i < config.lengthmz; i++) {
        for (int j = 0; j < config.numz; j++) {
            inp->barr[index2D(config.numz, i, j)] = 0;
            const float testmass = inp->mtab[index2D(config.numz, i, j)];
            const float nativelimit = nativecharge(testmass, 0);
            if (testmass < config.massub && testmass > config.masslb && (float) inp->nztab[j] < nativelimit + config.nativezub && (float) inp->nztab[j] > nativelimit +
                config.nativezlb) {
                if (neartest(inp->testmasses, testmass, config.mfilelen, config.mtabsig) == 1) {
                    inp->barr[index2D(config.numz, i, j)] = 1;
                }
            }
        }
    }
}

void TestMassListLimit(const Config config, Input *inp, const int *testmasspos) {

    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < config.lengthmz; i++) {
        for (int j = 0; j < config.numz; j++) {
            inp->barr[index2D(config.numz, i, j)] = 0;
            const float testmass = inp->mtab[index2D(config.numz, i, j)];
            const float nativelimit = nativecharge(testmass, 0);
            if (testmass < config.massub && testmass > config.masslb && (float) inp->nztab[j] < nativelimit + config.nativezub && (float) inp->nztab[j] > nativelimit +
                config.nativezlb) {
                for (int k = 0; k < config.mfilelen; k++) {
                    if (testmasspos[index2D(config.numz, k, j)] == i) {
                        inp->barr[index2D(config.numz, i, j)] = 1;
                    }
                }
            }
        }
    }
}

void TestMass(const Config config, Input *inp) {
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < config.lengthmz; i++) {
        for (int j = 0; j < config.numz; j++) {
            const float testmass = inp->mtab[index2D(config.numz, i, j)];
            const float nativelimit = nativecharge(testmass, 0);
            if (testmass < config.massub && testmass > config.masslb && (float) inp->nztab[j] < nativelimit + config.nativezub && (float) inp->nztab[j] > nativelimit +
                config.nativezlb) {
                inp->barr[index2D(config.numz, i, j)] = 1;
            } else {
                inp->barr[index2D(config.numz, i, j)] = 0;
            }
        }
    }
}


void SetLimits(const Config config, Input *inp) {
    //If there is a mass file read, it will only allow masses close to those masses within some config.mtabsig window.
    if (config.mflag == 1 && config.limitflag == 0) {
        TestMassListWindowed(config, inp);
    }
    //If there is a mass file read and the mass table window (config.mtabsig) is 0, it will only write intensities at the m/z values closest to the m/z values read in from the mfile.
    else if (config.mflag == 1 && config.limitflag == 1) {
        //Determines the indexes of each test mass from mfile in m/z space
        int *testmasspos = malloc(sizeof(int) * config.mfilelen * config.numz);
        if (testmasspos) {
            for (int i = 0; i < config.mfilelen; i++) {
                for (int j = 0; j < config.numz; j++) {
                    const float mztest = calcmz(inp->testmasses[i], inp->nztab[j], config.adductmass);
                       // (inp->testmasses[i] + config.adductmass * inp->nztab[j]) / (float) inp->nztab[j];
                    testmasspos[index2D(config.numz, i, j)] = nearfast(inp->dataMZ, mztest, config.lengthmz);
                }
            }
            TestMassListLimit(config, inp, testmasspos);
        }
        free(testmasspos);
    }
    //Normally, write the intensity values if the values fall within the mass upperbound and lower bound
    else {
        TestMass(config, inp);
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


void KillB(const float *I, char *B, const float intthresh, const int lengthmz, const int numz) {
    for (int i = 0; i < lengthmz; i++) {
        for (int j = 0; j < numz; j++) {
            if (I[i] <= intthresh) { B[index2D(numz, i, j)] = 0; }
        }
    }
}

// .......................
//
// Peak Shape Functions
//
// .......................


void MakeSparseBlur(const int numclose, char *barr, const int *closezind,
                    const int *closemind, int *closeind, const float *closeval, float *closearray, const Config config, const Input *inp) {
    int lengthmz = config.lengthmz;
    int numz = config.numz;
    float molig = config.molig;

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
                    mzsig = 2.0f * fabsf(inp->dataMZ[i2] - inp->dataMZ[i1]);
                    if (mzsig > config.massbins || mzsig == 0) { mzsig = config.massbins * 2; }
                }
                float newthreshold = mzsig * 2;

                for (int k = 0; k < numclose; k++) {
                    //Find the z index and test if it is within the charge range
                    int indz = (int) (j + closezind[k]);
                    if (indz < 0 || indz >= numz || (inp->nztab[j] + closezind[k]) == 0) {
                        closeind[index3D(numz, numclose, i, j, k)] = -1;
                        closearray[index3D(numz, numclose, i, j, k)] = 0;
                    } else {
                        //Find the nearest m/z value and test if it is close enough to the predicted one and within appropriate ranges
                        float point = calcmz(inp->mtab[index2D(numz, i, j)] + (float) closemind[k] * molig, inp->nztab[j] + closezind[k], config.adductmass);
                        if (point < inp->dataMZ[0] - newthreshold || point > inp->dataMZ[lengthmz - 1] + newthreshold) {
                            closeind[index3D(numz, numclose, i, j, k)] = -1;
                            closearray[index3D(numz, numclose, i, j, k)] = 0;
                        } else {
                            int ind = nearfast(inp->dataMZ, point, lengthmz);
                            float closepoint = inp->dataMZ[ind];
                            int newind = index2D(numz, ind, indz);
                            if (barr[newind] == 1 && fabsf(point - closepoint) < newthreshold) {
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


void MakePeakShape2D(const Config config, Decon *decon, int maxlength, const float *dataMZ, const int makereverse, const int inflateflag) {


    float mzsig = fabsf(config.mzsig);
    if (inflateflag == 1) {
        mzsig *= config.peakshapeinflate;
    }
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < config.lengthmz; i++) {
        int start = decon->starttab[i];
        int end = decon->endtab[i];
        for (int j = start; j <= end; j++) {
            int j2 = fixk(j, config.lengthmz);
            decon->mzdist[index2D(maxlength, i, j2 - start)] = mzpeakshape(dataMZ[i], dataMZ[j2], mzsig, config.psfun);
            if (makereverse == 1) {
                decon->rmzdist[index2D(maxlength, i, j2 - start)] = mzpeakshape(dataMZ[j2], dataMZ[i], mzsig, config.psfun);
            }
        }
    }
}

void MakePeakShape1D(const Config config, Decon * decon, const float *dataMZ, const int makereverse, const int inflateflag) {
    float binsize = dataMZ[1] - dataMZ[0];
    float newrange = config.psmzthresh / binsize;

    float mzsig = fabsf(config.mzsig);
    if (inflateflag == 1) {
        mzsig *= config.peakshapeinflate;
    }

    int n;
    for (n = (int) -newrange; n < (int) newrange; n++) {
        decon->mzdist[indexmod(config.lengthmz, 0, n)] = mzpeakshape(0, (float) n * binsize, mzsig, config.psfun);
        if (makereverse == 1) { decon->rmzdist[indexmod(config.lengthmz, 0, n)] = mzpeakshape((float) n * binsize, 0, mzsig, config.psfun); }
    }
    printf("\nNotice: Assuming linearized data. \n\n");
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

int SetUpPeakShape(Config config, Input inp, Decon *decon, const int silent, const int verbose) {
    //Create a list of start and end values to box in arrays based on the above threshold
    decon->starttab = (int *) calloc(config.lengthmz, sizeof(int));
    decon->endtab = (int *) calloc(config.lengthmz, sizeof(int));
    if (decon->starttab == NULL || decon->endtab == NULL) {
        fprintf(stderr, "Error allocating memory for starttab or endtab in SetUpPeakShape.\n");
        exit(11);
    }

    int maxlength;
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
        if (decon->mzdist == NULL) {
            fprintf(stderr, "Error allocating memory for mzdist in SetUpPeakShape.\n");
            exit(11);
        }
        if (pslen * sizeof(float) / 1000000000 > 4) {
            printf("Danger: Your data may crash the memory. Consider setting the Peak FWHM to 0.\n");
        }
        int makereverse = 0;
        if (config.mzsig < 0 || config.beta < 0) {
            makereverse = 1;
            decon->rmzdist = (float *) calloc(pslen, sizeof(float));
        } else { decon->rmzdist = (float *) calloc(0, sizeof(float)); }
        if (decon->rmzdist == NULL) {
            fprintf(stderr, "Error allocating memory for rmzdist in SetUpPeakShape.\n");
            exit(11);
        }
        printf("Test: %d\n", config.psfun);
        //Calculates the distance between mz values as a 2D or 3D matrix

        if (config.speedyflag == 0) {
            if (verbose == 1) { printf("Making Peak Shape 2D\n"); }
            MakePeakShape2D(config, decon, maxlength, inp.dataMZ, makereverse, 1);
        } else {
            if (verbose == 1) { printf("Making Peak Shape 1D\n"); }
            //Calculates peak shape as a 1D list centered at the first element for circular convolutions
            MakePeakShape1D(config, decon, inp.dataMZ, makereverse, 0);
        }
        if (silent == 0) { printf("mzdist set: %f\t maxlength: %d\n", decon->mzdist[0], maxlength); }
    } else {
        decon->mzdist = (float *) calloc(0, sizeof(float));
        maxlength = 0;
    }
    return maxlength;
}

float Reconvolve(const Config config, const int maxlength, Decon *decon, const char *barr) {
    float newblurmax = 0;
    if (config.speedyflag == 0) {
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < config.lengthmz; i++) {
            for (int j = 0; j < config.numz; j++) {
                float cv = 0;
                if (barr[index2D(config.numz, i, j)] == 1) {
                    for (int k = decon->starttab[i]; k <= decon->endtab[i]; k++) {
                        int k2 = fixk(k, config.lengthmz);
                        if (decon->blur[index2D(config.numz, k2, j)] != 0) {
                            int start = decon->starttab[k2];
                            cv += decon->blur[index2D(config.numz, k2, j)] * decon->mzdist[index2D(maxlength, k2, i - start)];
                        }
                    }
                }
                decon->newblur[index2D(config.numz, i, j)] = cv;
                if (cv > newblurmax) {
                    newblurmax = cv;
                }
            }
        }
    } else {
        #pragma omp parallel for schedule(auto)
        for (int i = 0; i < config.lengthmz; i++) {
            for (int j = 0; j < config.numz; j++) {
                float cv = 0;
                if (barr[index2D(config.numz, i, j)] == 1) {
                    for (int k = decon->starttab[i]; k <= decon->endtab[i]; k++) {
                        if (decon->blur[index2D(config.numz, k, j)] != 0) {
                            cv += decon->blur[index2D(config.numz, k, j)] * decon->mzdist[indexmod(config.lengthmz, k, i)];
                        }
                    }
                }
                decon->newblur[index2D(config.numz, i, j)] = cv;
                if (cv > newblurmax) {
                    newblurmax = cv;
                }
            }
        }
    }
    return newblurmax;
}


void charge_scaling(float *blur, const int *nztab, const int lengthmz, const int numz) {
    for (int i = 0; i < lengthmz; i++) {
        for (int j = 0; j < numz; j++) {
            const int charge = nztab[j];
            const float z = (float) charge;
            if (z != 0) { blur[index2D(numz, i, j)] /= z; }
        }
    }
}

// .......................
//
// Transform Functions
//
// .......................

void IntegrateTransform(const Config config, Decon *decon, const float *mtab, float massmax, float massmin) {

    const float *blur;
    if (config.rawflag == 1 || config.rawflag == 3) {
        blur = decon->blur;
    }
    else if (config.rawflag == 0 || config.rawflag == 2) {
        blur = decon->newblur;
    }
    else {
        return; // Do nothing
    }

    for (int i = 0; i < config.lengthmz; i++) {
        for (int j = 0; j < config.numz; j++) {
            float testmass = mtab[index2D(config.numz, i, j)];
            if (testmass < massmax && testmass > massmin) {
                int index = nearfast(decon->massaxis, testmass, decon->mlen);
                float newval = blur[index2D(config.numz, i, j)];
                if (decon->massaxis[index] == testmass) {
                    decon->massaxisval[index] += newval;
                    decon->massgrid[index2D(config.numz, index, j)] += newval;
                }

                if (decon->massaxis[index] < testmass && index < decon->mlen - 2) {
                    int index2 = index + 1;
                    float interpos = LinearInterpolatePosition(decon->massaxis[index], decon->massaxis[index2], testmass);
                    decon->massaxisval[index] += (1.0f - interpos) * newval;
                    decon->massgrid[index2D(config.numz, index, j)] += (1.0f - interpos) * newval;
                    decon->massaxisval[index2] += (interpos) * newval;
                    decon->massgrid[index2D(config.numz, index2, j)] += (interpos) * newval;
                }

                if (decon->massaxis[index] > testmass && index > 0) {
                    int index2 = index - 1;
                    float interpos = LinearInterpolatePosition(decon->massaxis[index], decon->massaxis[index2], testmass);
                    decon->massaxisval[index] += (1 - interpos) * newval;
                    decon->massgrid[index2D(config.numz, index, j)] += (1 - interpos) * newval;
                    decon->massaxisval[index2] += (interpos) * newval;
                    decon->massgrid[index2D(config.numz, index2, j)] += (interpos) * newval;
                }
            }
        }
    }
}

void InterpolateTransform(const Config config, Decon *decon, const Input *inp) {
    const float *blur;
    if (config.rawflag == 1 || config.rawflag == 3) {
        blur = decon->blur;
    }
    else if (config.rawflag == 0 || config.rawflag == 2) {
        blur = decon->newblur;
    }
    else {
        return; // Do nothing
    }

    float startmzval = inp->dataMZ[0];
    float endmzval = inp->dataMZ[config.lengthmz - 1];
    //#pragma omp parallel for schedule(auto)
    for (int i = 0; i < decon->mlen; i++) {
        float val = 0;

        for (int j = 0; j < config.numz; j++) {

            float mztest = calcmz(decon->massaxis[i], inp->nztab[j], config.adductmass);

            if (mztest > startmzval && mztest < endmzval) {
                int index = nearfast(inp->dataMZ, mztest, config.lengthmz);
                int index2 = index;
                float newval = 0;
                if (inp->dataMZ[index] == mztest) {
                    newval = blur[index2D(config.numz, index, j)];
                    val += newval;
                    decon->massgrid[index2D(config.numz, i, j)] = newval;
                } else {
                    if (inp->dataMZ[index] > mztest && index > 1 && index < config.lengthmz - 1) {
                        index2 = index;
                        index = index - 1;
                    } else if (inp->dataMZ[index] < mztest && index < config.lengthmz - 2 && index > 0) {
                        index2 = index + 1;
                    }
                    if (index2 > index && (inp->dataMZ[index2] - inp->dataMZ[index]) != 0) {
                        float mu = (mztest - inp->dataMZ[index]) / (inp->dataMZ[index2] - inp->dataMZ[index]);
                        float y0 = blur[index2D(config.numz, index - 1, j)];
                        float y1 = blur[index2D(config.numz, index, j)];
                        float y2 = blur[index2D(config.numz, index2, j)];
                        float y3 = blur[index2D(config.numz, index2 + 1, j)];
                        newval = clip(CubicInterpolate(y0, y1, y2, y3, mu), 0);
                        //newval=CRSplineInterpolate(y0,y1,y2,y3,mu);
                        //newval=LinearInterpolate(y1,y2,mu);
                        val += newval;
                        decon->massgrid[index2D(config.numz, i, j)] = newval;
                    }
                }
            }
        }
        decon->massaxisval[i] = val;
    }
}

void SmartTransform(const Config config, Decon *decon, const Input *inp) {
    const float *blur;
    if (config.rawflag == 1 || config.rawflag == 3) {
        blur = decon->blur;
    }
    else if (config.rawflag == 0 || config.rawflag == 2) {
        blur = decon->newblur;
    }
    else {
        return; // Do nothing
    }

    float startmzval = inp->dataMZ[0];
    float endmzval = inp->dataMZ[config.lengthmz - 1];
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i < decon->mlen; i++) {
        float val = 0;
        for (int j = 0; j < config.numz; j++) {
            int z = inp->nztab[j];
            float mtest = decon->massaxis[i];
            float mztest = calcmz(mtest, z, config.adductmass);

            float mzlower;
            float mlower;
            if (i > 0) {
                mlower = decon->massaxis[i - 1];
                mzlower = calcmz(mlower, z, config.adductmass);
            } else {
                mzlower = mztest;
                mlower = mtest;
            }

            float mzupper;
            float mupper;
            if (i < decon->mlen - 1) {
                mupper = decon->massaxis[i + 1];
                mzupper = calcmz(mupper, z, config.adductmass);
            } else {
                mzupper = mztest;
                mupper = mtest;
            }


            if (mzupper > startmzval && mzlower < endmzval) {
                int index = nearfast(inp->dataMZ, mztest, config.lengthmz);
                int index1 = nearfast(inp->dataMZ, mzlower, config.lengthmz);
                int index2 = nearfast(inp->dataMZ, mzupper, config.lengthmz);
                float imz = inp->dataMZ[index];
                float newval = 0;
                if (index2 - index1 < 5) {
                    if (imz == mztest) {
                        newval = clip(blur[index2D(config.numz, index, j)], 0);
                        val += newval;
                        decon->massgrid[index2D(config.numz, i, j)] = newval;
                    } else {
                        int edge = 0;
                        index2 = index;
                        if (imz > mztest) {
                            index = index - 1;
                        } else if (imz < mztest) {
                            index2 = index + 1;
                        }

                        if (index < 1 || index2 >= config.lengthmz - 1) {
                            edge = 1;
                        }
                        if (index < 0 || index2 >= config.lengthmz) {
                            edge = 2;
                        }


                        if (index2 > index && (inp->dataMZ[index2] - inp->dataMZ[index]) != 0 && edge == 0) {
                            float mu = (mztest - inp->dataMZ[index]) / (inp->dataMZ[index2] - inp->dataMZ[index]);
                            float y0 = blur[index2D(config.numz, index - 1, j)];
                            float y1 = blur[index2D(config.numz, index, j)];
                            float y2 = blur[index2D(config.numz, index2, j)];
                            float y3 = blur[index2D(config.numz, index2 + 1, j)];
                            newval = clip(CubicInterpolate(y0, y1, y2, y3, mu), 0);
                            //newval=clip(CRSplineInterpolate(y0,y1,y2,y3,mu), 0);
                            //newval=clip(LinearInterpolate(y1,y2,mu),0);
                            val += newval;
                            decon->massgrid[index2D(config.numz, i, j)] = newval;
                            //printf("0\n");
                        } else if (edge == 1 && (inp->dataMZ[index2] - inp->dataMZ[index]) != 0) {
                            float mu = (mztest - inp->dataMZ[index]) / (inp->dataMZ[index2] - inp->dataMZ[index]);
                            float y1 = blur[index2D(config.numz, index, j)];
                            float y2 = blur[index2D(config.numz, index2, j)];
                            newval = clip(LinearInterpolate(y1, y2, mu), 0);
                            val += newval;
                            decon->massgrid[index2D(config.numz, i, j)] = newval;
                            //printf("1 %d %d %f %f %f\n", index, index2, mztest, imz, newval, massaxis[i]);
                        } else if (edge == 2 && (inp->dataMZ[index2] - inp->dataMZ[index]) != 0) {
                            if (index2 == 0) {
                                index = 0;
                                index2 = 1;
                            }
                            if (index == config.lengthmz - 1) {
                                index = config.lengthmz - 1;
                                index2 = config.lengthmz - 2;
                            }
                            float mu = (mztest - inp->dataMZ[index]) / (inp->dataMZ[index] - inp->dataMZ[index2]);
                            float y1 = blur[index2D(config.numz, index, j)];
                            float y2 = 0;
                            newval = clip(LinearInterpolate(y1, y2, mu), 0);
                            val += newval;
                            decon->massgrid[index2D(config.numz, i, j)] = newval;
                            //printf("2\n");
                        }
                    }
                } else {
                    //printf("Integrating\n");
                    float num = 0;
                    for (int k = index1; k < index2 + 1; k++) {
                        float kmz = inp->dataMZ[k];
                        float km = (kmz - config.adductmass) * (float) z;
                        float scale;

                        if (mztest < kmz && km < mupper) {
                            //scale = LinearInterpolatePosition(mzupper, mztest, kmz);
                            scale = LinearInterpolatePosition(mupper, mtest, km);
                        } else if (kmz < mztest && km > mlower) {
                            //scale = LinearInterpolatePosition(mzlower, mztest, kmz);
                            scale = LinearInterpolatePosition(mlower, mtest, km);
                        } else if (kmz == mztest) { scale = 1; } else { scale = 0; }

                        newval += scale * blur[index2D(config.numz, k, j)];
                        num += scale;
                    }
                    if (num != 0) { newval /= num; }
                    newval = clip(newval, 0);
                    val += newval;
                    decon->massgrid[index2D(config.numz, i, j)] = newval;
                }
            }
        }
        decon->massaxisval[i] = val;
    }
}



// .....................
//
// Double Dec
//
// .......................


// Integrates data in kernel to match sampling in data. Returns the new kernel length
int integrate_dd(const float *kernel_x, const float *kernel_y, const int kernellen, const float *data_x,
                 const int datalen, float **kernel_x_new, float **kernel_y_new) {
    if (kernellen <= 1 || datalen <= 1) {
        return kernellen;
    }
    const float diff = data_x[1] - data_x[0]; // kernel sampling needs to match this
    const float kdiff = kernel_x[1] - kernel_x[0]; // the original kernel sampling
    const int newlen = (int) ((kernel_x[kernellen - 1] - kernel_x[0]) / diff) + 1; // new number of points in kernel
    // these are the newly allocated arrays for the kernel
    int truelen;
    if (newlen > datalen) {
        truelen = newlen;
    } else {
        truelen = datalen;
    }
    *kernel_x_new = (float *) calloc(truelen, sizeof(float));
    *kernel_y_new = (float *) calloc(truelen, sizeof(float));
    if (*kernel_x_new == NULL || *kernel_y_new == NULL) {
        fprintf(stderr, "Error allocating memory for kernel_x_new or kernel_y_new in integrate_dd.\n");
        exit(11);
    }

    if (*kernel_x_new && *kernel_y_new) {
        float current_x = kernel_x[0];
        float current_xl = kernel_x[0];
        float current_xr = kernel_x[0] + (diff / 2);
        int current_index = 0;
        for (int i = 0; i < newlen; i++) {
            float y_val = 0;
            for (int j = current_index; j < kernellen; j++) {
                // For the first value, add area to the left of the point
                if (j == current_index && j != 0 && kernel_x[j] >= current_xl && kernel_x[j] < current_xr) {
                    const float left_mu = LinearInterpolatePosition(kernel_x[j - 1], kernel_x[j], current_xl);
                    const float left_y = LinearInterpolate(kernel_y[j - 1], kernel_y[j], left_mu);
                    y_val += (left_y + kernel_y[j]) * (kernel_x[j] - current_xl) / 2.0f;
                }
                // Next, add the area to the right of the point (it's either to the next point, to the
                // boundary, or we're at the last point)
                if (kernel_x[j] >= current_xl && kernel_x[j] < current_xr && (j + 1) < kernellen &&
                    kernel_x[j + 1] < current_xr) {
                    y_val += (kernel_y[j] + kernel_y[j + 1]) * kdiff / 2.0f;
                } else if (kernel_x[j] >= current_xl && kernel_x[j] < current_xr &&
                           (j + 1) < kernellen && kernel_x[j + 1] >= current_xr) {
                    const float right_mu = LinearInterpolatePosition(kernel_x[j], kernel_x[j + 1], current_xr);
                    const float right_y = LinearInterpolate(kernel_y[j], kernel_y[j + 1], right_mu);
                    y_val += (kernel_y[j] + right_y) * (current_xr - kernel_x[j]) / 2.0f;
                } else if (kernel_x[j] >= current_xr || (j + 1) >= kernellen) {
                    current_index = j;
                    break;
                }
            }
            (*kernel_x_new)[i] = current_x;
            (*kernel_y_new)[i] = y_val; // Should probably divide by diff. CHECK! --actually, I don't think it matters
            current_x += diff;
            current_xl = current_xr;
            current_xr += diff;
        }
        // free(kernel_x_new);
        // free(kernel_y_new);
    }
    return newlen;
}

// (Linear) Interpolates data in kernel to match sampling in data. Returns the new kernel length
int interpolate_dd(const float *kernel_x, const float *kernel_y, const int kernellen, const float *data_x,
                   int datalen, float **kernel_x_new, float **kernel_y_new) {
    if (kernellen <= 1 || datalen <= 1) {
        return kernellen;
    }
    const float diff = data_x[1] - data_x[0]; // kernel sampling needs to match this
    const int newlen = (int) ((kernel_x[kernellen - 1] - kernel_x[0]) / diff) + 1; // new number of points in kernel
    // these are the newly allocated arrays for the kernel
    int truelen;
    if (newlen > datalen) {
        truelen = newlen;
    } else {
        truelen = datalen;
    }
    *kernel_x_new = (float *) calloc(truelen, sizeof(float));
    *kernel_y_new = (float *) calloc(truelen, sizeof(float));
    if (*kernel_x_new && *kernel_y_new) {
        float current_x = kernel_x[0];
        for (int i = 0; i < newlen; i++) {
            const int nearest_index = nearfast(kernel_x, current_x, kernellen);
            if (kernel_x[nearest_index] == current_x) {
                (*kernel_y_new)[i] = kernel_y[nearest_index];
            } else if (kernel_x[nearest_index] < current_x) {
                if ((nearest_index + 1) < kernellen) {
                    const float mu = LinearInterpolatePosition(kernel_x[nearest_index], kernel_x[nearest_index + 1],
                                                           current_x);
                    const float y_val = LinearInterpolate(kernel_y[nearest_index], kernel_y[nearest_index + 1], mu);
                    (*kernel_y_new)[i] = y_val;
                } else {
                    // this should never be the case
                    (*kernel_y_new)[i] = kernel_y[nearest_index];
                }
            } else if (kernel_x[nearest_index] > current_x) {
                if (nearest_index > 0) {
                    const float mu = LinearInterpolatePosition(kernel_x[nearest_index - 1], kernel_x[nearest_index],
                                                           current_x);
                    const float y_val = LinearInterpolate(kernel_y[nearest_index - 1], kernel_y[nearest_index], mu);
                    (*kernel_y_new)[i] = y_val;
                } else {
                    // this should also never happen
                    (*kernel_y_new)[i] = kernel_y[nearest_index];
                }
            }
            (*kernel_x_new)[i] = current_x;
            current_x += diff;
        }
        // free(kernel_x_new);
        // free(kernel_y_new);
    }
    return newlen;
}

void cconv2fast(float *a, float *b, float *c, const int length) {
    // We don't seem to have complex.h, unfortunately
    // fftwf_complex is a double[2] of real (0) and imaginary (1)
    const int complen = (length / 2) + 1;
    fftwf_complex *A = fftwf_malloc(complen * sizeof(fftwf_complex));
    fftwf_complex *B = fftwf_malloc(complen * sizeof(fftwf_complex));
    fftwf_complex *C = fftwf_malloc(complen * sizeof(fftwf_complex));

    // Possible issue: create plan after we've created inputs (params)
    fftwf_plan p1 = fftwf_plan_dft_r2c_1d(length, a, A, FFTW_ESTIMATE);
    fftwf_plan p2 = fftwf_plan_dft_r2c_1d(length, b, B, FFTW_ESTIMATE);
    fftwf_execute(p1);
    fftwf_execute(p2);

    // A * B = (ac - bd) + i(ad + bc)
    for (int j = 0; j < complen; j++) {
        C[j][0] = (A[j][0] * B[j][0]) - (A[j][1] * B[j][1]);
        C[j][1] = (A[j][0] * B[j][1]) + (A[j][1] * B[j][0]);
    }

    fftwf_plan p3 = fftwf_plan_dft_c2r_1d(length, C, c, FFTW_ESTIMATE);
    fftwf_execute(p3); // Will destroy input array C, but that's okay

    fftwf_destroy_plan(p1);
    fftwf_destroy_plan(p2);
    fftwf_destroy_plan(p3);
    fftwf_free(A);
    fftwf_free(B);
    fftwf_free(C);
}

void complex_mult(const fftwf_complex *A, const fftwf_complex *B, fftwf_complex *product_ft, const int complen) {
    // A * B = (ac - bd) + i(ad + bc)
    for (int j = 0; j < complen; j++) {
        product_ft[j][0] = (A[j][0] * B[j][0]) - (A[j][1] * B[j][1]);
        product_ft[j][1] = (A[j][0] * B[j][1]) + (A[j][1] * B[j][0]);
    }
}

void dd_deconv2(float *kernel_y, const float *data_y, const int length, float *output) {
    // Create flipped point spread function kernel_star
    float *kernel_star = calloc(length, sizeof(float));
    // Create estimate for solution
    float *estimate = calloc(length, sizeof(float));
    // Allocate arrays for convolutions
    float *conv1 = calloc(length, sizeof(float));
    float *conv2 = calloc(length, sizeof(float));
    if (conv1 && conv2 && estimate && kernel_star) {
        for (int i = 0; i < length; i++) {
            kernel_star[i] = kernel_y[length - i - 1];
        }

        for (int i = 0; i < length; i++) {
            estimate[i] = data_y[i];
        }

        const int complen = (length / 2) + 1;
        fftwf_complex *kernel_ft = fftwf_malloc(complen * sizeof(fftwf_complex));
        fftwf_complex *kernel_star_ft = fftwf_malloc(complen * sizeof(fftwf_complex));
        fftwf_complex *estimate_ft = fftwf_malloc(complen * sizeof(fftwf_complex));
        fftwf_complex *conv1_ft = fftwf_malloc(complen * sizeof(fftwf_complex));
        fftwf_complex *product_ft = fftwf_malloc(complen * sizeof(fftwf_complex));

        // fftwf_MEASURE takes a few seconds and overwrites input arrays, so we use ESTIMATE
        fftwf_plan pk = fftwf_plan_dft_r2c_1d(length, kernel_y, kernel_ft, FFTW_ESTIMATE); // for kernel_y
        fftwf_plan pks = fftwf_plan_dft_r2c_1d(length, kernel_star, kernel_star_ft, FFTW_ESTIMATE); // for kernel_star
        fftwf_plan p1 = fftwf_plan_dft_r2c_1d(length, estimate, estimate_ft, FFTW_ESTIMATE); // for estimate
        fftwf_plan p2 = fftwf_plan_dft_r2c_1d(length, conv1, conv1_ft, FFTW_ESTIMATE); // for conv1
        fftwf_plan p1r = fftwf_plan_dft_c2r_1d(length, product_ft, conv1, FFTW_ESTIMATE); // to conv1
        fftwf_plan p2r = fftwf_plan_dft_c2r_1d(length, product_ft, conv2, FFTW_ESTIMATE); // to conv2

        // We only need to do the transforms for kernel and kernel* once
        fftwf_execute(pk);
        fftwf_execute(pks);

        // Perform iterations
        int j = 0;
        float diff = 1.0f;
        while (j < 50 && diff > 0.0001) {
            // Thresholds same as in Python
            // Find new estimate
            fftwf_execute(p1);
            complex_mult(kernel_ft, estimate_ft, product_ft, complen);
            fftwf_execute(p1r);
            for (int k = 0; k < length; k++) {
                if (conv1[k] != 0) {
                    conv1[k] = data_y[k] / conv1[k];
                }
            }
            // cconv2fast(conv1, kernel_star, conv2, length);
            fftwf_execute(p2);
            complex_mult(conv1_ft, kernel_star_ft, product_ft, complen);
            fftwf_execute(p2r);
            for (int k = 0; k < length; k++) {
                conv2[k] = conv2[k] * estimate[k]; // Store new estimate in conv2
            }
            // Find how much the estimate changed
            float sum_diff = 0.0f;
            float sum_est = 0.0f;
            for (int k = 0; k < length; k++) {
                sum_diff += powf((estimate[k] - conv2[k]), 2.0f);
                sum_est += estimate[k];
            }
            diff = sum_diff / sum_est;
            // Set new estimate as estimate
            for (int k = 0; k < length; k++) {
                estimate[k] = conv2[k];
            }
            j++;
        }

        // estimate now contains our deconvolution
        // Normalize and "return"
        float estimate_max = 0.0f;
        for (int i = 0; i < length; i++) {
            if (estimate_max < estimate[i]) {
                estimate_max = estimate[i];
            }
        }
        for (int i = 0; i < length; i++) {
            output[i] = estimate[i] / estimate_max;
        }

        // Free arrays
        free(kernel_star);
        free(estimate);
        free(conv1);
        free(conv2);
        fftwf_destroy_plan(pk);
        fftwf_destroy_plan(pks);
        fftwf_destroy_plan(p1);
        fftwf_destroy_plan(p2);
        fftwf_destroy_plan(p1r);
        fftwf_destroy_plan(p2r);
        fftwf_free(kernel_ft);
        fftwf_free(kernel_star_ft);
        fftwf_free(estimate_ft);
        fftwf_free(conv1_ft);
        fftwf_free(product_ft);
    }
}


void DoubleDecon(const Config *config, Decon *decon) {
    if (decon->massaxis == NULL || decon->massaxisval == NULL) {
        fprintf(stderr, "Error: Config or Decon struct is NULL in DoubleDecon.\n");
        exit(1);
    }

    // Get larger length
    int true_length;
    const int kernel_length = getfilelength(config->kernel);
    const int data_length = decon->mlen;
    if (kernel_length > data_length) {
        true_length = kernel_length;
    } else {
        true_length = data_length;
    }

    // Read in kernel file
    float *kernel_x_init = calloc(true_length, sizeof(float));
    float *kernel_y_init = calloc(true_length, sizeof(float));
    if (kernel_x_init == NULL || kernel_y_init == NULL) {
        fprintf(stderr, "Error allocating memory for kernel_x_init or kernel_y_init in DoubleDecon.\n");
        exit(1);
    }

    readfile(config->kernel, kernel_length, kernel_x_init, kernel_y_init);
    // printf("Kernel file read.\n");

    // Enforce the same sampling on the kernel
    if (kernel_length > 1 && data_length > 1) {
        const float diff = decon->massaxis[1] - decon->massaxis[0]; // kernel sampling needs to match this
        const float kdiff = kernel_x_init[1] - kernel_x_init[0]; // the original kernel sampling
        if (diff != kdiff) {
            const int newlen = (int) ((kernel_x_init[kernel_length - 1] - kernel_x_init[0]) / diff) + 1;
            if (newlen > data_length) {
                true_length = newlen;
            } else {
                true_length = data_length;
            }
        }
    }

    // Read in data (i.e. copy from decon struct)
    float *data_x = calloc(true_length, sizeof(float));
    float *data_y = calloc(true_length, sizeof(float));
    if (data_x == NULL || data_y == NULL) {
        fprintf(stderr, "Error allocating memory for data_x or data_y in DoubleDecon.\n");
        exit(1);
    }

    for (int i = 0; i < data_length; i++) {
        data_x[i] = decon->massaxis[i];
        data_y[i] = decon->massaxisval[i];
    }
    NormMax(data_y, data_length);

    // Integrate or interpolate if necessary...
    float *kernel_x = NULL;
    float *kernel_y = NULL;
    int kernel_length2 = true_length;

    // printf("true length is %d\n", kernel_length2);
    if (kernel_length > 1 && data_length > 1 &&
        (data_x[1] - data_x[0]) > (kernel_x_init[1] - kernel_x_init[0])) {
        printf("Integrating\n");
        kernel_length2 = integrate_dd(kernel_x_init, kernel_y_init, kernel_length, data_x,
                                       true_length, &kernel_x, &kernel_y);
    } else if (kernel_length > 1 && data_length > 1 &&
               (data_x[1] - data_x[0]) < (kernel_x_init[1] - kernel_x_init[0])) {
        printf("Interpolating\n");
        kernel_length2 = interpolate_dd(kernel_x_init, kernel_y_init, kernel_length, data_x,
                                        true_length, &kernel_x, &kernel_y);
    } else {
        printf("Sampling is OK\n");
        kernel_x = kernel_x_init;
        kernel_y = kernel_y_init;
    }

    // ...and find max and normalize
    float max_kernel_y = 0.0f;
    int max_kernel_i = 0;
    for (int i = 0; i < kernel_length2; i++) {
        if (kernel_y[i] > max_kernel_y) {
            max_kernel_y = kernel_y[i];
            max_kernel_i = i;
        }
    }

    for (int i = 0; i < kernel_length2; i++) {
        // Normalize
        kernel_y[i] = kernel_y[i] / max_kernel_y;
    }

    // Pad x-axis for the shorter one (is padding kernel even necessary???)
    if (data_length < true_length) {
        // Pad data_x
        const float diff = data_x[1] - data_x[0];
        float last = data_x[data_length - 1];
        for (int i = data_length; i < true_length; i++) {
            data_x[i] = last + diff;
            last = data_x[i];
        }
    } else if (kernel_length < true_length && kernel_length > 2) {
        // Pad kernel_x
        const float diff = kernel_x[1] - kernel_x[0];
        float last = kernel_x[kernel_length - 1];
        for (int i = kernel_length; i < true_length; i++) {
            kernel_x[i] = last + diff;
            last = kernel_x[i];
        }
    }

    // Prepare kernel
    float *real_kernel_y = calloc(true_length, sizeof(float));
    if (real_kernel_y == NULL) {
        fprintf(stderr, "Error allocating memory for real_kernel_y in DoubleDecon.\n");
        exit(1);
    }

    const int part1_length = true_length - max_kernel_i;
    for (int i = 0; i < part1_length; i++) {
        real_kernel_y[i] = kernel_y[max_kernel_i + i];
    }
    for (int i = 0; i < max_kernel_i; i++) {
        real_kernel_y[part1_length + i] = kernel_y[i];
    }

    // Run Richardson-Lucy deconvolution
    float *doubledec = calloc(true_length, sizeof(float));
    if (doubledec == NULL) {
        fprintf(stderr, "Error allocating memory for doubledec in DoubleDecon.\n");
        exit(1);
    }
    dd_deconv2(real_kernel_y, data_y, true_length, doubledec);

    // printf("Copying results to Decon struct.\n");
    int lb = nearfast(data_x, config->masslb, true_length);
    if (data_x[lb] < config->masslb) lb++;
    const int ub = nearfast(data_x, config->massub, true_length);
    if (data_x[ub] > config->massub) lb--;
    const int write_length = ub - lb + 1;
    if (write_length > decon->mlen) {
        // printf("Warning: new length exceeds previous mlen.\n");
        free(decon->massaxis);
        free(decon->massaxisval);
        decon->massaxis = (float *) calloc(write_length, sizeof(float));
        decon->massaxisval = (float *) calloc(write_length, sizeof(float));
        if (decon->massaxis == NULL || decon->massaxisval == NULL) {
            fprintf(stderr, "Error: Config or Decon struct is NULL in DoubleDecon.\n");
            exit(1);
        }
    }

    // Copy results to the Decon struct
    for (int i = 0; i < write_length; i++) {
        decon->massaxis[i] = data_x[i + lb];
        decon->massaxisval[i] = doubledec[i + lb];
    }

    decon->mlen = write_length;
    // printf("Results copied to Decon.\n");

    // Free memory
    if (kernel_x && kernel_y) {
        free(kernel_x);
        free(kernel_y);
    }

    free(real_kernel_y);
    free(data_x);
    free(data_y);
    free(doubledec);

    free(kernel_x_init);
    free(kernel_y_init);

}

