//
// Created by mm96978 on 7/7/2025.
//

#ifndef UDCORE_H
#define UDCORE_H


#include "udstruct.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "udtools.h"
#include <fftw3.h>
#include "udio.h"

#ifdef UNIDEC_BUILD_DLL
    #ifdef _WIN32
        #define UDCORE_EXPORTS EXTERN __declspec(dllexport)
    #else
        #define UDCORE_EXPORTS EXTERN __attribute__((__visibility__("default")))
    #endif
#else
    #define UDCORE_EXPORTS EXTERN
#endif


void check_ratios(Config config, Input inp, char * barr, int numclose,
              const int *__restrict closeind,
              const float *__restrict blur);

void adjust_ratios(Config config, char * barr, int numclose,
              const int * closeind, float * blur);

void blur_it(const IntraDecon *intra, float * newblur, const float * blur);

void ApplyBlur(const Config *config, const IntraDecon *intra, float *output, const float *input);
void blur_it_mean(const IntraDecon *intra, float * newblur, const float * blur, float zerolog);
// Shared-grid setup reserves one ln+1 log buffer per active worker and maps missing neighbors to ln.
void blur_it_mean_three_batches(const IntraDecon *intra, float *output, const float *input,
                                int batch_count, float zerolog);

void blur_it_hybrid1(const IntraDecon *intra, int lengthmz, int numz, float * newblur, const float * blur,
                     float zerolog);

void blur_it_hybrid2(const IntraDecon *intra, int lengthmz, int numz, float * newblur, const float * blur,
                     float zerolog);

void blur_baseline(float *baseline, int lengthmz, const float *dataMZ, float mzsig, int mult,
                   int filterwidth);

void midblur_baseline(float *baseline, int lengthmz, const float *dataMZ, float mzsig, int mult);

void convolve_simp(int lengthmz, int maxlength, const int *__restrict starttab,
                   const int *__restrict endtab, const float *__restrict mzdist,
                   const float *__restrict deltas, float *__restrict denom, int speedyflag);

void sum_deltas(int lengthmz, int numz, const float *__restrict blur, const char *__restrict barr,
                float *deltas);

void apply_ratios(const int lengthmz, const int numz, const float *__restrict blur, const char *__restrict barr,
                  const float *__restrict denom, float *blur2);

void deconvolve_baseline(const int lengthmz, const float *dataMZ, const float *dataInt, float *baseline,
                         const float mzsig);

float deconvolve_iteration_speedy(const Config *config, Decon *decon, IntraDecon *intra,
                                   const float *dataMZ);

UDCORE_EXPORTS void softargmax(float *blur, const int lengthmz, const int numz, const float beta);
void softmax_peakwidth(const Config config, const Decon decon, float *blur, const char * barr, const float beta);

UDCORE_EXPORTS void point_smoothing(float *blur, float *scratch, float *sums, const char *barr,
                                   const int lengthmz, const int numz, const int width);
void point_smoothing_to_scratch(const float *blur, float *scratch, float *sums, const char *barr,
                                const int lengthmz, const int numz, const int width);

void highest_n_chargestates(float *blur, float *scratch, int lengthmz, int numz, int n, float zcutpercent);
void clip_minor_chargestates(float *blur, int lengthmz, int numz, float zcutoff, float zcutpercent);
void suppression_satelite(float *blur, float *scratch, int lengthmz, int numz, int n);
void suppression_harmonic(float *blur, float *scratch, int lengthmz, int numz, const int *ztab);
int setup_suppression_ztab(int *ztab, const float *zaxis, int numz);
UDCORE_EXPORTS void apply_suppressions(float *blur, float *scratch, int lengthmz, int numz,
                                      int satellite_n, int harmonic, const int *ztab, int top_n,
                                      float top_x, float suppression_percent);

float getfitdatspeedy(float *fitdat, const float *blur, const int lengthmz, const int numz,
                      const int maxlength, const float maxint,
                      const int *starttab,
                      const int *endtab, const float *mzdist, const int speedyflag);

float errfunspeedy(Config config, Decon decon, const float *dataInt, float *rsquared, float *fitdatnorm);

void CalcMasses(const Config *config, Input *inp);
void TestMassListWindowed(const Config config, Input *inp);
void TestMassListLimit(const Config config, Input *inp, const int *testmasspos);
void TestMass(const Config config, Input *inp);
void SetLimits(const Config config, Input *inp);
void ignorezeros(char *barr, const float *dataInt, const int lengthmz, const int numz);
void KillB(const float *I, char *B, const float intthresh, const int lengthmz, const int numz);

void MakeSparseBlur(const int numclose, char *barr, const int *closezind,
                    const int *closemind, int *closeind, const float *closeval, float *closearray, const Config config, const Input *inp);
void MakePeakShape2D(const Config config, Decon *decon, const Input *inp, int makereverse,
                     const int inflateflag);
void MakePeakShape1D(const Config config, Decon * decon, const float *dataMZ, int makereverse, const int inflateflag, const int silent);
int SetStartsEnds(const Config config, const Input *inp, int *starttab, int *endtab);
int SetUpPeakShape(Config config, Input inp, Decon *decon, const int silent, const int verbose);
int PrepareBlur(Config config, Input inp, IntraDecon *intra, int silent);
void SetUpBlur(Config config, Input inp, IntraDecon *intra, const int silent);
float Reconvolve(const Config config, Decon *decon, const char *barr);
void charge_scaling(float *blur, const int *nztab, const int lengthmz, const int numz);

void IntegrateTransform(const Config config, Decon *decon, const float *mtab, float massmax, float massmin);
void InterpolateTransform(const Config config, Decon *decon, const Input *inp);
void SmartTransform(const Config config, Decon *decon, const Input *inp);


void DoubleDecon(const Config *config, Decon *decon);

#endif //UDCORE_H
