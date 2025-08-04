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


void blur_it(const int lengthmz,
             const int numz,
             const int numclose,
             const int *__restrict closeind,
             const float *__restrict closearray,
             float *__restrict newblur,
             const float *__restrict blur,
             const char *__restrict barr);

void blur_it_mean(const int lengthmz,
                  const int numz,
                  const int numclose,
                  const int *__restrict closeind,
                  float *__restrict newblur,
                  const float *__restrict blur,
                  const char *__restrict barr,
                  const float *__restrict closearray,
                  const float zerolog);

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
                     const float zerolog);

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
                     const float zerolog);

void blur_baseline(float *baseline, const int lengthmz, const float *dataMZ, const float mzsig, int mult,
                   int filterwidth);

void midblur_baseline(float *baseline, const int lengthmz, const float *dataMZ, const float mzsig, int mult);

void convolve_simp(const int lengthmz, const int maxlength, const int *starttab, const int *endtab, const float *mzdist,
                   const float *deltas, float *denom, const int speedyflag);

void sum_deltas(const int lengthmz, const int numz, const float *__restrict blur, const char *__restrict barr,
                float *deltas);

void apply_ratios(const int lengthmz, const int numz, const float *__restrict blur, const char *__restrict barr,
                  const float *__restrict denom, float *blur2);

void deconvolve_baseline(const int lengthmz, const float *dataMZ, const float *dataInt, float *baseline,
                         const float mzsig);

float deconvolve_iteration_speedy(const int lengthmz, const int numz, const int maxlength, const float *__restrict blur,
                                  float *__restrict blur2,
                                  const char *__restrict barr, const int aggressiveflag,
                                  const float *__restrict dataInt,
                                  const int *starttab, const int *endtab, const float *mzdist, const float *rmzdist,
                                  const int speedyflag, const int baselineflag, float *baseline,
                                  float *noise, const float mzsig, const float *dataMZ, const int filterwidth,
                                  const float psig);

void softargmax(float *blur, const int lengthmz, const int numz, const float beta);
void point_smoothing(float *blur, const char *barr, const int lengthmz, const int numz, const int width);
float getfitdatspeedy(float *fitdat, const float *blur, const int lengthmz, const int numz,
                      const int maxlength, const float maxint,
                      const int *starttab,
                      const int *endtab, const float *mzdist, const int speedyflag);

float errfunspeedy(Config config, Decon decon, const char *barr, const float *dataInt, const int maxlength,
                   const int *starttab, const int *endtab,
                   const float *mzdist, float *rsquared);

void CalcMasses(const Config *config, Input *inp);
void TestMassListWindowed(const Config config, Input *inp);
void TestMassListLimit(const Config config, Input *inp, const int *testmasspos);
void TestMass(const Config config, Input *inp);
void SetLimits(const Config config, Input *inp);
void ignorezeros(char *barr, const float *dataInt, const int lengthmz, const int numz);
void KillB(const float *I, char *B, const float intthresh, const int lengthmz, const int numz);

void MakeSparseBlur(const int numclose, char *barr, const int *closezind,
                    const int *closemind, int *closeind, const float *closeval, float *closearray, const Config config, const Input *inp);
void MakePeakShape2D(const Config config, Decon *decon, int maxlength, const float *dataMZ, int makereverse,
                     const int inflateflag);
void MakePeakShape1D(const Config config, Decon * decon, const float *dataMZ, int makereverse, const int inflateflag);
int SetStartsEnds(const Config config, const Input *inp, int *starttab, int *endtab);
int SetUpPeakShape(Config config, Input inp, Decon *decon, const int silent, const int verbose);
float Reconvolve(const Config config, const int maxlength, Decon *decon, const char *barr);
void charge_scaling(float *blur, const int *nztab, const int lengthmz, const int numz);

void IntegrateTransform(const Config config, Decon *decon, const float *mtab, float massmax, float massmin);
void InterpolateTransform(const Config config, Decon *decon, const Input *inp);
void SmartTransform(const Config config, Decon *decon, const Input *inp);


void DoubleDecon(const Config *config, Decon *decon);

#endif //UDCORE_H
