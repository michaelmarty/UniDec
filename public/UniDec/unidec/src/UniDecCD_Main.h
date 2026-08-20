/*
* UniDecCD_Main.h
*
*  Created on : 16 August 2022 (at home with covid :( )
* Author : Michael.Marty
*/

//
// 
// Copyright 2022 University of Arizona
//
//

#ifndef UNIDEC_CD_MAIN_H
#define UNIDEC_CD_MAIN_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "udcore.h"
//#include <fftw3.h>


void blur_it_CD(float *output, const float *input, const int *upinds, const int *loinds, const int length, const float floor);

void setup_blur_z(int *zupind, int *zloind, const float *mzdat, const float *zdat, const int lines, const float adductmass, const float mzranges[4], const int size[3]);

void setup_blur_m(int *mupind, int *mloind, const float *mzdat, const float *zdat, const int lines, const float adductmass, const float mzranges[4], const int size[3], const float molig);

void precompute_fft2D(const float *list2, const int *size, fftwf_complex *out2);

void fftconvolve2D_precomputed(float *corr, const float *list1, const fftwf_complex *out2, const int *size, fftwf_plan p1, fftwf_plan p3, fftwf_complex *in1, const fftwf_complex *out1);

void complex_conjugate(const fftwf_complex *in, fftwf_complex *out, int length);

int run_unidec_CD(int argc, char *argv[], Config config);


#endif

