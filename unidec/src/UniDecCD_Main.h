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
#include "UniDecIM.h"


void blur_it_CD(float *output, const float *input, const int *upinds, const int *loinds, const int length, const float floor);

void setup_blur_z(int *zupind, int *zloind, const float *mzdat, const float *zdat, const int lines, const float adductmass, const float mzranges[4], const int size[3]);

void setup_blur_m(int *mupind, int *mloind, const float *mzdat, const float *zdat, const int lines, const float adductmass, const float mzranges[4], const int size[3], const float molig);

void precompute_fft2D(float *list2, int *size, fftw_complex *out2);

void fftconvolve2D_precomputed(float *corr, float *list1, fftw_complex *out2, int *size, fftw_plan p1, fftw_plan p3, fftw_complex *in1, fftw_complex *out1);

void complex_conjugate(fftw_complex *in, fftw_complex *out, int length);

int run_unidec_CD(int argc, char *argv[], Config config);


#endif

