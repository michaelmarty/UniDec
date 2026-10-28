/*
* UniDecIM.h
*
*  Created on : 23 December 2016
* Author : Michael.Marty
*/
//
// Copyright 2015 University of Oxford
// Copyright 2016 University of Arizona
//
//
#pragma once

#if HAVE_MALLOC_H
#include <malloc.h>
#else
#include <stdlib.h>
#endif


#ifndef UNIDEC_IM_H
#define UNIDEC_IM_H

#include <stdio.h>
#include <math.h>
#include <fftw3.h>

#include "udtools.h"
#include "udio.h"

void fftconvolve2D(float *corr, const float *list1, const float *list2, const int *size);
void convolve3D(float *out, const float *in, const float *peakshape, const int *size);

float calcCCS(const float mass, const int z, const float dt, const float ccsconst, const float hmass, const float to, const int type);
float calcCCSTwaveLog(float mass, int z, float dt, float tcal1, float tcal2, float hmass, float edc);
float calcCCSTwaveLinear(float mass, int z, float dt, float tcal1, float tcal2, float hmass, float edc);
float calcCCSTwavePower(float mass, int z, float dt, float tcal1, float tcal2, float hmass, float edc);
float calcDtTwaveLog(float mass, int z, float ccs, float tcal1, float tcal2, float hmass, float edc);
float calcDtTwaveLinear(float mass, int z, float ccs, float tcal1, float tcal2, float hmass, float edc);
float calcDtTwavePower(float mass, int z, float ccs, float tcal1, float tcal2, float hmass, float edc);
float calcDt(const float mass, const int z, const float ccs, const float ccsconst, const float hmass, const float to, const int type);
float calcCCSSLIMpoly3(float mass, int z, float dt, float tcal1, float tcal2, float tcal3, float tcal4, float hmass, float edc);
float calcDtSLIMpoly3(float mass, int z, float ccs, float tcal1, float tcal2, float tcal3, float tcal4, float hmass, float edc);
float calcCCSSLIMpoly2(float mass, int z, float dt, float tcal1, float tcal2, float tcal3, float hmass, float edc);
float calcDtSLIMpoly2(float mass, int z, float ccs, float tcal1, float tcal2, float tcal3, float hmass, float edc);

void blur_it_IM(const int *size, float *blur, float *newblur, const int *closetab, const int *barr, float csig);
void sumdeltas(const int *size, float *deltas, const float *newblur);

void bayes(const int *size, const float*denom, float *blur, const float *newblur, const float *dataInt, const int *barr);
void getranges(const int *size, const float *blur, const float *masstab, const float *ccstab, float *range, const int *barr);
void makeaxis(float *axis, int length, float min, float binsize);

void readfilemanual(char *infile, int lengthmz, float *array1, float *array2, float *array3, float *array4, float *array5);
void writemfileres(char *outfile, char *suffix, const float *array1, const float *array2, const float *array3, const float *array4, const float *array5, const float *array6, int length);
void writemzgrid(char *outfile, char *suffix, const float *blur, const int *size);
float errfun(int length, const float *dataInt, const float *fitdat);

void writezslice(const int *size, char *outfile, char *suffix, const int *ztab, const float *array, int k);
float nativeCCS(float mass, float fudge, float gasmass);
float bilinearinterpolation(const int *size, int indm, int indc, int k, const float *mzext, const float *dtext, float tempmz, float tempdt, int rawflag, const float *newblur, const float *blur);
float cubicInterpolate_IM(float p[4], float x);
float bicubicInterpolate_IM(float p[4][4], float x, float y);
float bicubicinterpolation(const int *size, int indm, int indc, int k, const float *mzext, const float *dtext, float tempmz, float tempdt, int rawflag, const float *newblur, const float *blur);

void MFileInt(int mfilelen, const float *massaxis, const float *massaxisval, const float *testmasses, float *testmassint, int maaxle);
void MFileCCS(int mfilelen, const float *massaxis, const float *ccsaxis, const float *massccsgrid, const float *testmasses, float *testmassCCSavg, float *testmassCCSstddev, int maaxle, int ccaxle);
void MFileZ(int mfilelen, const float *massaxis, const int *ztab, const float *masszgrid, const float *testmasses, float *testmassZavg, float *testmassZstddev, int maaxle, int numz);
void KillB_IM(const float *IntArray, int *B, const int * size, float intthresh);
void ManualAssign_IM(char *manualfile, const int * size, const float *mzdat, const float * dtdat, const int *ztab, int *barr);
void TwaveError(int twaveflag);

#endif