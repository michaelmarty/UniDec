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

void blur_it_IM(int *size, float *blur, float *newblur, int *closetab, int *barr, float csig);
void sumdeltas(int *size, float *deltas, float *newblur);

void bayes(int *size, float*denom, float *blur, float *newblur, float *dataInt, int *barr);
void getranges(int *size, float *blur, float *masstab, float *ccstab, float *range, int *barr);
void makeaxis(float *axis, int length, float min, float binsize);

void readfilemanual(char *infile, int lengthmz, float *array1, float *array2, float *array3, float *array4, float *array5);
void writemfileres(char *outfile, char *suffix, float *array1, float *array2, float *array3, float *array4, float *array5, float *array6, int length);
void writemzgrid(char *outfile, char *suffix, float *mzext, float *dtext, int *ztab, float *blur, int *size);
float errfun(int length, float *dataInt, float *fitdat);

void writezslice(int *size, char *outfile, char *suffix, float *massaxis, float *ccsaxis, int *ztab, float *array, int k);
float nativeCCS(float mass, float fudge, float gasmass);
float bilinearinterpolation(const int *size, int indm, int indc, int k, const float *mzext, const float *dtext, float tempmz, float tempdt, int rawflag, const float *newblur, const float *blur);
float cubicInterpolate_IM(float p[4], float x);
float bicubicInterpolate_IM(float p[4][4], float x, float y);
float bicubicinterpolation(const int *size, int indm, int indc, int k, const float *mzext, const float *dtext, float tempmz, float tempdt, int rawflag, const float *newblur, const float *blur);

void MFileInt(int mfilelen, float *massaxis, float *massaxisval, float *testmasses, float *testmassint, float maaxle);
void MFileCCS(int mfilelen, const float *massaxis, const float *ccsaxis, const float *massccsgrid, const float *testmasses, float *testmassCCSavg, float *testmassCCSstddev, int maaxle, int ccaxle);
void MFileZ(int mfilelen, const float *massaxis, const int *ztab, const float *masszgrid, const float *testmasses, float *testmassZavg, float *testmassZstddev, int maaxle, int numz);
void KillB_IM(const float *IntArray, int *B, const int * size, float intthresh);
void ManualAssign_IM(char *manualfile, const int * size, const float *mzdat, const float * dtdat, const int *ztab, int *barr);
void TwaveError(int twaveflag);

#endif