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
//#include <complex.h>
//#include <omp.h>



void readfilemanual(char *infile, int lengthmz, float *array1, float *array2, float *array3, float *array4, float *array5);
int GetSize0(float *array, int lines);
int GetSize1(float *array, int lines);
void Extract(float *mzext, float *dtext, float *mzdat, float *dtdat, int *size);
float Gaus(float x, float x0, float sig);
float Lorentz(float x, float x0, float sig);
float SplitGL(float x, float y, float sig);
float Peak(const float x1, const float x2, const float y1, const float y2, const float sig1, const float sig2, const int psfun, const int zpsfun);
void GetPeaks(float *peak, int *Size, float *dH, float *dC, float sig1, float sig2, int psfun, int zpsfun);
void fftconvolve2D(float *corr, float *list1, float *list2, int *size);
void convolve3D(float *out, float *in, float *peakshape, int *size);
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
float calcMass(float mz, int z, float adductmass);
float calcMz(float mass, int z, float adductmass);
int nearint(int *array, int point, int numdat);
int nearfastpoint(float *dataMZ, float point, int numdat);
float euclid(float a, float b, float c, float d);
int nearunsorted_IM(float *testmasses, float *testtimes, float mz, float dt, int lengthtest);
void sum2D(int *size, float *sum, float *array, int axis);
void blur_it_IM(int *size, float *blur, float *newblur, int *closetab, int *barr, float csig);
void sumdeltas(int *size, float *deltas, float *newblur);
void sum1D(int *size, float *sum, float *array, int axis);
void bayes(int *size, float*denom, float *blur, float *newblur, float *dataInt, int *barr);
void getranges(int *size, float *blur, float *masstab, float *ccstab, float *range, int *barr);
void makeaxis(float *axis, int length, float min, float binsize);
void write1D(char *outfile, char *suffix, float *array, int length);
void write2D(char *outfile, char *suffix, float *array1, float *array2, int length);
void write3D(char *outfile, char *suffix, float *array1, float *array2, float *array3, int length);
void writemfileres(char *outfile, char *suffix, float *array1, float *array2, float *array3, float *array4, float *array5, float *array6, int length);
void writemzgrid(char *outfile, char *suffix, float *mzext, float *dtext, int *ztab, float *blur, int *size);
float errfun(int length, float *dataInt, float *fitdat);
float getmax(int length, float *data);
void normalize(int length, float *data, float norm);
void writezslice(int *size, char *outfile, char *suffix, float *massaxis, float *ccsaxis, int *ztab, float *array, int k);
float nativeCCS(float mass, float fudge, float gasmass);
float bilinearinterpolation(int *size, int indm, int indc, int k, float *mzext, float *dtext, float tempmz, float tempdt, int rawflag, float *newblur, float *blur);
float cubicInterpolate_IM(float p[4], float x);
float bicubicInterpolate_IM(float p[4][4], float x, float y);
float bicubicinterpolation(int *size, int indm, int indc, int k, float *mzext, float *dtext, float tempmz, float tempdt, int rawflag, float *newblur, float *blur);
float WeightedAverage(int length, float *xarray, float *warray);
float WeightedAverageInt(int length, int *xarray, float *warray);
float StdDev2(int length, float *xarray, float *warray, float wmean);
float StdDev2Int(int length, int *xarray, float *warray, float wmean);
void MFileInt(int mfilelen, float *massaxis, float *massaxisval, float *testmasses, float *testmassint, float maaxle);
void MFileCCS(int mfilelen, float *massaxis, float *ccsaxis, float *massccsgrid, float *testmasses, float *testmassCCSavg, float *testmassCCSstddev, int maaxle, int ccaxle);
void MFileZ(int mfilelen, float *massaxis, int *ztab, float *masszgrid, float *testmasses, float *testmassZavg, float *testmassZstddev, int maaxle, int numz);
void KillB_IM(float *IntArray, int *B, int * size, float intthresh);
void ManualAssign_IM(char *manualfile, int * size, float *mzdat, float * dtdat, int *ztab, int *barr);
void TwaveError(int twaveflag);

static int index4D(int *size, int r, int c, int d, int f)
{
    return r*size[1] * size[2] * size[3] + c*size[2] * size[3] + d*size[3] + f;
}

#endif