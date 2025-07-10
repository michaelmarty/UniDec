//
// Created by mm96978 on 7/7/2025.
//

#ifndef UDTOOLS_H
#define UDTOOLS_H
#include <stdlib.h>
#include <stdio.h>
#include <math.h>



void floatPrint(const float *array, const int length);
void IntPrint(const int *array, const int length);

float Average(const int length, const float *xarray);
float StdDev(const int length, const float *xarray, const float wmean);
float Max(const float *blur, const int length);
float Min(const float *blur, const int length);
float Sum(const float *blur, const int length);
int ArgMax(const float *blur, const int lengthmz);
void ApplyCutoff(float *array, const float cutoff, const int lengthmz);
void NormMax(float *data, const int length);
void NormSum(float *data, const int length);
void Normalize(const int length, float *data, const float norm);

float ndis(const float x, const float y, const float sig);
int mod(const int a, const int b);
float clip(const float x, const float cutoff);
float calcmz(const float mass, const int z, const float adductmass);
float calcmass(const float mz, const int z, const float adductmass);
float mzpeakshape(const float x, const float y, const float sig, const int psfun);
int compare_function(const void *a, const void *b);
float nativecharge(const float mass, const float fudge);

int nearunsorted(const float *testmasses, const float point, const int lengthtest);
int neartest(const float *testmasses, const float point, const int lengthtest, const float cutoff);
int nearfast(const float *dataMZ, const float point, const int numdat);

float LinearInterpolate(const float y1, const float y2, const float mu);
float LinearInterpolatePosition(const float x1, const float x2, const float x);
float CubicInterpolate(const float y0, const float y1, const float y2, const float y3, const float mu);


inline int index2D(const int ncols, const int r, const int c) {
    return r * ncols + c;
}

inline int indexmod(const int length, const int r, const int c) {
    return mod((c - r), length);
}

inline int index3D(const int ncols, const int nrows, const int r, const int c, const int d) {
    return r * ncols * nrows + c * nrows + d;
}

static int index4D(const int *size, const int r, const int c, const int d, const int f)
{
    return r*size[1] * size[2] * size[3] + c*size[2] * size[3] + d*size[3] + f;
}

inline int fixk(int k, int lengthmz) {
    k = abs(k);
    if (k >= lengthmz) { k = 2 * lengthmz - k - 2; }
    //if (k < 0) { k = 0; }
    return k;
}


#endif //UDTOOLS_H