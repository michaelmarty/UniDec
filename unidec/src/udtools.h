//
// Created by mm96978 on 7/7/2025.
//

#ifndef UDTOOLS_H
#define UDTOOLS_H

#endif //UDTOOLS_H

void floatPrint(const float *array, const int length);
void IntPrint(const int *array, const int length);

float Average(const int length, const float *xarray);
float Max(const float *blur, const int length);
float Min(const float *blur, const int length);
float Sum(const float *blur, const int length);
int ArgMax(const float *blur, const int lengthmz);
void ApplyCutoff(float *array, const float cutoff, const int lengthmz);

float ndis(const float x, const float y, const float sig);
int mod(const int a, const int b);
float clip(const float x, const float cutoff);
float calcmz(const float mass, const int z, const float adductmass);
float mzpeakshape(const float x, const float y, const float sig, const int psfun);

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
