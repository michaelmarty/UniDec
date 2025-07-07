//
// Created by Michael Marty on 7/7/2025.
//

#include "udtools.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <string.h>

// ............
//
// Print Functions
//
// ............


//Print a float array
void floatPrint(const float *array, const int length) {
    for (int i = 0; i < length; i++) {
        printf("%f\n", array[i]);
    }
}

//Print an int array
void IntPrint(const int *array, const int length) {
    for (int i = 0; i < length; i++) {
        printf("%d\n", array[i]);
    }
}

// ............
//
// Array Functions
//
// ............

//Calculate Average
float Average(const int length, const float *xarray) {
    float temp1 = 0;
    float temp2 = (float) length;
    for (int i = 0; i < length; i++) {
        temp1 += xarray[i];
    }
    if (temp2 == 0) { return 0; }
    return temp1 / temp2;
}

//Find the maximum value in an array
float Max(const float *blur, const int length) {
    float blurmax = 0;
    for (int i = 0; i < length; i++) {
        if (blur[i] > blurmax) {
            blurmax = blur[i];
        }
    }
    return blurmax;
}

//Find the minimum value in an array
float Min(const float *blur, const int length) {
    float blurmin = blur[0];
    for (int i = 0; i < length; i++) {
        if (blur[i] < blurmin) {
            blurmin = blur[i];
        }
    }
    return blurmin;
}


//Find the Sum in an array
float Sum(const float *blur, const int length) {
    float sum = 0;
    for (int i = 0; i < length; i++) {
        sum += blur[i];
    }
    return sum;
}

//Find the index of the maximum value in an array
int ArgMax(const float *blur, const int lengthmz) {
    float max = blur[0];
    int pos = 0;
    for (int i = 0; i < lengthmz; i++) {
        if (blur[i] > max) {
            max = blur[i];
            pos = i;
        }
    }
    return pos;
}

//Apply a cutoff to a 1D array
void ApplyCutoff(float *array, const float cutoff, const int lengthmz) {
    //#pragma omp parallel for private (i,j), schedule(auto)
    for (int i = 0; i < lengthmz; i++) {
        if (array[i] < cutoff) { array[i] = 0; }
    }
}

// ............
//
// Math Functions
//
// ............

// Function for normal distribution
float ndis(const float x, const float y, const float sig) {
    if (sig == 0) { return 0; }
    return 1.0f / (sig * 2.50663f) * expf(-(powf(x - y, 2.0f)) / (2.0f * sig * sig));
}

//Actual Modulus operator rather than Remainder operator %
int mod(const int a, const int b) {
    const int r = a % b;
    return r < 0 ? r + b : r;
}

// Clip a number to a cutoff value
float clip(const float x, const float cutoff) {
    if (x > cutoff) { return x; }
    return 0;
}

// Simple function to calculate mz from mass, z, and adductmass
float calcmz(const float mass, const int z, const float adductmass) {
    return (mass + (float) z * adductmass) / (float) z;
}


//Function for defining m/z peak shape. Other peak shape functions could be easily added here.
float mzpeakshape(const float x, const float y, const float sig, const int psfun) {
    if (sig == 0) {
        printf("Error: mzpeakshape sigma is 0\n");
        exit(103);
    }
    float result;
    if (psfun == 0) {
        result = expf(-(powf(x - y, 2.f)) / (2.f * sig * sig));
    } else if (psfun == 1) {
        result = powf((sig / 2.f), 2.f) / (powf((x - y), 2.f) + powf((sig / 2.f), 2.f));
    } else if (psfun == 2) {
        if (y < x) {
            result = expf(-(powf(x - y, 2.f)) / (2.f * sig * sig * 0.180337f));
        } else {
            result = (sig / 2.f) * (sig / 2.f) / (powf((x - y), 2.f) + powf((sig / 2), 2.f));
        }
    } else {
        printf("Invalid Peak Function");
        exit(14);
    }
    return result;
}

// ............
//
// Nearest Neighbor Functions
//
// ............

//Slow nearest unsorted.
int nearunsorted(const float *testmasses, const float point, const int lengthtest) {
    float minval = fabsf(point - testmasses[0]);
    int pos = 0;
    for (int i = 1; i < lengthtest; i++) {
        float difftest = fabsf(point - testmasses[i]);
        if (difftest < minval) {
            minval = difftest;
            pos = i;
        }
    }
    return pos;
}

//Slow way to test if two points are within a cutoff distance of each other. Works on unsorted list.
int neartest(const float *testmasses, const float point, const int lengthtest, const float cutoff) {
    float minval = fabsf(point - testmasses[0]);
    float val = testmasses[0];
    for (int i = 0; i < lengthtest; i++) {
        const float difftest = fabsf(point - testmasses[i]);
        if (difftest < minval) {
            minval = difftest;
            val = testmasses[i];
        }
    }
    int test = 0;
    if (fabsf(point - val) < cutoff) {
        test = 1;
    }
    return test;
}

//Fast way of finding the nearest data point in an ordered list.
int nearfast(const float *dataMZ, const float point, const int numdat) {
    int start = 0;
    int length = numdat - 1;
    int end = 0;
    int diff = length - start;
    while (diff > 1) {
        if (point < dataMZ[start + (length - start) / 2]) {
            length = start + (length - start) / 2;
        } else if (point == dataMZ[start + (length - start) / 2]) {
            end = start + (length - start) / 2;
            return end;
        } else if (point > dataMZ[start + (length - start) / 2]) {
            start = start + (length - start) / 2;
        }
        diff = length - start;
    }
    if (fabsf(point - dataMZ[start]) >= fabsf(point - dataMZ[length])) {
        end = length;
    } else {
        end = start;
    }
    return end;
}

// .............
//
// Spline and Interpolate Functions
//
// ............


//Perform a linear interpolation. I've left the code for other interpolation functions below but they don't seem to matter.
float LinearInterpolate(const float y1, const float y2, const float mu) {
    return (y1 * (1 - mu) + y2 * mu);
}

// //Calculate the position of a point between two points in a linear interpolation.
float LinearInterpolatePosition(const float x1, const float x2, const float x) {
    if (x2 - x1 == 0) { return 0; }
    return (x - x1) / (x2 - x1);
}

// Perform cubic interpolation
float CubicInterpolate(const float y0, const float y1, const float y2, const float y3, const float mu) {
    const float mu2 = mu * mu;
    const float a0 = y3 - y2 - y0 + y1;
    const float a1 = y0 - y1 - a0;
    const float a2 = y2 - y0;
    const float a3 = y1;
    return (a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3);
}

// Perform cubic spline interpolation
float CRSplineInterpolate(const float y0, const float y1, const float y2, const float y3, const float mu) {
    const float mu2 = mu * mu;
    const float a0 = -0.5f * y0 + 1.5f * y1 - 1.5f * y2 + 0.5f * y3;
    const float a1 = y0 - 2.5f * y1 + 2.0f * y2 - 0.5f * y3;
    const float a2 = -0.5f * y0 + 0.5f * y2;
    const float a3 = y1;
    return (a0 * mu * mu2 + a1 * mu2 + a2 * mu + a3);
}