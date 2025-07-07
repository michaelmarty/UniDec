//
// Created by Michael Marty on 7/7/2025.
//

#include "udtools.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include <string.h>


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

// Simple function to calculate mz from mass, z, and adductmass
float calcmz(const float mass, const int z, const float adductmass) {
    return (mass + (float) z * adductmass) / (float) z;
}


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

