//
// Created by Michael Marty on 7/7/2025.
//

#include "udtools.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


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