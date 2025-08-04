//
// Created by Michael Marty on 7/7/2025.
//

#include "udtools.h"

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


float WeightedAverage(const int length, const float *xarray, const float *warray)
{
    float temp1 = 0;
    float temp2 = 0;
    for (int i = 0; i<length; i++)
    {
        temp1 += xarray[i] * warray[i];
        temp2 += warray[i];
    }
    return temp1 / temp2;
}

float WeightedAverageInt(const int length, const int *xarray, const float *warray)
{
    float temp1 = 0;
    float temp2 = 0;
    for (int i = 0; i<length; i++)
    {
        temp1 += (float) xarray[i] * warray[i];
        temp2 += warray[i];
    }
    return temp1 / temp2;
}

//Calculate Standard Deviation
float StdDev(const int length, const float *xarray, const float wmean) {
    float temp1 = 0;
    float temp2 = 0;
    for (int i = 0; i < length; i++) {
        temp1 += powf(xarray[i] - wmean, 2.f);
        temp2 += 1;
    }
    if (temp2 == 0) { return 0; }
    return sqrtf(temp1 / temp2);
}


float WeightedStdDev(const int length, const float *xarray, const float *warray, const float wmean)
{
    float temp1 = 0;
    float temp2 = 0;
    for (int i = 0; i<length; i++)
    {
        temp1 += powf(xarray[i] - wmean, 2)*warray[i];
        temp2 += warray[i];
    }
    return sqrtf(temp1 / temp2);
}

float WeightedStdDevInt(const int length, const int *xarray, const float *warray, const float wmean)
{
    float temp1 = 0;
    float temp2 = 0;
    for (int i = 0; i<length; i++)
    {
        temp1 += powf((float) xarray[i] - wmean, 2)*warray[i];
        temp2 += warray[i];
    }
    return sqrtf(temp1 / temp2);
}

//Find the maximum value in an array
float Max(const float *blur, const int length) {
    float blurmax = blur[0];
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


void Sum1D(const int *size, float *sum, const float *array, const int axis)
{
    if (axis == 0) {
        for (int i = 0; i<size[0]; i++)
        {
            float temp = 0;
            for (int j = 0; j<size[1]; j++)
            {

                for (int k = 0; k<size[2]; k++)
                {
                    temp += array[index3D(size[1], size[2], i, j, k)];
                }

            }
            sum[i] = temp;
        }
    }
    else if (axis == 1)
    {
        for (int j = 0; j<size[1]; j++)
        {
            float temp = 0;
            for (int i = 0; i<size[0]; i++)
            {

                for (int k = 0; k<size[2]; k++)
                {
                    temp += array[index3D(size[1], size[2], i, j, k)];
                }

            }
            sum[j] = temp;
        }
    }
    else if (axis == 2)
    {
        for (int k = 0; k<size[2]; k++)
        {
            float temp = 0;
            for (int j = 0; j<size[1]; j++)
            {

                for (int i = 0; i<size[0]; i++)
                {
                    temp += array[index3D(size[1], size[2], i, j, k)];
                }

            }
            sum[k] = temp;
        }
    }
}


void Sum2D(const int *size, float *sum, const float *array, const int axis)
{
    if (axis == 2)
    {
        for (int i = 0; i<size[0]; i++)
        {
            for (int j = 0; j<size[1]; j++)
            {
                float temp = 0;
                for (int k = 0; k<size[2]; k++)
                {
                    temp += array[index3D(size[1], size[2], i, j, k)];
                }
                sum[index2D(size[1], i, j)] = temp;
            }
        }
    }
    else if (axis == 1)
    {
        for (int i = 0; i<size[0]; i++)
        {
            for (int k = 0; k<size[2]; k++)
            {
                float temp = 0;
                for (int j = 0; j<size[1]; j++)
                {
                    temp += array[index3D(size[1], size[2], i, j, k)];
                }
                sum[index2D(size[2], i, k)] = temp;
            }
        }
    }
    else if (axis == 0)
    {
        for (int k = 0; k<size[2]; k++)
        {
            for (int j = 0; j<size[1]; j++)
            {
                float temp = 0;
                for (int i = 0; i<size[0]; i++)
                {
                    temp += array[index3D(size[1], size[2], i, j, k)];
                }
                sum[index2D(size[1], k, j)] = temp;
            }
        }
    }
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

void NormMax(float *data, const int length) {
    const float norm = Max(data, length);
    if (norm > 0) {
        for (int i = 0; i < length; i++) {
            data[i] = data[i] / norm;
        }
    }
}

void NormSum(float *data, const int length) {
    const float norm = Sum(data, length);
    if (norm > 0) {
        for (int i = 0; i < length; i++) {
            data[i] = data[i] / norm;
        }
    }
}

void Normalize(const int length, float *data, const float norm)
{
    if (norm == 0) { printf("Tried to Normalize by Zero! Aborting!\n"); exit(50); }
    #pragma omp parallel for schedule(auto)
    for (int i = 0; i<length; i++)
    {
        data[i] = data[i] / norm;
    }
}


int GetSize0(const float *array, const int lines)
{
    int cnt = 1;
    float test = array[0];
    for (int i = 1; i<lines; i++) //look for the number of times dH decrements
    {
        if (array[i] > test)
        {
            test = array[i];
            cnt++;
        }
    }
    return cnt;
}

int GetSize1(const float *array, const int lines)
{
    int cnt = 1;
    int flag = 0;
    const float test = array[0];
    for (int i = 1; i<lines; i++)
    {
        if (array[i] > test && flag == 0)
        {
            cnt++;
        }
        else
        {
            flag = 1;
        }
    }
    return cnt;
}

// From two lists of data, extract the unique values for the first dimension (mz) and the second dimension (dt)
void PullXY(float *mzext, float *dtext, const float *mzdat, const float *dtdat, const int *size)
{
    //Extract mz data
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i<size[0]; i++)
    {
        mzext[i] = mzdat[i*size[1]];
    }
    //Extract dt data
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i<size[1]; i++)
    {
        dtext[i] = dtdat[i];
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


float Gaus(const float x, const float x0, const float sig)
{
    return expf(-powf(x - x0, 2.f) / (2.f*sig*sig));
}


//function to return a gaussian at a specific point
float Lorentz(const float x, const float x0, const float sig)
{
    return powf(sig / 2, 2) / (powf(x - x0, 2) + powf(sig / 2, 2));
}

float SplitGL(const float x, const float y, const float sig)
{
    if (y<x)
    {
        return expf(-(powf(x - y, 2)) / (2 * sig*sig*0.180337f));
    }
    else
    {
        return (sig / 2)*(sig / 2) / (powf((x - y), 2) + powf((sig / 2), 2));
    }
}


//function to return a 2D peak value at a given point
float PeakDist(const float x1, const float x2, const float y1, const float y2, const float sig1, const float sig2, const int psfun, const int zpsfun)
{
    float d2 = 0;
    float d1 = 0;
    if (sig2 == 0) {
        if (y1 == y2) { d2 = 1; }
        else { d2 = 0; }
    }
    else {
        if (zpsfun == 0) {
            d2 = Gaus(y1, y2, sig2);
        }
        if (zpsfun == 1) {
            d2 = Lorentz(y1, y2, sig2);
        }
        if (psfun == 2) {
            d2 = SplitGL(-1*y1, -1*y2, sig2);
        }
    }

    if (sig1 == 0) {
        if (x1 == x2) { d1 = 1; }
        else { d1 = 0; }
    }
    else {
        if (psfun == 0) {
            d1= Gaus(x1, x2, sig1);
        }
        if (psfun == 1) {
            d1 = Lorentz(x1, x2, sig1);
        }
        if (psfun == 2) {
            d1 = SplitGL(x1, x2, sig1);
        }
    }

    return d1 * d2;
}



void MakeKernel2D(float *peak, const int *Size, const float *dH, const float *dC, const float sig1,
    const float sig2, const int psfun, const int zpsfun)
{
    const float hmax = 2 * dH[Size[0] - 1] - dH[Size[0] - 2];
    const float cmax = 2 * dC[Size[1] - 1] - dC[Size[1] - 2];
    const float hmin = dH[0];
    const float cmin = dC[0];

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i<Size[0]; i++)
    {
        for (int j = 0; j<Size[1]; j++)
        {
            const float Hval = dH[i];
            const float Cval = dC[j];

            //put the peakj in the four corners
            const float v1 = PeakDist(hmin, Hval, cmin, Cval, sig1, sig2, psfun, zpsfun);
            const float v2 = PeakDist(hmax, Hval, cmax, Cval, sig1, sig2, psfun, zpsfun);
            const float v3 = PeakDist(hmin, Hval, cmax, Cval, sig1, sig2, psfun, zpsfun);
            const float v4 = PeakDist(hmax, Hval, cmin, Cval, sig1, sig2, psfun, zpsfun);

            const float val = v1 + v2 + v3 + v4;
            peak[index2D(Size[1], i, j)] = val; //store peak in complete peak list
        }
    }
    //printf("mzmax: %f\n",hmax);
    //printf("dtmax: %f\n",cmax);
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

float euclid(const float a, const float b, const float c, const float d)
{
    return sqrtf((a - b)*(a - b) + (c - d)*(c - d));
}

// Simple function to calculate mz from mass, z, and adductmass
float calcmz(const float mass, const int z, const float adductmass) {
    return (mass + (float) z * adductmass) / (float) z;
}

// Simple function to calculate mass from mz, z, and adductmass
float calcmass(const float mz, const int z, const float adductmass) {
    return (float) z * (mz - adductmass);
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

// Compare function for qsort
int compare_function(const void *a, const void *b) {
    const float *x = (float *) a;
    const float *y = (float *) b;
    if (*x < *y) { return -1; }
    if (*x > *y) { return 1; }
    return 0;
}

//Average native charge state from Champ
float nativecharge(const float mass, const float fudge) {
    return 0.0467f * powf(mass, 0.533f) + fudge;
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
        const float difftest = fabsf(point - testmasses[i]);
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

//Fast way of finding the nearest data point index in an ordered list.
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

// Return the nearest data point value in an ordered list.
float nearfastval(const float *dataMZ, const float point, const int numdat)
{
    const int index = nearfast(dataMZ, point, numdat);
    return dataMZ[index];
}


int nearint(const int *array, const int point, const int numdat)
{
    int temp = -1;
    for (int i = 0; i<numdat; i++)
    {
        if (array[i] == point) { temp = i; }
    }
    if (temp == -1) { printf("Problem finding integer in list\n"); }
    return temp;
}


//Slow nearest unsorted for 2D data. Returns the index of the closest point in the testmasses and testtimes arrays.
int nearunsorted2D(const float *testmasses, const float *testtimes, const float mz, const float dt, const int lengthtest)
{
    float minval = euclid(mz, testmasses[0], dt, testtimes[0]);
    int pos = 0;
    for (int i = 1; i<lengthtest; i++)
    {
        const float difftest = euclid(mz, testmasses[i], dt, testtimes[i]);
        if (difftest<minval)
        {
            minval = difftest;
            pos = i;
        }
    }
    return pos;
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