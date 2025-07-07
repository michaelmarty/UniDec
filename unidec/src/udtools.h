//
// Created by mm96978 on 7/7/2025.
//

#ifndef UDTOOLS_H
#define UDTOOLS_H

#endif //UDTOOLS_H

void floatPrint(const float *array, const int length);
void IntPrint(const int *array, const int length);
float Average(const int length, const float *xarray);
float ndis(const float x, const float y, const float sig);
int mod(const int a, const int b);
float calcmz(const float mass, const int z, const float adductmass);

int nearunsorted(const float *testmasses, const float point, const int lengthtest);
int neartest(const float *testmasses, const float point, const int lengthtest, const float cutoff);
int nearfast(const float *dataMZ, const float point, const int numdat);
