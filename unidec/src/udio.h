//
// Created by mm96978 on 7/7/2025.
//

#ifndef UDIO_H
#define UDIO_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdio.h>
#include "udtools.h"
#include "udstruct.h"
#include "h5io.h"

void readfile(const char *infile, int lengthmz, float *dataMZ, float *dataInt);
void readfile3(char *infile, int lengthmz, float *array1, float *array2, float *array3);
void readfile3bin(const char *infile, int lengthmz, float *array1, float *array2, float *array3);
void readfile1bin(const char *infile, int lengthmz, float *data);
void writefile1bin(const char *outstring, int lengthmz, const float *data);
void readmfile(const char *infile, int mfilelen, float *testmasses);
int getfilelength(const char *infile);
int getfilelengthbin(const char *infile, int size, int width);

void ReadInputs(int argc, char* argv[], Config* config, Input* inp);
void SetupZtab(const Config config, Input *inp);
void WritePeaks(const Config config, const Decon* decon);
void WriteDecon(const Config config, const Decon* decon, const Input* inp);
void ManualAssign(const float* dataMZ, char* barr, const int* nztab, Config config);

#endif //UDIO_H