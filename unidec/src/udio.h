//
// Created by mm96978 on 7/7/2025.
//

#ifndef UDIO_H
#define UDIO_H

#ifdef __cplusplus
    #define EXTERN extern "C"
#else
    #define EXTERN
#endif

// Define UNIDECLIB_EXPORTS ONLY inside isodeclib.h
#ifdef UNIDEC_BUILD_DLL
    #if defined(_WIN32) || defined(_WIN64)
        #define UDIO_EXPORTS EXTERN __declspec(dllexport)
        // #pragma message "Message: 1 "
    #else
        #define UDIO_EXPORTS EXTERN __attribute__((__visibility__("default")))
// #pragma message "Message: 2 "
    #endif
#else
    #if defined(_WIN32) || defined(_WIN64)
        #define UDIO_EXPORTS EXTERN // __declspec(dllimport)
// #pragma message "Message: 3 "
    #else
        #define UDIO_EXPORTS EXTERN
// #pragma message "Message: 4 "
    #endif
#endif

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
void write1D(char *outfile, char *suffix, const float *array, int length);
void write2D(char *outfile, char *suffix, const float *array1, const float *array2, int length);
void write3D(char *outfile, char *suffix, const float *array1, const float *array2, const float *array3, int length);
void writefile1bin(const char *outstring, int lengthmz, const float *data);
void readmfile(const char *infile, int mfilelen, float *testmasses);
int getfilelength(const char *infile);
int getfilelengthbin(const char *infile, int size, int width);

void ReadInputs(int argc, char* argv[], Config* config, Input* inp);
UDIO_EXPORTS void SetupZtab(const Config config, Input *inp);
void WritePeaks(const Config config, const Decon* decon);
void WriteDecon(const Config config, const Decon* decon, const Input* inp);
void ManualAssign(const float* dataMZ, char* barr, const int* nztab, Config config);

#endif //UDIO_H