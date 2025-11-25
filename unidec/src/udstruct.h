//
// Created by mm96978 on 7/7/2025.
//

#ifndef UDSTRUCT_H
#define UDSTRUCT_H


#ifdef __cplusplus
    #define EXTERN extern "C"
#else
    #define EXTERN
#endif

// 🔥 Define UDSTRUCT_EXPORTS ONLY inside isodeclib.h
#ifdef UNIDEC_BUILD_DLL
    #if defined(_WIN32) || defined(_WIN64)
        #define UDSTRUCT_EXPORTS EXTERN __declspec(dllexport)
        // #pragma message "Message: 1 "
    #else
        #define UDSTRUCT_EXPORTS EXTERN __attribute__((__visibility__("default")))
// #pragma message "Message: 2 "
    #endif
#else
    #if defined(_WIN32) || defined(_WIN64)
        #define UDSTRUCT_EXPORTS EXTERN // __declspec(dllimport)
// #pragma message "Message: 3 "
    #else
        #define UDSTRUCT_EXPORTS EXTERN
// #pragma message "Message: 4 "
    #endif
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
    #include <io.h>
#endif

#include "hdf5.h"

typedef struct Input Input;

struct Input {
    float *dataMZ;
    float *dataInt;
    float *testmasses;
    int *nztab;
    float *mtab;
    char *barr;
    float *fwhmlist;
};

typedef struct Decon Decon;

struct Decon {
    float *fitdat;
    float *baseline;
    float *noise;
    float *massgrid;
    float *massaxis;
    float *massaxisval;
    float *blur;
    float *newblur;
    float *peakx;
    float *peaky;
    float *dscores;
    int *starttab;
    int *endtab;
    float *mzdist;
    float *rmzdist;
    float error;
    float rsquared;
    int iterations;
    float uniscore;
    float conv;
    float threshold;
    int mlen;
    int plen;
    int scanindex;
    int maxlength;
};


typedef struct Config Config;

struct Config {
    char infile[500];
    char outfile[500];
    int numit;
    int numz;
    int endz;
    int startz;
    float zsig;
    float psig;
    float beta;
    float mzsig;
    float msig;
    float molig;
    float massub;
    float masslb;
    int psfun;
    int zpsfun;
    float psmzthresh;
    float mtabsig;
    char mfile[500];
    char manualfile[500];
    int mflag;
    float massbins;
    int limitflag;
    float psthresh;
    int speedyflag;
    int linflag;
    int aggressiveflag;
    float adductmass;
    int rawflag;
    float nativezub;
    float nativezlb;
    int poolflag;
    int manualflag;
    float intthresh;
    float peakshapeinflate;
    int fixedmassaxis;
    int isotopemode;
    int filetype;
    int imflag;
    //IM Parameters
    float dtsig;
    float csig;
    float ccsub;
    float ccslb;
    float ccsbins;
    float temp;
    float press;
    float volt;
    float tcal1;
    float tcal2;
    float tcal3;
    float tcal4;
    int twaveflag;
    float hmass;
    float to;
    float len;
    float edc;
    float nativeccsub;
    float nativeccslb;
    int baselineflag;
    int noiseflag;
    int zout;
    int metamode;
    float minmz;
    float maxmz;
    float mzres;
    float mzbins;
    float bsub;
    float datareduction;
    float peakwin;
    float peakthresh;
    float exwindow;
    int exchoice;
    int exchoicez;
    float exthresh;
    int exnorm;
    int exnormz;
    int peaknorm;
    int orbimode;
    int datanorm;
    //Experimental Parameters
    int filterwidth;
    float zerolog;
    int lengthmz;
    int mfilelen;
    int isolength;
    // DoubleDec Parameters
    int doubledec;
    char kernel[500];
    hid_t file_id;
    char dataset[1024];
    // Other
    int silent;
    int cdmsflag;
    int variablepw;
    float minratio;
};

UDSTRUCT_EXPORTS Input SetupInputs();
UDSTRUCT_EXPORTS void FreeInputs(Input inp);
UDSTRUCT_EXPORTS Decon SetupDecon();
UDSTRUCT_EXPORTS void FreeDecon(Decon decon);
UDSTRUCT_EXPORTS void SetDefaultConfig(Config *config);
UDSTRUCT_EXPORTS void PostImport(Config *config);
UDSTRUCT_EXPORTS Config LoadConfig(Config config, const char *filename);
UDSTRUCT_EXPORTS void PrintHelp();

#endif //UDSTRUCT_H