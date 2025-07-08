//
// Created by mm96978 on 7/7/2025.
//

#ifndef UDSTRUCT_H
#define UDSTRUCT_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <io.h>
#include "hdf5.h"

typedef struct Input Input;

struct Input {
    float *dataMZ;
    float *dataInt;
    float *testmasses;
    int *nztab;
    float *mtab;
    char *barr;
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
    float twaveflag;
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
    int mzbins;
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
    int silent;
    int cdmsflag;
};

Input SetupInputs();
void FreeInputs(Input inp);
Decon SetupDecon();
void FreeDecon(Decon decon);
Config SetDefaultConfig();
Config PostImport(Config config);
Config LoadConfig(Config config, const char *filename);
void PrintHelp();

#endif //UDSTRUCT_H