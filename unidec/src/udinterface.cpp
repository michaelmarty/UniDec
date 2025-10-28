//
// Created by Michael Marty on 9/23/2025.
//

#include "udstruct.h"
#include "udinterface.hpp"
#include <Rcpp.h>

// Function to convert CPP InputVec to C Input struct
void ConvertToCInput(const InputVec * inputVec, Input * inputC) {
    inputC->dataMZ = const_cast<float*>(inputVec->dataMZ.data());
    inputC->dataInt = const_cast<float*>(inputVec->dataInt.data());
}

// Function to convert C Decon struct to CPP DeconVec
DeconVec ConvertToCPPDecon(const Decon * deconC) {
    DeconVec deconVec;
    // deconVec.fitdat.assign(deconC->fitdat, deconC->fitdat + deconC->mlen);
    // deconVec.baseline.assign(deconC->baseline, deconC->baseline + deconC->mlen);
    // deconVec.noise.assign(deconC->noise, deconC->noise + deconC->mlen);
    // deconVec.massgrid.assign(deconC->massgrid, deconC->massgrid + deconC->mlen);
    deconVec.massaxis.assign(deconC->massaxis, deconC->massaxis + deconC->mlen);
    deconVec.massaxisval.assign(deconC->massaxisval, deconC->massaxisval + deconC->mlen);
    // deconVec.blur.assign(deconC->blur, deconC->blur + deconC->mlen);
    // deconVec.newblur.assign(deconC->newblur, deconC->newblur + deconC->mlen);
    // deconVec.peakx.assign(deconC->peakx, deconC->peakx + deconC->plen);
    // deconVec.peaky.assign(deconC->peaky, deconC->peaky + deconC->plen);
    // deconVec.dscores.assign(deconC->dscores, deconC->dscores + deconC->plen);
    // deconVec.starttab.assign(deconC->starttab, deconC->starttab + deconC->plen);
    // deconVec.endtab.assign(deconC->endtab, deconC->endtab + deconC->plen);
    // deconVec.mzdist.assign(deconC->mzdist, deconC->mzdist + deconC->mlen);
    // deconVec.rmzdist.assign(deconC->rmzdist, deconC->rmzdist + deconC->mlen);

    // deconVec.error = deconC->error;
    // deconVec.rsquared = deconC->rsquared;
    // deconVec.iterations = deconC->iterations;
    // deconVec.uniscore = deconC->uniscore;
    // deconVec.conv = deconC->conv;
    // deconVec.threshold = deconC->threshold;
    // deconVec.mlen = deconC->mlen;
    // deconVec.plen = deconC->plen;
    // deconVec.scanindex = deconC->scanindex;

    return deconVec;
}

extern "C" int run_unidec_core(Config config, Input inp, Decon *decon, const int verbose, const int autotune);

// Main interface function to run deconvolution
DeconVec RunDeconvolution(Config * config, InputVec * inputVec, int verbose) {
    Input inputC = SetupInputs();
    ConvertToCInput(inputVec, &inputC);
    PostImport(config);
    Decon deconC;
    // Print startz from config for debugging
    SetupZtab(*config, &inputC);
    // Print lengthmz and ztab[1] for debugging
    // Rprintf( "Lengthmz: %d, Ztab[1]: %d\n", config->lengthmz, inputC.nztab[1]);
    run_unidec_core(*config, inputC, &deconC, verbose, 0);
    DeconVec deconVec = ConvertToCPPDecon(&deconC);
    return deconVec;
}

