//
// Created by Michael Marty on 9/23/2025.
//

#include "udinterface.hpp"

// Function to convert CPP InputVec to C Input struct
Input ConvertToCInput(const InputVec& inputVec) {
    Input inputC;
    inputC.dataMZ = const_cast<float*>(inputVec.dataMZ.data());
    inputC.dataInt = const_cast<float*>(inputVec.dataInt.data());
    inputC.testmasses = const_cast<float*>(inputVec.testmasses.data());
    inputC.nztab = const_cast<int*>(inputVec.nztab.data());
    inputC.mtab = const_cast<float*>(inputVec.mtab.data());
    // Allocate a mutable buffer and copy the string data
    std::vector<char> barrBuffer(inputVec.barr.begin(), inputVec.barr.end());
    barrBuffer.push_back('\0'); // Ensure null-termination
    inputC.barr = barrBuffer.data();
    return inputC;
}

// Function to convert C Decon struct to CPP DeconVec
DeconVec ConvertToCPPDecon(const Decon& deconC) {
    DeconVec deconVec;
    deconVec.fitdat.assign(deconC.fitdat, deconC.fitdat + deconC.mlen);
    deconVec.baseline.assign(deconC.baseline, deconC.baseline + deconC.mlen);
    deconVec.noise.assign(deconC.noise, deconC.noise + deconC.mlen);
    deconVec.massgrid.assign(deconC.massgrid, deconC.massgrid + deconC.mlen);
    deconVec.massaxis.assign(deconC.massaxis, deconC.massaxis + deconC.mlen);
    deconVec.massaxisval.assign(deconC.massaxisval, deconC.massaxisval + deconC.mlen);
    deconVec.blur.assign(deconC.blur, deconC.blur + deconC.mlen);
    deconVec.newblur.assign(deconC.newblur, deconC.newblur + deconC.mlen);
    deconVec.peakx.assign(deconC.peakx, deconC.peakx + deconC.plen);
    deconVec.peaky.assign(deconC.peaky, deconC.peaky + deconC.plen);
    deconVec.dscores.assign(deconC.dscores, deconC.dscores + deconC.plen);
    deconVec.starttab.assign(deconC.starttab, deconC.starttab + deconC.plen);
    deconVec.endtab.assign(deconC.endtab, deconC.endtab + deconC.plen);
    deconVec.mzdist.assign(deconC.mzdist, deconC.mzdist + deconC.mlen);
    deconVec.rmzdist.assign(deconC.rmzdist, deconC.rmzdist + deconC.mlen);

    deconVec.error = deconC.error;
    deconVec.rsquared = deconC.rsquared;
    deconVec.iterations = deconC.iterations;
    deconVec.uniscore = deconC.uniscore;
    deconVec.conv = deconC.conv;
    deconVec.threshold = deconC.threshold;
    deconVec.mlen = deconC.mlen;
    deconVec.plen = deconC.plen;
    deconVec.scanindex = deconC.scanindex;

    return deconVec;
}

extern "C" int run_unidec_core(Config config, Input inp, Decon *decon, const int verbose, const int autotune);

// Main interface function to run deconvolution
DeconVec RunDeconvolution(const Config& config, const InputVec& inputVec, int silent, int verbose) {
    // Convert InputVec to C Input struct
    Input inputC = ConvertToCInput(inputVec);
    Decon deconC;
    // UNIDECLIB_EXPORTS int run_unidec_core(Config config, Input inp, Decon *decon, const int verbose, const int autotune);
    run_unidec_core(config, inputC, &deconC, verbose, 0);
    // Convert C Decon struct to DeconVec
    DeconVec deconVec = ConvertToCPPDecon(deconC);
    return deconVec;
}

