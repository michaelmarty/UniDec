//
// Created by Michael Marty on 9/23/2025.
//

#ifndef UDINTERFACE_H
#define UDINTERFACE_H

#include "udstruct.h"
#include "udmain.h"
#include <vector>

// Create vector class to mirror the Input C struct

class InputVec {
public:
    std::vector<float> dataMZ;
    std::vector<float> dataInt;
    std::vector<float> testmasses;
    std::vector<int> nztab;
    std::vector<float> mtab;
    std::vector<bool> barr;
};

// Create Vector class to mirror the Decon C struct
class DeconVec {
public:
    std::vector<float> fitdat;
    std::vector<float> baseline;
    std::vector<float> noise;
    std::vector<float> massgrid;
    std::vector<float> massaxis;
    std::vector<float> massaxisval;
    std::vector<float> blur;
    std::vector<float> newblur;
    std::vector<float> peakx;
    std::vector<float> peaky;
    std::vector<float> dscores;
    std::vector<int> starttab;
    std::vector<int> endtab;
    std::vector<float> mzdist;
    std::vector<float> rmzdist;
    float error = 0.0f;
    float rsquared = 0.0f;
    int iterations = 0;
    float uniscore = 0.0f;
    float conv = 0.0f;
    float threshold = 0.0f;
    int mlen = 0;
    int plen = 0;
    int scanindex = 0;
};

// Expose this function to DLL DeconVec RunDeconvolution(const Config& config, const InputVec& inputVec, int silent, int verbose)
DeconVec RunDeconvolution(Config *config, InputVec *inputVec, int verbose);

#endif //UDINTERFACE_H
