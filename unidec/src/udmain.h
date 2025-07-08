//
// Created by Michael Marty on 7/8/2025.
//

#ifndef UDMAIN_H
#define UDMAIN_H

#ifdef __cplusplus
    #define EXTERN extern "C"
#else
    #define EXTERN
#endif

// 🔥 Define UNIDECLIB_EXPORTS ONLY inside isodeclib.h
#ifdef UNIDEC_BUILD_DLL
    #if defined(_WIN32) || defined(_WIN64)
        #define UNIDECLIB_EXPORTS EXTERN __declspec(dllexport)
        // #pragma message "Message: 1 "
    #else
        #define UNIDECLIB_EXPORTS EXTERN __attribute__((__visibility__("default")))
// #pragma message "Message: 2 "
    #endif
#else
    #if defined(_WIN32) || defined(_WIN64)
        #define UNIDECLIB_EXPORTS EXTERN // __declspec(dllimport)
// #pragma message "Message: 3 "
    #else
        #define UNIDECLIB_EXPORTS EXTERN
// #pragma message "Message: 4 "
    #endif
#endif



#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include "hdf5.h"
#include "udcore.h"
#include "udtools.h"
#include "udio.h"
#include "h5io.h"
#include "UD_score.h"
#include "udstruct.h"

UNIDECLIB_EXPORTS Decon ExitToBlank(const Config config, Decon decon);
UNIDECLIB_EXPORTS Decon MainDeconvolution(const Config config, const Input inp, const int silent, const int verbose);
UNIDECLIB_EXPORTS void RunAutotune(Config *config, const Input *inp, Decon *decon);
UNIDECLIB_EXPORTS int run_unidec_core(Config config, Input inp, Decon *decon, const int verbose, const int autotune);
UNIDECLIB_EXPORTS int run_unidec(int argc, char *argv[], Config config);

#endif //UDMAIN_H
