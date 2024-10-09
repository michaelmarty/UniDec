#include "isodec_exports.h"

#ifndef ISODECLIB_LIBRARY_H
#define ISODECLIB_LIBRARY_H

extern ISODEC_EXPORTS void run(void);
extern ISODEC_EXPORTS int encode(const double* cmz, const float* cint, const int n, float * emat, const struct IsoConfig config, const struct IsoSettings settings);
extern ISODEC_EXPORTS void predict_charge(const double* cmz, const float* cint, const int n, const char* fname, int* charge);
extern ISODEC_EXPORTS int process_spectrum(const double* cmz, const float* cint, int n, const char* fname, struct MatchedPeak * matchedpeaks, struct IsoSettings settings);
extern ISODEC_EXPORTS int process_spectrum_default(const double* cmz, const float* cint, const int n, const char* fname, struct MatchedPeak * matchedpeaks);

#endif //ISODECLIB_LIBRARY_H
