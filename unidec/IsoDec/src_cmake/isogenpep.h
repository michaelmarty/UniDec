#ifndef ISOGENPEP_H
#define ISOGENPEP_H

#ifdef __cplusplus
    #define EXTERN extern "C"
#else
    #define EXTERN
#endif

#ifdef ISOGEN_BUILD_DLL
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENPEP_EXPORTS EXTERN __declspec(dllexport)
    #else
        #define ISOGENPEP_EXPORTS EXTERN __attribute__((__visibility__("default")))
    #endif
#else
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENPEP_EXPORTS EXTERN __declspec(dllimport)
    #else
        #define ISOGENPEP_EXPORTS EXTERN
    #endif
#endif


// Function declarations
ISOGENPEP_EXPORTS float fft_pep_mass_to_dist(float mass, float* isodist, int isolen, int offset);
ISOGENPEP_EXPORTS float fft_pep_seq_to_dist(const char* sequence, float* isodist, int isolen, int offset);
ISOGENPEP_EXPORTS float pep_mass_to_dist_fitting(float mass, float *isodist, int isolen, int offset);
ISOGENPEP_EXPORTS float nn_pep_seq_to_dist(const char* seq, float* isodist, int offset);
ISOGENPEP_EXPORTS float nn_pep_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset);
#endif
