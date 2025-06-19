#ifndef ISOGENRNA_H
#define ISOGENRNA_H

#ifdef __cplusplus
    #define EXTERN extern "C"
#else
    #define EXTERN
#endif

#ifdef ISOGEN_BUILD_DLL
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENRNA_EXPORTS EXTERN __declspec(dllexport)
    #else
        #define ISOGENRNA_EXPORTS EXTERN __attribute__((__visibility__("default")))
    #endif
#else
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENRNA_EXPORTS EXTERN __declspec(dllimport)
    #else
        #define ISOGENRNA_EXPORTS EXTERN
    #endif
#endif

// Declare functions implemented in isogenrna.c
ISOGENRNA_EXPORTS float fft_rna_mass_to_dist(float mass, float* isodist, int isolen, int offset);
ISOGENRNA_EXPORTS float fft_rna_seq_to_dist(const char* sequence, float* isodist, int isolen, int offset);
ISOGENRNA_EXPORTS float nn_rna_seq_to_dist(const char* seq, float* isodist, int offset);
ISOGENRNA_EXPORTS float nn_rna_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset);
#endif
