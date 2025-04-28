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
ISOGENRNA_EXPORTS float* rna_mass_to_list(float initialMass);
ISOGENRNA_EXPORTS float fft_rna_list_to_dist(const float isolist[5], const int length, float *isodist);
ISOGENRNA_EXPORTS float* fft_rna_seq_to_dist(char* sequence);
ISOGENRNA_EXPORTS float* rnaToVector(const char* rna);
ISOGENRNA_EXPORTS float rnaVectorToMass(float* rnaVector);
ISOGENRNA_EXPORTS float fft_rna_mass_to_dist(float mass, float* isovals, int isolen, int offset);

#endif
