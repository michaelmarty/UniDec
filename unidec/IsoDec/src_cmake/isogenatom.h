#ifndef ISOGENATOM_H
#define ISOGENATOM_H

#ifdef __cplusplus
    #define EXTERN extern "C"
#else
    #define EXTERN
#endif

#ifdef ISOGEN_BUILD_DLL
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENATOM_EXPORTS EXTERN __declspec(dllexport)
    #else
        #define ISOGENATOM_EXPORTS EXTERN __attribute__((__visibility__("default")))
    #endif
#else
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENATOM_EXPORTS EXTERN __declspec(dllimport)
    #else
        #define ISOGENATOM_EXPORTS EXTERN
    #endif
#endif

ISOGENATOM_EXPORTS float isodist_from_formula_vector(const int* formula, const int formulalen, float* isodist, const int isolen);
ISOGENATOM_EXPORTS float isodist_from_formula(const char *formula, float *isodist, int isolen);
ISOGENATOM_EXPORTS void formula_to_vector(const char* formula, int* formulalist);
ISOGENATOM_EXPORTS float* get_formula_dist(char* sequence);
ISOGENATOM_EXPORTS void fft_atom_sequence_to_dist(const char* sequence, float* isodist, int isolen);

#endif
