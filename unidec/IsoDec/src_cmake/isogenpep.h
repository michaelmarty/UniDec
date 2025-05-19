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
ISOGENPEP_EXPORTS void mass_to_formula_averaging(float mass, int* formula);
ISOGENPEP_EXPORTS float fft_pep_list_to_dist(const int isolist[5], int length, float *isodist);
ISOGENPEP_EXPORTS float isodist_from_peptide_sequence(const char *sequence, float *isodist, int isolen);
ISOGENPEP_EXPORTS float isodist_from_nt_sequence(const char *sequence, float *isodist, int isolen);
ISOGENPEP_EXPORTS float fft_pep_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset);
ISOGENPEP_EXPORTS float isomike(float mass, float *isodist, int isolen, int offset);
ISOGENPEP_EXPORTS float isogen_fancy(float mass, float *isodist, int isolen, int offset, char* type);
ISOGENPEP_EXPORTS float isodist_from_averagine_mass(const float mass, float* isodist, const int isolen, const int offset);
ISOGENPEP_EXPORTS float fft_pep_seq_to_dist(const char* sequence, float* isodist, const int isolen);

#endif
