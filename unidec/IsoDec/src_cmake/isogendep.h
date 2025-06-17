#ifndef ISOGENDEP_LIBRARY_H
#define ISOGENDEP_LIBRARY_H

#ifdef __cplusplus
    #define EXTERN extern "C"
#else
    #define EXTERN
#endif

// 🔥 Define ISOGENDEP_EXPORTS ONLY inside isogendep.h
#ifdef ISOGEN_BUILD_DLL
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENDEP_EXPORTS EXTERN __declspec(dllexport)
    #else
        #define ISOGENDEP_EXPORTS EXTERN __attribute__((__visibility__("default")))
    #endif
#else
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENDEP_EXPORTS EXTERN __declspec(dllimport)
    #else
        #define ISOGENDEP_EXPORTS EXTERN
    #endif
#endif


// Ensure other headers are found at compile time
#include "fftw3.h"
#include "isogenrna.h"

struct IsoGenWeights {
    float *w1;
    float *b1;
    float *w2;
    float *b2;
    float *w3;
    float *b3;
    int vl1;
    int vl2;
    int vl3;
    int vl4;
    int ml1;
    int ml2;
    int ml3;
    int tot;
};

// Function Declarations
ISOGENDEP_EXPORTS void setup_ft(const int number, fftw_complex* outft, const int length, const int ftlen);
ISOGENDEP_EXPORTS void convolve_all(const int isolist[11], fftw_complex* allft, const fftw_complex* cft, const fftw_complex* hft,
                  const fftw_complex* nft, const fftw_complex* oft, const fftw_complex* sft, const fftw_complex* feft,
                  const fftw_complex* caft, const fftw_complex* kft, const fftw_complex* nift, const fftw_complex* znft,
                  const fftw_complex* mgft, const int length);
ISOGENDEP_EXPORTS void complex_power(const fftw_complex in, const int power, double* outreal, double* outimag);
ISOGENDEP_EXPORTS void complex_multiplication(const double ar, const double ai, const double br,
                                              const double bi, double* oreal, double* oimag);
ISOGENDEP_EXPORTS double normalize_isodist(double* buffer, int length);
ISOGENDEP_EXPORTS void FreeIsogenWeights(const struct IsoGenWeights weights);
ISOGENDEP_EXPORTS struct IsoGenWeights LoadWeights(const struct IsoGenWeights weights, const unsigned char *model_weights);
ISOGENDEP_EXPORTS struct IsoGenWeights SetupWeights(const int veclen, const int outlen);
ISOGENDEP_EXPORTS void neural_net(const float *vector, float *isodist, const struct IsoGenWeights weights);
ISOGENDEP_EXPORTS float softsign(const float x);
ISOGENDEP_EXPORTS float sigmoid(const float x);
ISOGENDEP_EXPORTS void mass_to_vector(const float mass, float *vector);
ISOGENDEP_EXPORTS void isogenmass_nn(const float mass, float *isodist, int isolen, const char* analyte);
ISOGENDEP_EXPORTS void isogen_atom(const char* sequence, float* isodist, int isolen);
ISOGENDEP_EXPORTS float fft_list_to_dist(const int isolist[11], const int length, float* isodist);

extern const double fftarray[3597][2];
extern const int arraylen;
extern const char* element_names[];
extern const int isotope_numbers[];
extern double isotope_masses[];
extern double isotope_abundances[];
extern char* amino_acid_names[];
extern int amino_acid_vectors[][5];
extern const int dna_vectors[][5];
extern int num_simp_elements;
extern const int numelements;
extern char* rna_names[];
extern char* dna_names[];
extern const char pepOrder[];
extern const char *elements[];


#endif // ISOGENDEP_LIBRARY_H
