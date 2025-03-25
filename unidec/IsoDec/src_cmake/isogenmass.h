#ifndef ISOGENMASS_H
#define ISOGENMASS_H

#ifdef __cplusplus
    #define EXTERN extern "C"
#else
    #define EXTERN
#endif

// 🔥 Define ISOGENMASS_EXPORTS ONLY inside isogenmass.h
#ifdef ISOGEN_BUILD_DLL
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENMASS_EXPORTS EXTERN __declspec(dllexport)
    #else
        #define ISOGENMASS_EXPORTS EXTERN __attribute__((__visibility__("default")))
    #endif
#else
    #if defined(_WIN32) || defined(_WIN64)
        #define ISOGENMASS_EXPORTS EXTERN __declspec(dllimport)
    #else
        #define ISOGENMASS_EXPORTS EXTERN
    #endif
#endif


ISOGENMASS_EXPORTS void isogenmass(float mass, float *isodist);
ISOGENMASS_EXPORTS float isomike(float mass, float *isodist, int isolen, int offset);
ISOGENMASS_EXPORTS float isogenmass_fancy(float mass, float *isodist, int isolen, int offset);

// ✅ Structure Declaration
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

#endif // ISOGENMASS_H
