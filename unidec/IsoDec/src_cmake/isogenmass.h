//
// Created by marty on 10/17/2024.
//
#ifndef ISOGEN_H
#define ISOGEN_H

# if defined(__cplusplus)
#  define EXTERN extern "C"
# else
#  define EXTERN
# endif

# ifndef OS_Windows
    # define OS_Windows 0
    # if defined(_WIN32) || defined(_WIN64)
        // Set OS_Windows to 1
        #	undef OS_Windows
        #	define OS_Windows 1
    # endif
# endif

#ifdef ISOGEN_STATIC_DEFINE
#  define ISOGEN_EXPORTS
#else
#  ifndef ISOGEN_EXPORTS
#    ifdef isogen_EXPORTS
# if OS_Windows
        /* We are building this library */
#      define ISOGEN_EXPORTS EXTERN __declspec(dllexport)
#    else
        /* We are building this library */
# 	define ISOGEN_EXPORTS EXTERN __attribute__((__visibility__("default")))
#    endif
#    else
# if OS_Windows
        /* We are using this library */
#      define ISOGEN_EXPORTS EXTERN __declspec(dllimport)
#    else
        /* We are using this library */
# 	define ISOGEN_EXPORTS EXTERN
#   endif
#    endif
#  endif

#endif

#endif //ISOGEN_H

#ifndef ISOGENMASS_H
#define ISOGENMASS_H

ISOGEN_EXPORTS void isogenmass(float mass, float *isodist);
ISOGEN_EXPORTS float isomike(float mass, float *isodist, int isolen, int offset);
ISOGEN_EXPORTS float isogenmass_fancy(float mass, float *isodist, int isolen, int offset);

// Structure for the weights and biases of the neural network
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

#endif //ISOGENMASS_H



