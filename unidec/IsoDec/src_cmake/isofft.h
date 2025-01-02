#ifndef ISOFFT_LIBRARY_H
#define ISOFFT_LIBRARY_H

#ifndef ISOFFT_EXPORTS_H
#define ISOFFT_EXPORTS_H

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

#ifdef ISOFFTLIB_STATIC_DEFINE
#  define ISOFFT_EXPORTS
#else
#  ifndef ISOFFT_EXPORTS
#    ifdef ISOFFTlib_EXPORTS
# if OS_Windows
        /* We are building this library */
#      define ISOFFT_EXPORTS EXTERN __declspec(dllexport)
#    else
        /* We are building this library */
# 	define ISOFFT_EXPORTS EXTERN __attribute__((__visibility__("default")))
#    endif
#    else
# if OS_Windows
        /* We are using this library */
#      define ISOFFT_EXPORTS EXTERN __declspec(dllimport)
#    else
        /* We are using this library */
# 	define ISOFFT_EXPORTS EXTERN
#   endif
#    endif
#  endif

#endif

#endif


ISOFFT_EXPORTS void test_averagine(float mass);
ISOFFT_EXPORTS float isodist_from_averagine_mass(float mass, float * isodist, int isolen, int offset);
ISOFFT_EXPORTS float fftcalc_isodist(const int isolist[5], int length, float *isodist);
ISOFFT_EXPORTS void test_formula(const char *formula);
ISOFFT_EXPORTS float isodist_from_formula(const char *formula, float *isodist, int isolen);
ISOFFT_EXPORTS void test_peptide_sequence(const char *sequence);
ISOFFT_EXPORTS float isodist_from_peptide_sequence(const char *sequence, float *isodist, int isolen);
ISOFFT_EXPORTS void test_nt_sequence(const char *sequence);
ISOFFT_EXPORTS float isodist_from_nt_sequence(const char *sequence, float *isodist, int isolen);

#endif //ISOFFT_LIBRARY_H



//// "Returns" discrete Fourier transform of input array. Unused.
//// Use fftw library instead
//void discretefouriertransform(double* input, double** output, int length) {
//	double pi = 3.1415926535;
//	for (int i = 0; i < length; i++) {
//		output[i][0] = 0.0; // Real term
//		output[i][1] = 0.0; // Imaginary term
//		for (int j = 0; j < length; j++) {
//			double inner = 2.0 * pi / length * i * j;
//			output[i][0] += input[j] * cos(inner);
//			output[i][1] += input[j] * (-1.0) * sin(inner);
//		}
//	}
//}
//
//// "Returns" discrete inverse Fourier transform of input. Unused.
//// Use fftw library instead
//void inversefouriertransform(double** input, double* output, int length) {
//	double pi = 3.1415926535;
//	for (int i = 0; i < length; i++) {
//		output[i] = 0.0;	// Real term
//		double imag = 0.0;	// Imaginary term
//		for (int j = 0; j < length; j++) {	// (ac - bd) + i(ad + bc)
//			double inner = 2.0 * pi / length * i * j;
//			double c = cos(inner);
//			double d = sin(inner);
//			output[i] += (input[j][0] * c) - (input[j][1] * d);
//			imag += (input[j][0] * d) + (input[j][1] * c);
//		}
//		output[i] = output[i] / ((double)length);
//		imag = imag / ((double)length);
//		printf("Imaginary term is %.2f", imag);
//	}
//}

// Gives convolution of functions a and b. Unused.
// Use cconv2fast instead
/*
void cconv2(double* a, double* b, double* c, int length) {
	double** A = (double **) malloc(length * sizeof(double*));	// (a + bi)
	double** B = (double**) malloc(length * sizeof(double*));	// (c + di)
	double** C = (double**) malloc(length * sizeof(double*));
	for (int i = 0; i < length; i++) {
		A[i] = (double*) calloc(2, sizeof(double));
		B[i] = (double*) calloc(2, sizeof(double));
		C[i] = (double*) calloc(2, sizeof(double));
	}

	discretefouriertransform(a, A, length);
	discretefouriertransform(b, B, length);
	// A * B = (ac - bd) + i(ad + bc)
	for (int j = 0; j < length; j++) {
		C[j][0] = (A[j][0] * B[j][0]) - (A[j][1] * B[j][1]);
		C[j][1] = (A[j][0] * B[j][1]) + (A[j][1] * B[j][0]);
	}
	inversefouriertransform(C, c, length);
	for (int k = 0; k < length; k++) {
		if (c[k] < 0) {
			c[k] = c[k] * -1.0;
		}
		free(A[k]);
		free(B[k]);
		free(C[k]);
	}
	free(A);
	free(B);
	free(C);
}
*/