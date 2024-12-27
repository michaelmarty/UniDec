#ifndef ISODECLIB_LIBRARY_H
#define ISODECLIB_LIBRARY_H

#ifndef ISODEC_EXPORTS_H
#define ISODEC_EXPORTS_H

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

#ifdef ISODECLIB_STATIC_DEFINE
#  define ISODEC_EXPORTS
#else
#  ifndef ISODEC_EXPORTS
#    ifdef isodeclib_EXPORTS
# if OS_Windows
		/* We are building this library */
#      define ISODEC_EXPORTS EXTERN __declspec(dllexport)
#    else
		/* We are building this library */
# 	define ISODEC_EXPORTS EXTERN __attribute__((__visibility__("default")))
#    endif
#    else
# if OS_Windows
        /* We are using this library */
#      define ISODEC_EXPORTS EXTERN __declspec(dllimport)
#    else
		/* We are using this library */
# 	define ISODEC_EXPORTS EXTERN
#   endif
#    endif
#  endif

#endif

#endif /* ISODEC_EXPORTS_H */

//  "isodec_exports.h"


// Structure for matched peak
struct MatchedPeak
{
	float mz;
	int z;
	float monoiso;
	float peakmass;
	float avgmass;
	float area;
	float peakint;
	int matchedindsiso[64];
	int matchedindsexp[64];
	float isomz[64];
	float isodist[64];
	float isomass[64];
	float monoisos[16];
	int startindex;
	int endindex;
	float score;
	int realisolength;
};

// Structure for the config object. Mostly neural net parameters. The settings structure has the parameters for the peak detection and isotope distribution.
struct IsoConfig
{
	int verbose; // Verbose output
	int pres; // Precision of encoding matrix
	int maxz; // Maximum charge state
	int elen; // Encoding Length
	int l1; // Weights 1 Length
	int l2; // Bias 1 Length
	int l3; // Weights 2 Length
	int l4; // Bias 2 Length
	int dlen; // Data Length
};

// Structure for the settings object. Parameters for peak detection and isotope distribution checking.
struct IsoSettings
{
	int phaseres; // Precision of encoding matrix, USER
	int verbose; // Verbose output
	int peakwindow; // Peak Detection Window, USER
	float peakthresh; // Peak Detection Threshold, USER
	int minpeaks; // Minimum Peaks for an allowed peak
	float css_thresh; // Minimum cosine similarity score for isotope distribution, USER
	float matchtol; // Match Tolerance for peak detection in ppm, USER
	int maxshift; // Maximum shift allowed for isotope distribution, USER
	float mzwindow[2]; // MZ Window for isotope distribution, USER
	float plusoneintwindow[2]; // Plus One Intensity range. Will be used for charge state 1
	int knockdown_rounds; // Number of knockdown rounds, USER
	float min_score_diff; // Minimum score difference for isotope distribution to allow missed monoisotopic peaks
	float minareacovered; // Minimum area covered by isotope distribution. Use in or with css_thresh
	int isolength; // Isotope Distribution Length
	double mass_diff_c; // Mass difference between isotopes
	float adductmass; // Adduct Mass, USER FOR NEGATIVE MODE
	int minusoneaszero; // Use set the -1 isotope as 0 to help force better alignments
	float isotopethreshold; // Threshold for isotope distribution. Will remove relative intensities below this.
	float datathreshold; // Threshold for data. Will remove relative intensities below this relative to max intensity in each cluster, USER
	float zscore_threshold; //Ratio above which a secondary charge state prediction will be returned.
};

// Structure for the weights and biases of the neural network
struct Weights {
	float *w1;
	float *b1;
	float *w2;
	float *b2;
};

ISODEC_EXPORTS void run(char *filename, char *outfile, const char *weightfile);
ISODEC_EXPORTS int encode(const double* cmz, const float* cint, int n, float * emat, struct IsoConfig config, struct IsoSettings settings);
ISODEC_EXPORTS void predict_charge(const double* cmz, const float* cint, int n, const char* fname, int* charge);
ISODEC_EXPORTS int process_spectrum(const double* cmz, const float* cint, int n, const char* fname, struct MatchedPeak * matchedpeaks, struct IsoSettings settings);
ISODEC_EXPORTS int process_spectrum_default(const double* cmz, const float* cint, int n, const char* fname, struct MatchedPeak * matchedpeaks);

#endif //ISODECLIB_LIBRARY_H
