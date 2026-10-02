#ifndef ISODECLIB_LIBRARY_H
#define ISODECLIB_LIBRARY_H

#ifdef __cplusplus
	#define EXTERN extern "C"
#else
	#define EXTERN
#endif

// 🔥 Define ISODECLIB_EXPORTS ONLY inside isodeclib.h
#ifdef ISODEC_BUILD_DLL
	#if defined(_WIN32) || defined(_WIN64)
		#define ISODECLIB_EXPORTS EXTERN __declspec(dllexport)
	#else
		#define ISODECLIB_EXPORTS EXTERN __attribute__((__visibility__("default")))
	#endif
#else
	#if defined(_WIN32) || defined(_WIN64)
		#define ISODECLIB_EXPORTS EXTERN __declspec(dllimport)
	#else
		#define ISODECLIB_EXPORTS EXTERN
	#endif
#endif

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

#define ISODEC_FRAGMENT_ISOTOPES 128
struct FragmentHit {
    int fragment_index;
    int charge;
    int left;
    int right;
    int isotope_count;
    int match_count;
    double score;
    double scale;
    int isotope_indexes[ISODEC_FRAGMENT_ISOTOPES];
    int centroid_indexes[ISODEC_FRAGMENT_ISOTOPES];
};

// The caller owns the output through free_fragment_hits(). Returns 0 on success.
ISODECLIB_EXPORTS int match_fragment_batch(
    const double* mz, const double* intensity, int spectrum_count,
    const double* masses, const float* envelopes, int fragment_count, int isolen,
    int max_charge, double adduct, double isotope_threshold, double ppm_tolerance,
    int min_peaks, double cosine_threshold, double min_area, int minus_one_as_zero,
    struct FragmentHit** output, int* output_count);
ISODECLIB_EXPORTS void free_fragment_hits(struct FragmentHit* hits);

// Match theoretical m/z arrays against exported centroid windows using the
// mass-grouping selection rule. Offsets are exclusive; output indexes are
// relative to each window. Caller owns all input and output buffers.
// Returns 0 on success and -1 for invalid input or insufficient output space.
ISODECLIB_EXPORTS int match_group_centroids_batch(
    const double* mz, const double* intensity, int spectrum_count,
    const int* windows, const double* isotope_mz, const int* isotope_offsets,
    int peak_count, double ppm_tolerance,
    int* match_counts, int* match_indexes, int match_capacity);

struct MassGroupPeak {
    double monoiso, mz, peakint, matchedintensity, avgmass, rt;
    int charge, scan;
    const double* monoisos;
    int monoisos_count;
    const double* massdist;       // consecutive mass, intensity pairs
    int massdist_count;
    const double* decon_centroids; // consecutive mass, intensity pairs
    int centroid_count;
    int float32_intensity, float32_massdist;
};

struct MassGroupResult {
    double monoiso, lookup_mass, apexintensity, totalintensity, avgmass;
    double minrt, maxrt, apexrt;
    int minscan, maxscan, apexscan, seed_index, totalpeaks;
    double* monoisos;
    int monoisos_count;
    double* massdist;
    int massdist_count;
    double* decon_centroids;
    int centroid_count;
};

#define ISODEC_MASS_GROUP_ABI_VERSION 1
ISODECLIB_EXPORTS int mass_group_abi_version(void);
// Results, nested arrays, and hit-to-group IDs are owned by the caller via
// free_mass_groups(). order: 0 = input order, 1 = descending matched intensity.
// Returns 0 on success, -1 for invalid input, count overflow, or allocation failure.
// With valid output pointers, failures reset them to NULL/0. Input counts must
// describe allocated buffers; pair counts are limited to INT_MAX/2. Python
// validates shapes before this boundary. Mass/intensity arithmetic uses double
// when the legacy float32 flags are zero (the current Python contract).
ISODECLIB_EXPORTS int group_mass_peaks_batch(
    const struct MassGroupPeak* peaks, int peak_count,
    double ppm_tolerance, int maxshift, double mass_diff_c,
    double cosine_threshold, int scan_tolerance, int order,
    struct MassGroupResult** groups, int* group_count, int** group_ids);
ISODECLIB_EXPORTS void free_mass_groups(
    struct MassGroupResult* groups, int group_count, int* group_ids);

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

ISODECLIB_EXPORTS void run(char *filename, char *outfile, const char *weightfile, const char* type);
ISODECLIB_EXPORTS int encode(const double* cmz, const float* cint, int n, float * emat, struct IsoConfig config, struct IsoSettings settings);
ISODECLIB_EXPORTS int predict_charge(const double* cmz, const float* cint, int n, const char* fname);
ISODECLIB_EXPORTS int process_spectrum(const double* cmz, const float* cint, int n, const char* fname, struct MatchedPeak * matchedpeaks, struct IsoSettings settings, const char* type);
ISODECLIB_EXPORTS int process_spectrum_default(const double* cmz, const float* cint, int n, const char* fname, struct MatchedPeak * matchedpeaks, const char* type);
ISODECLIB_EXPORTS struct IsoSettings DefaultSettings();

#endif //ISODECLIB_LIBRARY_H
