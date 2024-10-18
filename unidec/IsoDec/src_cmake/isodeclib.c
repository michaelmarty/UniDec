#define CRT_SECURE_NO_WARNINGS
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "isodeclib.h"
#include "phase_model_8.h"
#include "phase_model_4.h"
#include "isogenmass.h"
// #include <omp.h>

// Linux code to convert binary file to header file
// xxd -i phase_model_8.bin > phase_model_8.h

#if defined(__AVX__) || defined(__AVX2__) || defined(__AVX512F__)
// #include "mkl.h"
// void apply_relu(float * result, const int n) {
//     for (int i = 0; i < n; i++) {
//         if (result[i] < 0) { result[i] = 0; }
//     }
// }
//
// void matrix_vector_multiply(const float* matrix, const float* vector, const float* bias, float* result, const int N, const int M, const bool relu) {
//     // memcopy bias to result
//     memcpy(result, bias, N * sizeof(float));
//     cblas_sgemv (CblasRowMajor, CblasNoTrans, (MKL_INT) N, (MKL_INT) M, 1.0f, matrix, (MKL_INT) M, vector, (MKL_INT) 1, 1.0f, result, (MKL_INT) 1);
//     // relu
//     if (relu) {
//         apply_relu(result, N);
//     }
// }

#include <immintrin.h>

void matrix_vector_multiply(const float *matrix, const float *vector, const float *bias, float *result, const int N,
                            const int M, const bool relu) {
    for (int i = 0; i < N; i++) {
        float val = bias[i];
        __m256 sum = _mm256_setzero_ps(); // Initialize sum to zero
        int jend = 0;
        for (int j = 0; j <= M - 8; j += 8) {
            // Load 4 elements from the matrix and vector
            const __m256 mat = _mm256_loadu_ps(&matrix[i * M + j]);
            const __m256 vec = _mm256_loadu_ps(&vector[j]);
            // Perform element-wise multiplication and add to sum
            sum = _mm256_add_ps(sum, _mm256_mul_ps(mat, vec));
            //sum = _mm256_fmadd_ps(mat, vec, sum); // Fails in linux build
            jend = j;
        }
        // Horizontal addition of the 8 elements in sum
        float temp[8];
        _mm256_storeu_ps(temp, sum);
        val += temp[0] + temp[1] + temp[2] + temp[3] + temp[4] + temp[5] + temp[6] + temp[7];
        // Handle remaining elements
        for (int j = jend; j < M; j++) {
            val += matrix[i * M + j] * vector[j];
        }
        // ReLU
        if (relu) {
            if (val < 0) { val = 0; }
        }
        result[i] = val;
    }
}
# else
void matrix_vector_multiply(const float* matrix, const float* vector, const float* bias, float* result, const int N, const int M, const bool relu) {
    #pragma omp parallel for
    for (int i = 0; i < N; i++) {
        float val = bias[i];
        for (int j = 0; j < M; j++) {
            val += matrix[i * M + j] * vector[j];
        }
        // ReLU
        if (relu) {
            if (val < 0) { val = 0; }
        }
        result[i] = val;
    }
}
#endif


#pragma warning( disable : 4305 4244 )

// Sets up the config object and puts all the correct lengths into it
struct IsoConfig SetupConfig(const int maxz, const int pres) {
    struct IsoConfig config;
    config.pres = pres;
    config.maxz = maxz;
    config.elen = maxz * pres;
    config.l1 = config.elen * config.elen;
    config.l2 = config.elen;
    config.l3 = config.elen * maxz;
    config.l4 = maxz;
    config.verbose = 0;
    config.dlen = 0;
    return config;
}

// Default settings for the isodec library
struct IsoSettings DefaultSettings() {
    struct IsoSettings settings;
    settings.phaseres = 8;
    settings.verbose = 0;
    settings.peakwindow = 80;
    settings.matchtol = 5;
    settings.minpeaks = 3;
    settings.peakthresh = 0.0001f;
    settings.css_thresh = 0.7f;
    settings.maxshift = 3;
    settings.mzwindow[0] = -1.5f;
    settings.mzwindow[1] = 2.5f;
    settings.plusoneintwindow[0] = 0.1f;
    settings.plusoneintwindow[1] = 0.6f;
    settings.knockdown_rounds = 5;
    settings.min_score_diff = 0.1f;
    settings.minareacovered = 0.25f;
    settings.isolength = 64;
    settings.mass_diff_c = 1.0033;
    settings.adductmass = 1.007276467f;
    settings.minusoneaszero = 1;
    settings.isotopethreshold = 0.01f;
    settings.datathreshold = 0.05f;
    settings.zscore_threshold = 0.95f;
    return settings;
}



// Sets up the weights structure and allocates memory for the weights and biases
struct Weights SetuptWeights(const struct IsoConfig config) {
    struct Weights weights;
    // ReSharper disable once CppDFAMemoryLeak
    weights.w1 = (float *) calloc(config.l1, sizeof(float));
    // ReSharper disable once CppDFAMemoryLeak
    weights.b1 = (float *) calloc(config.elen, sizeof(float));
    // ReSharper disable once CppDFAMemoryLeak
    weights.w2 = (float *) calloc(config.l3, sizeof(float));
    // ReSharper disable once CppDFAMemoryLeak
    weights.b2 = (float *) calloc(config.maxz, sizeof(float));

    // Check that they are not null
    if (weights.w1 == NULL || weights.b1 == NULL || weights.w2 == NULL || weights.b2 == NULL) {
        printf("Error allocating memory for weights\n");
        exit(100);
    }

    return weights;
}

// Frees the memory allocated for the weights and biases
void FreeWeights(const struct Weights weights) {
    free(weights.w1);
    free(weights.b1);
    free(weights.w2);
    free(weights.b2);
}

// Indexing function for 2D arrays
int index2D(const int ncols, const int r, const int c) {
    return r * ncols + c;
}

// Finds the maximum value in an array
float Max(const float *array, const int length) {
    float max = 0;
    for (int i = 0; i < length; i++) {
        if (array[i] > max) {
            max = array[i];
        }
    }
    return max;
}

// Finds the index of the maximum value in an array
int ArgMax(const float *array, const int length) {
    float maxval = 0;
    int maxpos = 0;
    //#pragma omp parallel for reduction(max:maxval)
    for (int i = 0; i < length; i++) {
        if (array[i] > maxval) {
            maxval = array[i];
            maxpos = i;
        }
    }
    return maxpos;
}

float round_to_decimal(const float val, const int decimals) {
    const float factor = powf(10.0f, (float) decimals);
    return roundf(val * factor) / factor;
}


// Sums the values in an array
float Sum(const float *list, const int length) {
    float sum = 0;
    for (int i = 0; i < length; i++) {
        sum += list[i];
    }
    return sum;
}

struct Weights parse_weights(const float *weights, const struct IsoConfig config) {
    const struct Weights parsedweights = SetuptWeights(config);
    // Parse the weights into the structure
    memcpy(parsedweights.w1, weights, config.l1 * sizeof(float));
    memcpy(parsedweights.b1, weights + config.l1, config.l2 * sizeof(float));
    memcpy(parsedweights.w2, weights + config.l1 + config.l2, config.l3 * sizeof(float));
    memcpy(parsedweights.b2, weights + config.l1 + config.l2 + config.l3, config.l4 * sizeof(float));
    return parsedweights;
}

struct Weights load_default_weights(const struct IsoConfig config) {
    const int length = config.l1 + config.l2 + config.l3 + config.l4;
    // printf("Model Number: %d\n", config.pres);
    // printf("Lengths: %d %d %d\n", length, phase_model_8_bin_len, phase_model_4_bin_len);

    float *default_weights = calloc(length, sizeof(float));
    if (default_weights == NULL) {
        printf("Error allocating memory for default weights\n");
        exit(100);
    }
    // Copy from phase_model_8_bin or phase_model_4_bin into default_weights
    if (config.pres == 8) {
        memcpy(default_weights, phase_model_8_bin, length * sizeof(float));
    } else {
        memcpy(default_weights, phase_model_4_bin, length * sizeof(float));
    }

    const struct Weights parsedweights = parse_weights(default_weights, config);
    free(default_weights);
    if (config.verbose == 1) { printf("Loaded Weights from Default Model: %d\n", config.pres); }
    return parsedweights;
}

// Loads the weights from a file into the weights structure
struct Weights load_weights(const char *filename, const struct IsoConfig config) {
    // ReSharper disable once CppDeprecatedEntity
    FILE *file = fopen(filename, "rb");
    if (file == NULL) {
        printf("Error opening weights file: %s\n", filename);
        exit(1);
    }
    const int length = config.l1 + config.l2 + config.l3 + config.l4;
    float *weights = calloc(length, sizeof(float));

    // Read in Weights
    if (weights) {
        if (fread(weights, sizeof(float), length, file) == 1) {
            printf("Error reading weights from file\n");
            exit(1);
        }
        fclose(file);
    } else {
        printf("Error allocating memory in weight loading function\n");
        exit(1);
    }

    const struct Weights parsedweights = parse_weights(weights, config);
    free(weights);

    if (config.verbose == 1) { printf("Loaded Weights from file %s\n", filename); }
    //printf("Loaded Weights from file %s\n", filename);
    return parsedweights;
}

//
// Encodes the data into the encoding matrix
// extern "C" __declspec(dllexport)
int encode(const double *cmz, const float *cint, const int n, float *emat, const struct IsoConfig config,
           const struct IsoSettings settings) {
    //printf("Encoding\n");
    //#pragma omp parallel for
    const float dmax = Max(cint, n);
    const float thresh = dmax * settings.datathreshold;
    int count = 0;
    double modint = 1;
    for (int i = 0; i < n; i++) {
        const float val = cint[i];
        if (val < thresh) { continue; }
        const double rescale = cmz[i] / settings.mass_diff_c;
        count++;
        for (int j = 0; j < config.maxz; j++) {
            const double phase = modf(rescale * (j + 1), &modint);
            const int phaseindex = (int) floor(phase * (double) config.pres);
            emat[index2D(config.pres, j, phaseindex)] += val;
        }
    }

    if (count < 2) { return 0; }
    //Normalize
    const float max = Max(emat, config.elen);
    //#pragma omp parallel for
    for (int i = 0; i < config.elen; i++) {
        emat[i] = emat[i] / max;
    }
    return 1;
}


// Predicts the charge state of the data
int predict_nn(const float *emat, const struct Weights weights, float *h1, float *h2, const struct IsoConfig config,
               const float zscore_threshold, int *z2) {
    matrix_vector_multiply(weights.w1, emat, weights.b1, h1, config.elen, config.elen, 1);
    matrix_vector_multiply(weights.w2, h1, weights.b2, h2, config.maxz, config.elen, 0);
    //  float val = 0;
    // //Calculate h1
    //  #pragma omp parallel for private(val)
    //      for (int i = 0; i < config.elen; i++) {
    //          val = weights.b1[i];
    //          for (int j = 0; j < config.elen; j++) {
    //              val += weights.w1[index2D(config.elen, i, j)] * emat[j];
    //          }
    //
    //          // ReLU
    //          if (val < 0) { val = 0; }
    //
    //          h1[i] = val;
    //      }
    //      //Calculate h2
    //  #pragma omp parallel for private(val)
    //      for (int i = 0; i < config.maxz; i++) {
    //          val = weights.b2[i];
    //          for (int j = 0; j < config.elen; j++) {
    //              val += weights.w2[index2D(config.elen, i, j)] * h1[j];
    //          }
    //          h2[i] = val;
    //      }
    // Get Max Position
    const int maxpos = ArgMax(h2, config.maxz);

    // Get Max Value
    const float maxval = h2[maxpos];

    // Find any peaks above threshold
    const float threshold = maxval * zscore_threshold;
    int maxpos2 = 0;
    float current_best = 0;
    for (int i = 0; i < config.maxz; i++) {
        //printf("z:%i score:%f\n", i, h2[i]);
        if (h2[i] > threshold && i != maxpos && h2[i] > current_best) {
            maxpos2 = i;
            current_best = h2[i];
        }
    }

    *z2 = maxpos2;
    //printf("Predicted Charge: %d\n", maxpos);
    return maxpos;
}

// Simple charge prediction library for a single input. Inputs are mz data, intensities, length, weight file name, and charge output pointer
// extern "C" __declspec(dllexport)
void predict_charge(const double *cmz, const float *cint, const int n, const char *fname, int *charge) {
    // Set up defaults
    const struct IsoSettings settings = DefaultSettings();
    struct IsoConfig config = SetupConfig(50, settings.phaseres);

    // Set up arrays
    config.dlen = n;

    char filename[500];
    // ReSharper disable once CppDeprecatedEntity
    strcpy(filename, fname);

    float *emat = calloc(config.elen, sizeof(float));

    // Load Weights
    // ReSharper disable once CppDFAMemoryLeak
    const struct Weights weights = load_weights(filename, config);

    // Set up intermediate layers for NN
    float *h1 = calloc(config.elen, sizeof(float));
    float *h2 = calloc(config.maxz, sizeof(float));
    int z2 = 0;

    // Check that not null
    if (emat == NULL || h1 == NULL || h2 == NULL) {
        printf("Error allocating memory\n");
        exit(101);
    }

    // Encode the data and predict charge
    const int good = encode(cmz, cint, n, emat, config, settings);
    if (good == 0) {
        *charge = 0;
        printf("Error, not enough peaks above data threshold\n");
    } else {
        *charge = predict_nn(emat, weights, h1, h2, config, settings.zscore_threshold, &z2);
    }

    // Free memory
    free(emat);
    free(h1);
    free(h2);
    FreeWeights(weights);
}

//
//Fast way of finding the nearest data point in an ordered list.
int nearfast(const double *dataMZ, const float point, const int numdat) {
    int start = 0;
    int length = numdat - 1;
    int end = 0;
    int diff = length - start;
    while (diff > 1) {
        if (point < dataMZ[start + (length - start) / 2]) {
            length = start + (length - start) / 2;
        } else if (point == dataMZ[start + (length - start) / 2]) {
            end = start + (length - start) / 2;
            // length = start + (length - start) / 2;
            // start = start + (length - start) / 2;
            return end;
        } else if (point > dataMZ[start + (length - start) / 2]) {
            start = start + (length - start) / 2;
        }
        diff = length - start;
    }
    if (fabsf(point - (float) dataMZ[start]) >= fabsf(point - (float) dataMZ[length])) {
        end = length;
    } else {
        end = start;
    }
    return end;
}

//Fast way of finding the nearest data point in an ordered list.
int nearfast_int(const int *existing, const int point, const int numdat) {
    if (numdat == 0) { return 0; }

    int left = 0;
    int right = numdat - 1;
    while (left < right) {
        const int mid = left + (right - left) / 2;

        if (existing[mid] < point) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    return left;
}

void insert_index(const int index, int *zeroindices, const int numzeros) {
    const int insert = nearfast_int(zeroindices, index, numzeros);
    if (insert == numzeros - 1) {
        zeroindices[insert] = index;
    } else {
        for (int i = numzeros - 1; i >= insert; i--) {
            zeroindices[i + 1] = zeroindices[i];
        }
        zeroindices[insert] = index;
    }
}

int remove_zeros(const int *zeroindices, const int zeros, double *mz, float *inten, const int datalen) {
    for (int i = 0; i <= zeros - 2; i++) {
        const int begin = zeroindices[i] + 1;
        const int end = zeroindices[i + 1];
        const int shift = i + 1;
        for (int j = begin; j < end; j++) {
            mz[j - shift] = mz[j];
            inten[j - shift] = inten[j];
        }
    }

    //The end of the array must be handled separately.
    for (int i = zeroindices[zeros - 1] + 1; i < datalen; i++) {
        mz[i - zeros] = mz[i];
        inten[i - zeros] = inten[i];
    }


    return datalen - zeros;
}

// Function to check if a point is a peak. Used in peak detection
int is_peak_index(const double *dataMZ, const float *dataInt, const int lengthmz, const int window, const float thresh,
                  const int index) {
    const float xval = (float) dataMZ[index]; // Recast to float
    const float yval = dataInt[index];
    //Check if below threshold
    if (yval <= thresh) { return 0; }
    //Check if local max
    int start = index - window;
    int end = index + window + 1;

    if (start < 0) { start = 0; }
    if (end > lengthmz) { end = lengthmz; }

    //float newmax = 0;
    for (int i = start; i < end; i++) {
        if (i != index) {
            if (dataInt[i] > yval) {
                return 0;
            }
            if (dataInt[i] == yval) {
                if (dataMZ[i] < xval) {
                    return 0;
                }
            }
        }
    }

    return 1;
}

// Function to detect peaks in a spectrum. Returns the number of peaks found and fills the peakx and peaky arrays with the peak data
int peak_detect(const double *dataMZ, const float *dataInt, const int lengthmz, const int window, const float thresh,
                float *peakx, float *peaky) {
    int plen = 0;
    const float max = Max(dataInt, lengthmz);
    const float absthresh = thresh * max;
    for (int i = 0; i < lengthmz; i++) {
        if (is_peak_index(dataMZ, dataInt, lengthmz, window, absthresh, i) == 1) {
            //printf("Peak %d: %f %f\n", plen, dataMZ[i], dataInt[i]);
            peakx[plen] = (float) dataMZ[i];
            peaky[plen] = dataInt[i];
            plen++;
        }
    }
    return plen;
}

// Function to calculate the isotope intensity distribution for a given mass
int isotope_dist(const float mass, float *isovals, const float normfactor,
                 const struct IsoSettings settings) {

    int offset = 0;
    if (settings.minusoneaszero == 1) {offset = 1;}

    float max;
    if (mass > 60000) {
        printf("Warning: Mass is very high, may not be accurate: %f\n", mass);
        max = isomike(mass, isovals, settings.isolength, offset);
    }
    else {
        max = isogenmass_fancy(mass, isovals, settings.isolength, offset);
    }

    int realisolength = 0;
    if (max > 0) {
        for (int k = 0; k < settings.isolength; k++) {
            const float yval = isovals[k] / max;
            if (yval > settings.isotopethreshold || (settings.minusoneaszero == 1 && k == 0)) {
                isovals[k] = yval * normfactor;
                realisolength++;
            } else {
                isovals[k] = 0;
            }
        }
    }
    return realisolength;
}

// Sets the mz values for the isotopes
void iso_mz_vals(const float mass, float *mzvals, const int z, const struct IsoSettings settings,
                 const int realisolength) {
    int offset = 0;
    if (settings.minusoneaszero == 1) {
        offset = -1;
    }
    const float mass_diff_c = (float) settings.mass_diff_c; // Recast to float
    for (int i = 0; i < realisolength; i++) {
        const float zfloat = (float) z;
        const float newmass = mass + (float) (i + offset) * mass_diff_c;
        mzvals[i] = (newmass + zfloat * settings.adductmass) / zfloat;
        // Removed check here on if z is 0 because that is checked earlier
    }
}

// Function used in qsorting of arrays from low to high
int compare_float(const void *a, const void *b) {
    if (*(float *) a > *(float *) b) { return 1; }
    if (*(float *) a < *(float *) b) { return -1; }
    return 0;
}

// Function used in qsorting of arrays from high to low
int compare_float_reverse(const void *a, const void *b) {
    if (*(float *) a > *(float *) b) { return -1; }
    if (*(float *) a < *(float *) b) { return 1; }
    return 0;
}

// Function to adjust the isotope distribution to match the experimental peaks
int bestshift_adjust(const double *mz, const float *inten, const int length, float *isomz, float *isodist,
                     const int maxshift,
                     const int charge, struct MatchedPeak *matchpeaks, const int nmatched, const float peakmz,
                     const float mass,
                     const struct IsoSettings settings, int realisolength) {
    float z = (float) charge;
    float mass_diff_c = (float) settings.mass_diff_c;

    // Test a range of shift values and find the one with the most overlap
    int bestshift = 0;
    float bestscore = 0;

    float ppmtol = settings.matchtol * peakmz / 1000000;

    float *scores;
    float *monoisos;
    int lshifts = 2 * maxshift + 1;
    scores = (float *) calloc(lshifts, sizeof(float));
    monoisos = (float *) calloc(lshifts, sizeof(float));

    if (settings.verbose) { printf("Peak: %f %d\n", peakmz, charge); }

    for (int i = -maxshift; i < maxshift + 1; i++) {
        float shift = (float) i * mass_diff_c / z;
        float score = 0;

        float ab = 0;
        float a2 = 0;
        float b2 = 0;

        for (int k = 0; k < realisolength; k++) {
            float isoval = isodist[k];

            if (isoval == 0 && !(settings.minusoneaszero == 1 && k == 0)) { continue; }

            float shiftedmz = isomz[k] + shift;

            float matchedint = 0;

            for (int j = 0; j < length; j++) {
                float diff = (float) mz[j] - shiftedmz;
                if (fabsf(diff) < ppmtol) {
                    float expval = inten[j];
                    if (expval > matchedint) { matchedint = expval; }
                    if (expval > matchedint) { matchedint = expval; }
                }
                if (diff > ppmtol) { break; }
            }

            ab += isoval * matchedint;
            a2 += isoval * isoval;
            b2 += matchedint * matchedint;
        }
        //score += inten[j] * isodist[k] * inten[j] * isodist[k];
        float monoiso = 0;
        if (ab != 0 && a2 != 0 && b2 != 0) {
            score = sqrtf(ab * ab / (a2 * b2));
            if (settings.minusoneaszero == 1) {
                monoiso = (isomz[1] - settings.adductmass) * z + (float) i * mass_diff_c;
            } else { monoiso = (isomz[0] - settings.adductmass) * z + (float) i * mass_diff_c; }
        }

        monoisos[i + maxshift] = monoiso;
        scores[i + maxshift] = score;

        if (score > bestscore) {
            bestscore = score;
            bestshift = i;
        }
        if (settings.verbose) { printf("Score: %f Shift: %d\n", score, i); }
    }
    // If the best score is not good enough, return 0 and move on
    if (bestscore < settings.css_thresh) {
        free(scores);
        free(monoisos);
        return 0;
    }

    // If the minusoneaszero parameter is set, remove the 0 at the start of the sequence
    if (settings.minusoneaszero == 1) {
        for (int j = 1; j < realisolength; j++) {
            isodist[j - 1] = isodist[j];
            isomz[j - 1] = isomz[j];
        }
        realisolength--;
    }

    int nmatches = 0;
    float areamatch = 0;
    float areaisototal = 0;
    float areaexptotal = 0;
    float avgnum = 0;
    float peakint = 0;
    float int1 = -1;
    float int2 = -1;
    int *matchedindsiso;
    matchedindsiso = (int *) calloc(realisolength, sizeof(int));
    int *matchedindsexp;
    matchedindsexp = (int *) calloc(realisolength, sizeof(int));
    float *matchedints;
    matchedints = (float *) calloc(realisolength, sizeof(float));

    for (int j = 0; j < realisolength; j++) {
        float shiftmz = isomz[j] + (float) bestshift * mass_diff_c / z;
        // Adjust the mz values for the best shift
        isomz[j] = shiftmz;
        // Calculate the area of the isotope distribution and weighted average m/z with these values
        areaisototal += isodist[j];
        avgnum += isodist[j] * shiftmz;

        float isoval = isodist[j];
        if (isoval == 0) { continue; }

        float highestmatch = 0;
        int highestindex = -1;
        int mid_index = realisolength / 2;
        // Find matches and zero out the intensities
        for (int i = 0; i < length; i++) {
            float val = inten[i];
            if (val == 0) { continue; } // Skip if intensity is already zero (matched)
            if (j == mid_index) { areaexptotal += val; }

            // Check if it's a match
            float diff = (float) mz[i] - shiftmz;
            if (fabsf(diff) < ppmtol) {
                // If it's the first or second isotope, save the intensity for z=1 checking
                if (j == 0 && val > int1) {
                    int1 = val;
                }
                if (j == 1 && val > int2) {
                    int2 = val;
                }

                // Find max matched peak intensity
                if (val > peakint) { peakint = val; }
                if (val > highestmatch) {
                    highestmatch = val;
                    highestindex = i;
                }
            }
        }
        if (highestmatch > 0 && highestindex >= 0) {
            // Add to area matched and increment nmatches
            areamatch += isoval;
            matchedindsiso[nmatches] = j;
            matchedindsexp[nmatches] = highestindex;
            matchedints[nmatches] = highestmatch;
            nmatches++;
        }
        if (settings.verbose) {
            printf("Matched: %d %f %d %f %f %f\n", j, shiftmz, highestindex, mz[highestindex], ppmtol,
                   inten[highestindex]);
        }
    }


    // Check if the highest 3 peaks match
    bool highest3match = 0;
    if (nmatches >= 3) {
        float *sortedints;
        sortedints = (float *) calloc(length, sizeof(float));
        memcpy(sortedints, inten, length * sizeof(float));
        qsort(sortedints, length, sizeof(float), compare_float_reverse);

        qsort(matchedints, nmatches, sizeof(float), compare_float_reverse);

        bool t1 = sortedints[0] == matchedints[0];
        bool t2 = sortedints[1] == matchedints[1];
        bool t3 = sortedints[2] == matchedints[2];
        highest3match = t1 && t2 && t3;

        free(sortedints);
    }

    free(matchedints);

    //
    //  Check Peaks and Add to Structure if Good
    //

    int isgood = 0;

    // Calculate the area covered by the matched peaks
    float areacovered = 0;
    if (areaexptotal != 0) { areacovered = areamatch / areaexptotal; }

    // Check if the peaks clear the noise threshold
    if (areacovered > settings.minareacovered || highest3match) {
        // Then, check if the number of peaks matches is good enough
        if (nmatches >= settings.minpeaks) {
            isgood = 1;
        } else if (nmatches == 2) {
            // Special cases if there are only two peaks and the charge state is one
            if (z == 1) {
                // If the charge state is one, check if the ratio of the first two matched peaks is within a certain range
                // These are required above to be the monoisotopic peak and next isotope
                if (int1 > 0 && int2 > 0) {
                    float ratio = int2 / int1;
                    if (ratio > settings.plusoneintwindow[0] && ratio < settings.plusoneintwindow[1]) {
                        isgood = 1;
                    }
                }
            }
        }
    }

    if (matchpeaks == NULL) {
        printf("Error: Matched Peaks Structure is NULL\n");
        exit(103);
    }

    float current_peak_score = matchpeaks[nmatched].score;

    if (current_peak_score > 0) {
        if (bestscore > current_peak_score && isgood == 1) {
            isgood = 1;
        } else {
            isgood = 0;
        }
    }


    if (settings.verbose) {
        printf("NMatches: %d Areacovered: %f of %f Top3: %d Is Good: %d\n", nmatches, areacovered,
               settings.minareacovered, highest3match, isgood);
    }
    // If the peak is good, calculate some things and add it to the matched peaks structure
    if (isgood == 1) {
        float avgmz = avgnum / areaisototal;
        float avgmass = (avgmz - settings.adductmass) * z;
        float massshift = (float) bestshift * mass_diff_c;
        float monoiso = (isomz[0] - settings.adductmass) * z; //Don't need to shift because isomz is already shifted


        float scorethresh = bestscore - settings.min_score_diff;
        for (int k = 0; k < 16; k++) {
            float val = -1.0f;
            if (k < lshifts) {
                float score = scores[k];
                if (score > scorethresh) {
                    val = monoisos[k];
                }
            }
            matchpeaks[nmatched].monoisos[k] = val;
        }

        matchpeaks[nmatched].mz = peakmz + massshift / z;
        matchpeaks[nmatched].z = (int) z;
        matchpeaks[nmatched].monoiso = monoiso;
        matchpeaks[nmatched].peakmass = mass + massshift;
        matchpeaks[nmatched].avgmass = avgmass;
        matchpeaks[nmatched].area = areaisototal;
        matchpeaks[nmatched].peakint = peakint;
        matchpeaks[nmatched].score = bestscore;
        matchpeaks[nmatched].realisolength = realisolength;
        // Fill up arrays in struct
        for (int i = 0; i < realisolength; i++) { matchpeaks[nmatched].isomz[i] = isomz[i]; }
        for (int i = 0; i < realisolength; i++) { matchpeaks[nmatched].isodist[i] = isodist[i]; }
        for (int i = 0; i < realisolength; i++) { matchpeaks[nmatched].isomass[i] = monoiso + (float) i * mass_diff_c; }
        for (int i = 0; i < realisolength; i++) { matchpeaks[nmatched].matchedindsiso[i] = matchedindsiso[i]; }
        for (int i = 0; i < realisolength; i++) { matchpeaks[nmatched].matchedindsexp[i] = matchedindsexp[i]; }
    }

    // Free Memory
    free(scores);
    free(monoisos);
    free(matchedindsiso);
    free(matchedindsexp);
    // Return whether peak was good for tally downstream
    return isgood;
}

int optimize_shift(const int z, const double *mz, const float *inten, const int length, struct MatchedPeak *matchpeaks,
                   const int nmatched,
                   const float peakmz, const struct IsoSettings settings) {
    // Set some parameters for max shifts allowed, which depend on charge state
    int maxshift = settings.maxshift;
    if (z < 3) { maxshift = 1; } else if (z < 6) { maxshift = 2; }
    // Find the peak with the highest intensity and calculate mass
    const float maxval = Max(inten, length);
    //int maxpos = ArgMax(inten, length);
    //float peakmz = mz[maxpos];
    const float mass = (peakmz - settings.adductmass) * (float) z;

    // Simulate isotope intensity distribution for this mass
    float *isovals = calloc(settings.isolength, sizeof(float));
    // Check that it is not null
    if (isovals == NULL) { return 0; }

    const int realisolength = isotope_dist(mass, isovals, maxval, settings);

    // Calculate the mz values for the isotopes
    float *isomzs = calloc(settings.isolength, sizeof(float));
    // Check that it is not null
    if (isomzs == NULL) { return 0; }
    iso_mz_vals(mass, isomzs, z, settings, realisolength);

    // Correct shift of mass values to match up peakmz
    const int isomaxpos = ArgMax(isovals, realisolength);
    const float peakiso = isomzs[isomaxpos];

    // Calculate shift to align the most abundant peaks of the centroids and the simulated isotope distribution
    const float shift = peakmz - peakiso;
    for (int i = 0; i < realisolength; i++) {
        isomzs[i] = isomzs[i] + shift;
    }

    int isgood = 0;
    if (matchpeaks != NULL) {
        isgood = bestshift_adjust(mz, inten, length, isomzs, isovals, maxshift, z, matchpeaks,
                                  nmatched, peakmz, mass, settings, realisolength);
    }

    free(isovals);
    free(isomzs);
    return isgood;
}

void print_matches(const struct MatchedPeak *matchedpeaks, const int nmatched) {
    for (int i = 0; i < nmatched; i++) {
        printf("Peak %d: %f Charge: %d ", i, matchedpeaks[i].mz, matchedpeaks[i].z);
        printf("Monoisotopic Mass: %f ", matchedpeaks[i].monoiso);
        //printf("Peak Mass: %f\n", matchedpeaks[i].peakmass);
        //printf("Average Mass: %f\n", matchedpeaks[i].avgmass);
        printf("Area: %f\n", matchedpeaks[i].area);
        //printf("Peak Intensity: %f\n", matchedpeaks[i].peakint);
    }
}


// extern "C" __declspec(dllexport)
int process_spectrum(const double *cmz, const float *cint, int n, const char *fname, struct MatchedPeak *matchedpeaks,
                     struct IsoSettings settings) {
    // Set some parameters
    struct IsoConfig config = SetupConfig(50, settings.phaseres);
    config.dlen = n;
    config.verbose = settings.verbose;
    if (config.verbose == 1) { printf("Processing Spectrum. Phaseres: %d Length: %d\n", settings.phaseres, n); }

    //IsoSettings settings = DefaultSettings();

    // Load Weights
    // ReSharper disable once CppDFAMemoryLeak
    struct Weights weights;
    if (fname == NULL) {
        weights = load_default_weights(config);
    } else {
        // Set weights file name
        char filename[500];
        // ReSharper disable once CppDeprecatedEntity
        strcpy(filename, fname);
        weights = load_weights(filename, config);
    }


    // Set up memory for peaks
    float *peakx;
    float *peaky;
    peakx = (float *) calloc(n, sizeof(float));
    peaky = (float *) calloc(n, sizeof(float));

    // Create array of MatchPeaks for filling downstream
    //MatchedPeak* matchedpeaks;
    //matchedpeaks = (MatchedPeak*)calloc(n, sizeof(MatchedPeak));
    int nmatched = 0;

    // create copy of cint
    float *cintcopy;
    cintcopy = (float *) calloc(n, sizeof(float));
    memcpy(cintcopy, cint, n * sizeof(float));

    // create a copy of cmz
    double *cmzcopy;
    cmzcopy = (double *) calloc(n, sizeof(double));
    memcpy(cmzcopy, cmz, n * sizeof(double));

    //Create an empty list of zeroindices and track how many zeros have been added.
    int *zeroindices;
    zeroindices = (int *) calloc(n, sizeof(int));
    int zeros = 0;

    // Check if cintcopy or cmzcopy are null
    if (cintcopy == NULL || cmzcopy == NULL || matchedpeaks == NULL
        || peakx == NULL || peaky == NULL || zeroindices == NULL) {
        printf("Error: Could not allocate memory for cintcopy or cmzcopy\n");
        free(zeroindices);
        free(cmzcopy);
        free(cintcopy);
        free(peakx);
        free(peaky);
        FreeWeights(weights);
        exit(102);
    }

    // Loop through knockdown rounds
    for (int k = 0; k < settings.knockdown_rounds; k++) {
        // Adjust settings based on round
        if (k > 0) {
            settings.peakwindow = (int) roundf((float) settings.peakwindow * 0.5f);
            if (settings.peakwindow < 1) { settings.peakwindow = 1; }

            settings.peakthresh = settings.peakthresh * 0.5f;
            if (settings.peakthresh < 0.000001) { settings.peakthresh = 0.000001f; }
        }
        if (k > 5) {
            settings.css_thresh = settings.css_thresh * 0.9f;
            if (settings.css_thresh < 0.6) { settings.css_thresh = 0.6f; }
        }

        // Detect Peaks
        int plen = peak_detect(cmzcopy, cintcopy, n, settings.peakwindow, settings.peakthresh, peakx, peaky);
        //Should rewrite this to include indexes
        if (config.verbose == 1) { printf("Round %d: Peaks: %d Window: %d\n", k, plen, settings.peakwindow); }
        if (plen == 0) { break; }

        // Loop through all peaks. Find window. If there is enough peaks in the window, predict the charge.
        int ngood = 0;
        //# pragma omp parallel for reduction(+:zeros, ngood, nmatched) schedule(dynamic) shared(zeros, nmatched, matchedpeaks, nmatched, zeroindices) private(plen, peakx, peaky, cmzcopy, cintcopy, n, weights, config, settings)
        for (int i = 0; i < plen; i++) {
            //printf("Peak %d: %f %f \n", i, peakx[i], peaky[i]);
            const float peakmz = peakx[i];
            const int peakindex = nearfast(cmzcopy, peakmz, n);
            // Check if the peak is already zeroed out
            if (cintcopy[peakindex] == 0) continue;

            // Find the window around the peak
            const float lowmz = peakmz + settings.mzwindow[0];
            const float highmz = peakmz + settings.mzwindow[1];
            int start = nearfast(cmzcopy, lowmz, n);
            int end = nearfast(cmzcopy, highmz, n) + 1;
            // Correct window if there are two peaks selected in the same window
            if (i < plen - 1) {
                const float nextpeakmz = peakx[i + 1];
                const int nextpeakindex = nearfast(cmzcopy, nextpeakmz, n);
                if (nextpeakindex <= end) {
                    end = nextpeakindex - 1;
                }
            }
            //printf("Lowmz: %f Peak: %f Highmz: %f\n", lowmz, peakmz, highmz);

            const int l = end - start;
            // This should be settings.minpeaks, but it is possible for 2 peaks to be assigned to a +1 charge, so we need to check for at least 2 peaks
            if (l < 2) {
                if (cintcopy[peakindex] != 0) {
                    // If the peak is not good, zero out the peak intensity in the copy array and store the index
                    zeros++;
                    insert_index(peakindex, zeroindices, zeros);
                }
                continue;
            }

            // Declare memory for weights and biases and emat and mz and inten arrays
            float *emat = calloc(config.elen, sizeof(float));
            float *h1 = calloc(config.elen, sizeof(float));
            float *h2 = calloc(config.maxz, sizeof(float));
            double *mz = calloc(l, sizeof(double));
            float *inten = calloc(l, sizeof(float));
            // Check if they are null
            if (mz == NULL || inten == NULL || emat == NULL || h1 == NULL || h2 == NULL) {
                printf("Error: Could not allocate memory for mz or inten\n");
                exit(105);
            }

            // Note, need to think about the zeros in here
            for (int j = 0; j < l; j++) {
                mz[j] = cmzcopy[start + j];
                inten[j] = cintcopy[start + j];
            }

            // Encode the data, this will check if enough peaks have met the data threshold
            int good = encode(mz, inten, l, emat, config, settings);
            int z2 = 0;
            // If good, predict the charge
            int z = 0;
            if (good == 1) { z = predict_nn(emat, weights, h1, h2, config, settings.zscore_threshold, &z2); }

            if (z > 0) {
                //if (z2 > 0) { printf("Z1:%i Z2: %i\n", z, z2); }
                // Optimize the shift and add to matched peaks if good
                int isgood = optimize_shift(z, mz, inten, l, matchedpeaks, nmatched, peakmz, settings);

                if (z2 > 0) {
                    int isgood2 = optimize_shift(z2, mz, inten, l, matchedpeaks, nmatched, peakmz, settings);
                    if (isgood2 == 1) {
                        isgood = 1;
                    }
                }

                if (isgood == 1) {
                    // Finally, copy the new intensities to the original intensity array, which will zero out the matched peaks
                    for (int j = 0; j < matchedpeaks[nmatched].realisolength; j++) {
                        const int index = matchedpeaks[nmatched].matchedindsexp[j];
                        if (index == 0) {
                            if (j == 0) {
                                zeros++;
                                insert_index(index + start, zeroindices, zeros);
                            } else { break; }
                        } else if (inten[index] != 0) {
                            zeros++;
                            insert_index(index + start, zeroindices, zeros);
                        }
                    }

                    // Save the start and end indexes for the peak
                    start = nearfast(cmz, lowmz, config.dlen);
                    end = nearfast(cmz, highmz, config.dlen) + 1;
                    matchedpeaks[nmatched].startindex = start;
                    matchedpeaks[nmatched].endindex = end;

                    //printf("Matched %d: %f: Charge: %d\n", i, matchedpeaks[nmatched].mz, z);
                    nmatched++;
                    ngood++;
                } else {
                    if (cintcopy[peakindex] != 0) {
                        // If the peak is not good, zero out the peak intensity in the copy array and store the index
                        zeros++;
                        insert_index(peakindex, zeroindices, zeros);
                    }
                }
            }
            free(mz);
            free(inten);
            free(emat);
            free(h1);
            free(h2);
            // End of Peak Loop
        }
        if (config.verbose == 1) { printf("\nGood: %d\n", ngood); }
        if (ngood == 0) { break; }
        //printf("n=%d", n);
        n = remove_zeros(zeroindices, zeros, cmzcopy, cintcopy, n);
        zeros = 0;
        // End of Knockdown Round Loop
    }

    if (config.verbose == 1) {
        print_matches(matchedpeaks, nmatched);
        printf("%d Peaks Matched!\n", nmatched);
    }

    free(zeroindices);
    free(cintcopy);
    free(cmzcopy);
    free(peakx);
    free(peaky);
    FreeWeights(weights);
    if (config.verbose == 1) { printf("Done\n"); }
    return nmatched;
}


// count the number of lines we have in datafile
int getfilelength(const char *infile) {
    int l = 0;
    // ReSharper disable once CppDeprecatedEntity
    FILE *file_ptr = fopen(infile, "r");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    int read = -2;
    while (EOF != read) {
        char input[501];
        // ReSharper disable once CppDeprecatedEntity
        read = fscanf(file_ptr, "%500[^\n]%*c", input);
        if (read != 0 && read != EOF) {
            l += 1;
        } else {
            break;
        }
    }

    fclose(file_ptr);
    return l;
}


//Reads in x y file.
void readfile(char *infile, const int lengthmz, double *dataMZ, float *dataInt) {
    char x[500] = {0};
    char y[500] = {0};

    // ReSharper disable once CppDeprecatedEntity
    FILE *file_ptr = fopen(infile, "r");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    for (int i = 0; i < lengthmz; i++) {
        // ReSharper disable once CppDeprecatedEntity
        const int match = fscanf(file_ptr, "%s %s", x, y);
        if (match != 0) {
            dataMZ[i] = strtod(x, NULL);
            dataInt[i] = strtof(y, NULL);
        }
    }

    fclose(file_ptr);
}

// Write the matched peaks to a file
void write_matches(const struct MatchedPeak *matchedpeaks, const int nmatched, const char *outfile) {
    // ReSharper disable once CppDeprecatedEntity
    FILE *file_ptr = fopen(outfile, "w");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", outfile);
        exit(1);
    }
    if (matchedpeaks != NULL) {
        for (int i = 0; i < nmatched; i++) {
            fprintf(file_ptr, "%f %d %f \n", matchedpeaks[i].mz,
                    matchedpeaks[i].z, matchedpeaks[i].monoiso);
        }
    }

    fclose(file_ptr);
}

// extern "C" __declspec(dllexport)
int process_spectrum_default(const double *cmz, const float *cint, const int n, const char *fname,
                             struct MatchedPeak *matchedpeaks) {
    const struct IsoSettings settings = DefaultSettings();
    return process_spectrum(cmz, cint, n, fname, matchedpeaks, settings);
}

void run(char *filename, char *outfile, const char *weightfile) {
    // printf("Starting Dll. Threads: %d\n", omp_get_max_threads());

    // char filename[500];
    // // ReSharper disable once CppDeprecatedEntity
    // strcpy(filename, "C:\\Data\\IsoNN\\test2.txt");

    const int length = getfilelength(filename);
    //int length = 23;
    if (length == 0) {
        printf("Error: No Data in File\n");
        exit(1);
    }

    double *cmz = calloc(length, sizeof(double));

    float *cint = calloc(length, sizeof(float));

    // Check cmz and cint are not null
    if (cmz == NULL || cint == NULL) {
        printf("Error: Memory not allocated.\n");
        exit(1);
    }

    readfile(filename, length, cmz, cint);
    printf("Loaded File: %s\n", filename);
    printf("Length: %d\n", length);

    struct MatchedPeak *matchedpeaks = calloc(length, sizeof(struct MatchedPeak));
    printf("Processing Spectrum\n");
    const int nmatched = process_spectrum_default(cmz, cint, length, weightfile, matchedpeaks);

    printf("Writing Matches %s\n", outfile);
    write_matches(matchedpeaks, nmatched, outfile);

    printf("Done: %d\n", nmatched);

    free(cmz);
    free(cint);
    free(matchedpeaks);
}
