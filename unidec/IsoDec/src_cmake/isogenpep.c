#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "fftw3.h"
#include "isogendep.h"
#include "isogenpep.h"
#include "isogenpep_model.h"
#include "isogenprot_model.h"
#include "isogenmass_model_8.h"
#include "isogenmass_model_32.h"
#include "isogenmass_model_64.h"
#include "isogenmass_model_128.h"


double massavgine = 111.1254;
double avgine[11] = {7.7583, 4.9384, 1.3577, 1.4773, 0.0417, 0, 0, 0, 0, 0, 0};
int numaminoacids = 23;

int num_simp_elements = 5;

//{H, C, N, O, S, Fe, K, Ca, Ni, Zn, Mg}
const int aa_vectors[][11] = {
    {5,3,1,1,0,0,0,0,0,0,0},
    {5,3,1,1,1,0,0,0,0,0,0},
    {5,4,1,3,0,0,0,0,0,0,0},
    {7,5,2,3,0,0,0,0,0,0,0},
    {9,9,1,1,0,0,0,0,0,0,0},
    {3,2,1,1,0,0,0,0,0,0,0},
    {7,6,3,1,0,0,0,0,0,0,0},
    {11,6,1,1,0,0,0,0,0,0,0},
    {12,6,2,1,0,0,0,0,0,0,0},
    {11,6,1,1,0,0,0,0,0,0,0},
    {9,5,1,1,1,0,0,0,0,0,0},
    {6,4,2,2,0,0,0,0,0,0,0},
    {7,5,1,1,0,0,0,0,0,0,0},
    {8,5,2,2,0,0,0,0,0,0,0},
    {12,6,4,1,0,0,0,0,0,0,0},
    {7,5,1,2,0,0,0,0,0,0,0},
    {7,4,1,2,0,0,0,0,0,0,0},
    {9,5,1,1,0,0,0,0,0,0,0},
    {10,11,2,1,0,0,0,0,0,0,0},
    {9,9,1,2,0,0,0,0,0,0,0}
};


const char aaorder[] = "ACDEFGHIKLMNPQRSTVWY";

const char *pep_encoding_elements[] = {"H", "C", "N", "O", "S", "Fe", "K", "Ca", "Ni", "Zn", "Mg"};

#define ISO_LEN 32

// Isotope Parameters
float isoparams[10] = {
    1.00840852e+00f, 1.25318718e-03f, 2.37226341e+00f, 8.19178000e-04f, -4.37741951e-01f,
    6.64992972e-04f, 9.94230511e-01f, 4.64975237e-01f, 1.00529041e-02f, 5.81240792e-01f
};

// Function used in isotope distribution calculation
float isotopemid(const float mass, const float *isoparams) {
    const float a = isoparams[4];
    const float b = isoparams[5];
    const float c = isoparams[6];
    return a + b * powf(mass, c);
}

// Function used in isotope distribution calculation
float isotopesig(const float mass, const float *isoparams) {
    const float a = isoparams[7];
    const float b = isoparams[8];
    const float c = isoparams[9];
    return a + b * powf(mass, c);
}

// Function used in isotope distribution calculation
float isotopealpha(const float mass, const float *isoparams) {
    const float a = isoparams[0];
    const float b = isoparams[1];
    return a * expf(-mass * b);
}

// Function used in isotope distribution calculation
float isotopebeta(const float mass, const float *isoparams) {
    const float a = isoparams[2];
    const float b = isoparams[3];
    return a * expf(-mass * b);
}

// Averagine distribution generation by curve fitting.
float pep_mass_to_dist_fitting(const float mass, float * isodist, const int isolen, const int offset) {
    const float mid = isotopemid(mass, isoparams);
    const float sig = isotopesig(mass, isoparams);
    if (sig == 0) {
        printf("Error: Sigma Isotope Parameter is 0");
        exit(102);
    }
    const float alpha = isotopealpha(mass, isoparams);
    const float amp = 1.0f - alpha;
    const float beta = isotopebeta(mass, isoparams);

    float max = 0.0f;
    for (int k = 0; k < isolen-offset; k++) {
        const float kfloat = (float) k;
        const float e = alpha * expf(-kfloat * beta);
        const float g = amp / (sig * 2.50662827f) * expf(-powf(kfloat - mid, 2) / (2 * powf(sig, 2)));
        const float val = e + g;
        if(val> max) {max = val;}
        isodist[k+offset] = val;
    }
    return max;
}


int aa_to_index(char aa) {
    for (int i = 0;i<20;i++) {
        if (aa == aaorder[i]){return i;}
    }
    return -1;
}

//NN
void add_mod_to_vec(char* mod, float* vec) {

    int i = 0;
    printf("Mod:%s\n", mod);
    while (mod[i] != '\0') {
        if (!isupper(mod[i])) {
            printf("Error while parsing modification formula.");
            return;
        }

        char symbol[3] = {0};
        symbol[0] = mod[i++];
        if (islower(mod[i])) {
            symbol[1] = mod[i++];
        }

        float count = 0;
        while (isdigit(mod[i])) {
            count = count * 10 + (mod[i++] - '0');
        }

        if (count == 0) count = 1;

        int symbol_matched = 0;
        for (int j = 0; j < 11; j++) {
            if (strcmp(symbol, pep_encoding_elements[j]) == 0) {
                vec[j+20] += count;
                printf("Added Element: %s, Count %f\n", symbol, count);
                symbol_matched = 1;
                break;
            }
        }
        if (symbol_matched == 0){ printf("Unable to add Element: %s, Count %f\n", symbol, count); }

    }
}

//NN
int pep_seq_to_nnvector(const char* seq, float* vector) {
    int length = strlen(seq);

    int aas = 0;

    int in_mod = 0;
    char curr_mod[100] = {0};
    int mod_index = 0;

    for (int i = 0; i < length; i++) {
        if (in_mod == 0) {
            if (seq[i] == '[') {
                in_mod = 1;
                continue;
            }

            aas += 1;
            int aaindex = aa_to_index(seq[i]);
            if (aaindex != -1){ vector[aaindex] += 1.0f; }
        }
        else {
            if (seq[i] == ']') {
                in_mod = 0;
                curr_mod[mod_index] = '\0';

                add_mod_to_vec(curr_mod, vector);

                mod_index = 0;
                memset(curr_mod, 0, sizeof(curr_mod));
            }
            else {
                curr_mod[mod_index] = seq[i];
                mod_index++;
            }
        }
    }
    return aas;
}


float nn_pep_seq_to_dist(const char* seq, float* isodist, int offset){
    float* vector = (float *) calloc(31, sizeof(float));
    int aas = pep_seq_to_nnvector(seq, vector);

    if (aas > 200) {
        printf("Sequence contains too many amino acids (>200).");
        return -1.0f;
    }

    struct IsoGenWeights weights = SetupWeights(31, 64);
    if (aas >= 1 && aas <= 50) { weights = LoadWeights(weights, isogenpep_model_bin); }
    else{ weights = LoadWeights(weights, isogenprot_model_bin); }
    neural_net(vector, isodist, weights);
    free(vector);

    for (int i = 64 - offset - 1; i >= 0; i--)
    {
        isodist[i + offset] = isodist[i];
        if (i < offset) { isodist[i] = 0.0f; }
    }

    float maxval = 0.0f;
    for (int i = 0; i < 64; i++) {
        if (isodist[i] > maxval) {maxval = isodist[i];}
    }
    return maxval;
}

//fft
void add_mod_to_fftlist(char* mod, int* fftlist) {
    int i = 0;
    while (mod[i] != '\0') {
        if (!isupper(mod[i])) {
            printf("Error while parsing modification formula.");
            return;
        }

        char symbol[3] = {0};
        symbol[0] = mod[i++];
        if (islower(mod[i])) {
            symbol[1] = mod[i++];
        }

        int count = 0;
        while (isdigit(mod[i])) {
            count = count * 10 + (mod[i++] - '0');
        }

        if (count == 0) count = 1;

        int symbol_matched = 0;
        for (int j = 0; j < 11; j++) {
            if (strcmp(symbol, pep_encoding_elements[j]) == 0) {
                fftlist[j] += count;
                printf("Added Element: %s, Count %i\n", symbol, count);
                symbol_matched = 1;
                break;
            }
        }
        if (symbol_matched == 0){ printf("Unable to add Element: %s, Count %i\n", symbol, count); }

    }
}

//fft
void pep_seq_to_fftlist(const char* sequence, int* fftlist)
{
    // Initialize the formulalist to zero but add the elements of water for the terminii
    fftlist[0] = 2; // Hydrogen
    fftlist[1] = 0; // Carbon
    fftlist[2] = 0; // Nitrogen
    fftlist[3] = 1; // Oxygen
    fftlist[4] = 0; // Sulfur
    fftlist[5] = 0; // Iron
    fftlist[6] = 0; // Potassium
    fftlist[7] = 0; // Calcium
    fftlist[8] = 0; // Nickel
    fftlist[9] = 0; // Zinc
    fftlist[10] = 0; // Magnesium


    int length = strlen(sequence);

    int in_mod = 0;
    char curr_mod[100] = {0};
    int mod_index = 0;

    for (int i = 0; i < length; i++) {
        if (in_mod == 0) {
            if (sequence[i] == '[') {
                in_mod = 1;
                continue;
            }

            //Handle the case where the character corresponds to an amino acid.
            int aaindex = aa_to_index(sequence[i]);
            if (aaindex != -1) {
                for (int j = 0;j<num_simp_elements;j++) {
                    fftlist[j] += aa_vectors[aaindex][j];
                }
            }
        }
        else {
            if (sequence[i] == ']') {
                in_mod = 0;
                curr_mod[mod_index] = '\0';

                add_mod_to_fftlist(curr_mod, fftlist);

                mod_index = 0;
                memset(curr_mod, 0, sizeof(curr_mod));
            }
            else {
                curr_mod[mod_index] = sequence[i];
                mod_index++;
            }
        }
    }
}


//fft
void pep_mass_to_fftlist(const float mass, int* fftlist)
{
    int num = (int)round(mass/massavgine);
    for (int i = 0; i < 5; i++)
    {
        fftlist[i] = (int)round(num * avgine[i]);
    }
    for (int i =5;i<11;i++) {
        fftlist[i] = 0;
    }
}


//fft
float fft_pep_seq_to_dist(const char* sequence, float* isodist, const int isolen, const int offset)
{
    int* formulalist = (int*)calloc(11, sizeof(int));
    // Check for null
    if (formulalist == NULL)
    {
        printf("Error: Could not allocate memory for formulalist\n");
        return 1;
    }
    pep_seq_to_fftlist(sequence, formulalist);
    float maxval = fft_list_to_dist(formulalist, isolen, isodist);

    for (int i = isolen - offset - 1; i >= 0; i--)
    {
        isodist[i + offset] = isodist[i];
        if (i < offset) { isodist[i] = 0.0f; }
    }

    if (maxval > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= maxval;
        }
    }

    free(formulalist);
    return maxval;
}

//fft
float fft_pep_mass_to_dist(const float mass, float *isodist, const int isolen, const int offset)
{
    int* fftlist = (int*)calloc(11, sizeof(int));
    pep_mass_to_fftlist(mass, fftlist);
    float max_val = fft_list_to_dist(fftlist, isolen, isodist);

    for (int i = isolen - offset - 1; i >= 0; i--)
    {
        isodist[i + offset] = isodist[i];
        if (i < offset) { isodist[i] = 0.0f; }
    }
    if (max_val > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= max_val;
        }
    }
    return max_val;
}


//nn
float nn_pep_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset) {
    float* vector = (float*)calloc(5, sizeof(float));
    if (vector == NULL) {
        printf("Error: Could not allocate memory for vector\n");
    }

    mass_to_vector(mass, vector);

    struct IsoGenWeights weights = SetupWeights(5, isolen);
    if (isolen == 8){ weights = LoadWeights(weights, isogenmass_model_8_bin); }
    else if ( isolen == 32 ){ weights = LoadWeights(weights, isogenmass_model_32_bin); }
    else if ( isolen == 64 ){ weights = LoadWeights(weights, isogenmass_model_64_bin); }
    else if ( isolen == 128 ){ weights = LoadWeights(weights, isogenmass_model_128_bin); }
    else {
        printf("Unsupported distribution length.");
        return -1.0f;
    }

    neural_net(vector, isodist, weights);
    free(vector);

    for (int i = isolen - offset - 1; i >= 0; i--)
    {
        isodist[i + offset] = isodist[i];
        if (i < offset) { isodist[i] = 0.0f; }
    }

    float maxval = 0.0f;
    for (int i = 0; i < isolen; i++) {
        if (isodist[i] > maxval) {maxval = isodist[i];}
    }

    for (int i = 0; i < isolen; i++) {
        isodist[i] /= maxval;
    }

    return maxval;
}

