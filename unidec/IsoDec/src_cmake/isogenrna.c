#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fftw3.h"
#include "isogendep.h"
#include "isogen_rnaveragine_model32.h"
#include "isogen_rnaveragine_model64.h"
#include "isogen_rnaveragine_model128.h"
#include "isogenrna_model_64.h"
#include "isogenrna_model_128.h"




const char rnaOrder[] = "ACGU";

//All rna averagine values assume 1 phospho per nucleotide
const float rnaveragine_mass = 320.283814;
const double rnaveragine_comp_numerical[] = {9.50, 10.75, 3.75, 7.0, 0};

static const unsigned char nt_lookup[256] = {
    ['A']=1, ['C']=2, ['G']=3, ['U']=4,
    ['a']=1, ['c']=2, ['g']=3, ['u']=4
};

const int rna_vectors[][5] = {
    {10,11,5,6,0},
    {9,11,3,7,0},
    {10,11,5,7,0},
    {9,10,2,8,0}
};

const double rnaveragine_coeff[] = {0.029661, 0.033564, 0.01171, 0.02186, 0};


void rna_mass_to_list(float initialMass, int* fftlist)
//Mass -> List 5 length list of number of {C, H, N, O, S}
{
    //Addition of an extra phospho here improves the calculation because
    //the rna averagine composition includes 1 phospho per nucleotide.
    float m = initialMass + 95.9534;

    for (int i = 0; i < 5; i++)
    {
        fftlist[i] = (int)round(m * rnaveragine_coeff[i]);
    }

    //Correct for one less phospho than monomers
    fftlist[3] -= 4;
}


int fft_rna_len_to_isolen(const int rna_len) {
    if (rna_len < 200) {
        return 64;
    }
    if (rna_len < 501) {
        return 128;
    }
    return 1024;
}


int nn_rna_mass_to_isolen(const float mass) {
    if (mass < 24000) {
        return 32;
    }
    if (mass < 65000) {
        return 64;
    }
    if (mass < 165000) {
        return 128;
    }
    return -1;
}


int fft_rna_mass_to_isolen(float mass) {
    if (mass < 65000) {
        return 64;
    }
    if (mass < 165000) {
        return 128;
    }
    return 1024;
}


int rna_seq_to_vector(const char* seq, float* vector)
//RNA -> Vector 4 length list of number of A C G U
{
    int len = 0;
    for (const char *p = seq; *p != '\0'; p++) {
        char c = *p;

        unsigned char idx = nt_lookup[c];

        if (idx)
            vector[idx-1]++;
        else {
            printf("Unexpected nucleotide: %c\n", *seq);
            return -1;
        }

        len++;
    }
    return len;
}


int rna_seq_to_fftlist(const char* seq, int* fftlist)
{
    // Initialize the formulalist to zero but add the elements of water for the terminii
    fftlist[0] = 0; // Carbon
    fftlist[1] = 0; // Hydrogen
    fftlist[2] = 0; // Nitrogen
    fftlist[3] = 0; // Oxygen
    fftlist[4] = 0; // Sulfur

    int len = 0;

    for (const char *p = seq; *p != '\0'; p++) {
        char c = *p;

        unsigned char idx = nt_lookup[c];

        if (idx) {
            int nt_index = idx - 1;
            fftlist[0] += rna_vectors[nt_index][0];
            fftlist[1] += rna_vectors[nt_index][1];
            fftlist[2] += rna_vectors[nt_index][2];
            fftlist[3] += rna_vectors[nt_index][3];
        }
        else {
            printf("Unexpected nucleotide: %c\n", c);
            return -1;
        }


        len++;
    }

    //The nucleotides in rna_vectors contain a 1 phospho
    //One phospho must be removed to correct for this.
    fftlist[3] -= 4;

    return len;
}

//fft
float fft_rna_mass_to_dist(const float mass, float *isodist, const int isolen, const int offset)
{
    int* fftlist = (int*)calloc(5, sizeof(int));
    if (fftlist == NULL) {
        printf("Could not allocate memory for fftlist.");
        return -1;
    }

    rna_mass_to_list(mass, fftlist);

    int fft_isolen = fft_rna_mass_to_isolen(mass);

    float* fft_isodist = (float*)calloc(fft_isolen, sizeof(float));

    float max_val = fft_list_to_dist(fftlist, fft_isolen, fft_isodist);

    if (fft_isolen < isolen) {
        for (int i = fft_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = fft_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = fft_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }

    free(fft_isodist);
    free(fftlist);

    if (max_val > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= max_val;
        }
    }

    return max_val;
}

//fft
float fft_rna_seq_to_dist(const char* sequence, float* isodist, const int isolen, const int offset)
{
    int* formulalist = (int*)calloc(5, sizeof(int));
    // Check for null
    if (formulalist == NULL)
    {
        printf("Error: Could not allocate memory for formulalist\n");
        return -1;
    }

    int len = rna_seq_to_fftlist(sequence, formulalist);
    if (len == -1) {
        return -1;
    }

    int fft_isolen = fft_rna_len_to_isolen(len);

    float* fft_isodist = (float*)calloc(fft_isolen, sizeof(float));

    float maxval = fft_list_to_dist(formulalist, fft_isolen, fft_isodist);

    if (fft_isolen < isolen) {
        for (int i = fft_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = fft_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = fft_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }

    free(fft_isodist);
    free(formulalist);

    if (maxval > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= maxval;
        }
    }

    return maxval;
}

//nn
float nn_rna_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset) {
    float* vector = (float*)calloc(5, sizeof(float));
    if (vector == NULL) {
        printf("Error: Could not allocate memory for vector\n");
    }

    mass_to_vector(mass, vector);

    int nn_isolen = nn_rna_mass_to_isolen(mass);

    if (nn_isolen == -1) {
        printf("Error: Mass outside of allowed NN mass range: %f\n", mass);
        return -1;
    }

    float* nn_isodist = (float*)calloc(nn_isolen, sizeof(float));

    struct IsoGenWeights weights = SetupWeights(5, nn_isolen);
    if (nn_isolen == 32){ weights = LoadWeights(weights, isogen_rnaveragine_model32_bin); }
    else if ( nn_isolen == 64 ){ weights = LoadWeights(weights, isogen_rnaveragine_model64_bin); }
    else { weights = LoadWeights(weights, isogen_rnaveragine_model128_bin); }

    neural_net(vector, nn_isodist, weights);
    free(vector);
    FreeIsogenWeights(weights);

    if (nn_isolen < isolen) {
        for (int i = nn_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = nn_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = nn_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    free(nn_isodist);

    float maxval = 0.0f;
    for (int i = 0; i < isolen; i++) {
        if (isodist[i] > maxval) {maxval = isodist[i];}
    }

    for (int i = 0; i < isolen; i++) {
        isodist[i] /= maxval;
    }
    return maxval;
}

//nn
float nn_rna_seq_to_dist(const char* seq, float* isodist, int isolen, int offset) {
    float* vector = calloc(4, sizeof(float));
    if (vector == NULL) {
        printf("Error: Could not allocate memory for vector\n");
    }

    int len = rna_seq_to_vector(seq, vector);

    if (len == -1) {
        return -1;
    }

    struct IsoGenWeights weights;
    float* nn_isodist;
    int nn_isolen;


    if (len >= 1 && len <= 200) {
        weights = SetupWeights(4, 64);
        weights = LoadWeights(weights, isogenrna_model_64_bin);
        nn_isolen = 64;
        nn_isodist = (float*)calloc(nn_isolen, sizeof(float));
    }
    if (len >= 201 && len <= 500) {
        weights = SetupWeights(4, 128);
        weights = LoadWeights(weights, isogenrna_model_128_bin);
        nn_isolen = 128;
        nn_isodist = (float*)calloc(nn_isolen, sizeof(float));
    }

    neural_net(vector, nn_isodist, weights);
    free(vector);
    FreeIsogenWeights(weights);

    if (nn_isolen < isolen) {
        for (int i = nn_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = nn_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = nn_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    free(nn_isodist);


    float maxval = 0.0f;
    for (int i = 0; i < isolen; i++) {
        if (isodist[i] > maxval) {maxval = isodist[i];}
    }

    for (int i = 0;i< isolen; i++) {
        isodist[i] /= maxval;
    }
    return maxval;
}

//BRAIN
float brain_rna_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset) {
    int* brain_list = (int*)calloc(5, sizeof(int));
    rna_mass_to_list(mass, brain_list);

    int brain_isolen = fft_rna_mass_to_isolen(mass);

    float* brain_isodist = (float*)calloc(brain_isolen, sizeof(float));

    float max_val = brain_list_to_dist(brain_list, brain_isolen, brain_isodist);

    if (brain_isolen < isolen) {
        for (int i = brain_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = brain_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = brain_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }

    free(brain_isodist);
    free(brain_list);

    if (max_val > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= max_val;
        }
    }

    return max_val;
}

//BRAIN
float brain_rna_seq_to_dist(const char* sequence, float* isodist, const int isolen, const int offset) {
    int* formulalist = (int*)calloc(5, sizeof(int));
    // Check for null
    if (formulalist == NULL)
    {
        printf("Error: Could not allocate memory for formulalist\n");
        return 1;
    }
    int len = rna_seq_to_fftlist(sequence, formulalist);

    int brain_isolen = fft_rna_len_to_isolen(len);

    float* brain_isodist = (float*)calloc(brain_isolen, sizeof(float));

    float maxval = brain_list_to_dist(formulalist, brain_isolen, brain_isodist);

    if (brain_isolen < isolen) {
        for (int i = brain_isolen - offset - 1; i >= 0; i--) {
            isodist[i+offset] = brain_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }
    else {
        for (int i = isolen - offset - 1; i >= 0; i--) {
            isodist[i + offset] = brain_isodist[i];
            if (i < offset) { isodist[i] = 0.0f; }
        }
    }

    free(brain_isodist);
    free(formulalist);

    if (maxval > 0.0f) {
        for (int i = 0; i < isolen; i++) {
            isodist[i] /= maxval;
        }
    }

    return maxval;
}