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
const float rnaveragine_mass = 320.28;
const double rnaveragine_comp_numerical[] = {9.50, 10.75, 3.75, 7.0, 0};

const int rna_vectors[][5] = {
    {10,11,5,2,0},
    {9,11,3,3,0},
    {10,11,5,3,0},
    {9,10,2,4,0}
};


void rna_mass_to_list(float initialMass, int* fftlist)
//Mass -> List 5 length list of number of {C, H, N, O, S}
{
    float rnaMass = rnaveragine_mass;
    //Addition of an extra phospho here simplifies the monomer number calculation
    float valuePer = (initialMass + 95.9534) / rnaMass;

    for (int i = 0; i < 5; i++)
    {
        fftlist[i] = (int)round(rnaveragine_comp_numerical[i] * valuePer);
    }

    //Correct for additional one less phospho than monomers
    fftlist[3] -= 4;
}



float fft_rna_mass_to_dist(float mass, float* isodist, int isolen, int offset)
// Mass to dist
{
    int* fftlist = (int*)calloc(5, sizeof(int));
    rna_mass_to_list(mass, fftlist);

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
    free(fftlist);
    return max_val;
}


int nt_to_index(char nt) {
    for (int i = 0;i<4;i++) {
        if (rnaOrder[i] == nt) {
            return i;
        }
    }
    return -1;
}


int rna_seq_to_vector(const char* seq, float* vector)
//RNA -> Vector 4 length list of number of A C G U
{
    int len = strlen(seq);
    for (int i = 0; i < len; i++)
    {
        int nt_index = nt_to_index(seq[i]);
        if (nt_index == -1) {
            printf("Unexpected nucleotide in sequence: %c\n", seq[i]);
        }
        else {
            vector[nt_index] += 1.0f;
        }
    }
    return len;
}


void rna_seq_to_fftlist(const char* sequence, int* fftlist)
{
    // Initialize the formulalist to zero but add the elements of water for the terminii
    fftlist[0] = 0; // Carbon
    fftlist[1] = 0; // Hydrogen
    fftlist[2] = 0; // Nitrogen
    fftlist[3] = 0; // Oxygen
    fftlist[4] = 0; // Sulfur


    int seq_len = strlen(sequence);
    for (int i = 0; i < seq_len; i++) {
        int nt_index = nt_to_index(sequence[i]);
        if (nt_index == -1) {
            printf("Unexpected nucleotide in sequence:%c\n", sequence[i]);
        }
        else {
            for (int j = 0;j<5;j++) {
                fftlist[j] += rna_vectors[nt_index][j];
            }
        }
    }

    fftlist[3] += (seq_len - 1) * 4;
}


float fft_rna_seq_to_dist(const char* sequence, float* isodist, const int isolen, const int offset)
{
    int* fftlist = (int*)calloc(5, sizeof(int));
    // Check for null
    if (fftlist == NULL)
    {
        printf("Error: Could not allocate memory for formulalist\n");
        return 1;
    }
    rna_seq_to_fftlist(sequence, fftlist);

    float maxval = fft_list_to_dist(fftlist, isolen, isodist);
    free(fftlist);

    for (int i = isolen - offset - 1; i >= 0; i--)
    {
        isodist[i + offset] = isodist[i];
        if (i < offset) { isodist[i] = 0.0f; }
    }

    for (int i = 0;i<isolen;i++) {
        isodist[i] /= maxval;
    }

    return maxval;
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


float nn_rna_seq_to_dist(const char* seq, float* isodist, int isolen, int offset) {
    float* vector = calloc(4, sizeof(float));
    if (vector == NULL) {
        printf("Error: Could not allocate memory for vector\n");
    }

    int len = rna_seq_to_vector(seq, vector);

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
    if (len > 500) {
        printf("Error: Sequence length outside of allowed NN range: %i\n", len);
        return -1;
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