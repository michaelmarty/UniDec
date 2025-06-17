#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "fftw3.h"
#include "isogendep.h"
#include "isogen_rnaveragine_model32.h"
#include "isogen_rnaveragine_model64.h"
#include "isogen_rnaveragine_model128.h"
#include "isogenrna_model_64.h"



const char rnaOrder[] = "ACGU";
const float rnaveragine_mass = 321.29163925;
const double rnaveragine_comp_numerical[] = {10.70, 9.50, 6.09, 3.70, 0, 0, 0, 0, 0, 0, 0};


const int rna_vectors[][11] = {
    {12,10,5,5,0,0,0,0,0,0,0},
    {12,9,6,3,0,0,0,0,0,0,0},
    {12,10,6,5,0,0,0,0,0,0,0},
    {11,9,7,2,0,0,0,0,0,0,0}
};



void rna_mass_to_list(float initialMass, int* fftlist)
//Mass -> List 11 length list of number of {H, C, N, O, S, Fe, K, Ca, Ni, Zn, Mg}
{
    float rnaMass = rnaveragine_mass;
    float valuePer = initialMass / rnaMass;

    for (int i = 0; i < 5; i++)
    {
        fftlist[i] = (int)round(rnaveragine_comp_numerical[i] * valuePer);
    }

    for (int i = 5;i<11;i++) {
        fftlist[i] =  0;
    }
}

// float* fft_rna_seq_to_dist(char* sequence)
// //String to isodist
// {
//     float* isolist = rnaToVector(sequence);
//     float mass = rnaVectorToMass(isolist);
//     float* rnaCountList = rna_mass_to_list(mass);
//     return fft_rna_list_to_dist(rnaCountList);
//
// }



float fft_rna_mass_to_dist(float mass, float* isodist, int isolen, int offset)
// Mass to dist
{
    int* fftlist = (int*)calloc(11, sizeof(int));
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


void rna_seq_to_vector(const char* seq, float* vector)
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
}


void rna_seq_to_fftlist(const char* sequence, int* fftlist)
{
    // Initialize the formulalist to zero but add the elements of water for the terminii
    fftlist[0] = 0; // Hydrogen
    fftlist[1] = 0; // Carbon
    fftlist[2] = 0; // Nitrogen
    fftlist[3] = 0; // Oxygen
    fftlist[4] = 0; // Sulfur
    fftlist[5] = 0; // Iron
    fftlist[6] = 0; // Potassium
    fftlist[7] = 0; // Calcium
    fftlist[8] = 0; // Nickel
    fftlist[9] = 0; // Zinc
    fftlist[10] = 0; // Magnesium

    int seq_len = strlen(sequence);
    for (int i = 0; i < seq_len; i++) {
        int nt_index = nt_to_index(sequence[i]);
        if (nt_index == -1) {
            printf("Unexpected nucleotide in sequence:%c\n", sequence[i]);
        }
        else {
            for (int j = 0;j<11;j++) {
                fftlist[j] += rna_vectors[nt_index][j];
            }
        }
    }
}


float fft_rna_seq_to_dist(const char* sequence, float* isodist, const int isolen, const int offset)
{
    int* fftlist = (int*)calloc(11, sizeof(int));
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


float nn_rna_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset) {
    float* vector = (float*)calloc(5, sizeof(float));
    if (vector == NULL) {
        printf("Error: Could not allocate memory for vector\n");
    }

    mass_to_vector(mass, vector);

    struct IsoGenWeights weights = SetupWeights(5, isolen);
    if (isolen == 32){ weights = LoadWeights(weights, isogen_rnaveragine_model32_bin); }
    else if ( isolen == 64 ){ weights = LoadWeights(weights, isogen_rnaveragine_model64_bin); }
    else if ( isolen == 128 ){ weights = LoadWeights(weights, isogen_rnaveragine_model128_bin); }
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


float nn_rna_seq_to_dist(const char* seq, float* isodist, int offset) {
    float* vector = calloc(4, sizeof(float));
    rna_seq_to_vector(seq, vector);


    struct IsoGenWeights weights = SetupWeights(4, 64);
    LoadWeights(weights, isogenrna_model_64_bin);


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

    for (int i = 0; i < 64; i++) {
        isodist[i] /= maxval;
    }

    return maxval;
}