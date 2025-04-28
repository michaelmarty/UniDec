#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "fftw3.h"
#include "isogendep.h"
#include "isogenpep.h"



double massavgine = 111.1254;
double avgine[5] = {4.9384, 7.7583, 1.3577, 1.4773, 0.0417};
int numaminoacids = 23;


#define ISO_LEN 32

float isogen_fancy(const float mass, float *isodist, const int isolen, const int offset, char* type) {
    // Run Isogen mass
    isogen_nn(mass, isodist, isolen, type);
    // Return the max value
    float max = 0.0f;
    for (int i = isolen-offset-1; i >= 0 ; i--) {
        isodist[i+offset] = isodist[i];
        if (isodist[i] > max) {
            max = isodist[i];
        }
        if (i<=offset){isodist[i]=0.0f;}
    }
    return max;
}

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

float isomike(const float mass, float * isodist, const int isolen, const int offset) {
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
void formula_vector_from_peptide_sequence(const char* sequence, int* vector)
{
    // Initialize the formulalist to zero but add the elements of water for the terminii
    vector[0] = 2; // Hydrogen
    vector[1] = 0; // Carbon
    vector[2] = 0; // Nitrogen
    vector[3] = 1; // Oxygen
    vector[4] = 0; // Sulfur

    for (int i = 0; i < strlen(sequence); i++)
    {
        char aa = sequence[i];
        int aa_index = -1;
        for (int j = 0; j < numaminoacids; j++)
        {
            if (aa == amino_acid_names[j][0])
            {
                aa_index = j;
                break;
            }
        }
        if (aa_index == -1)
        {
            continue;
        }
        for (int j = 0; j < num_simp_elements; j++)
        {
            vector[j] += amino_acid_vectors[aa_index][j];
        }
    }



    for (int i = 0; i < num_simp_elements; i++)
    {
        if (vector[i] == 0) { continue; }

    }
}

//fft
float fft_pep_seq_to_dist(const char* sequence, float* isodist, const int isolen)
{
    int* formulalist = (int*)calloc(numelements, sizeof(int));
    // Check for null
    if (formulalist == NULL)
    {
        printf("Error: Could not allocate memory for formulalist\n");
        return 1;
    }
    formula_vector_from_peptide_sequence(sequence, formulalist);
    float maxval = fft_pep_list_to_dist(formulalist, isolen, isodist);
    free(formulalist);
    return maxval;
}

//fft
void mass_to_formula_averaging(const float mass, int* forumla)
{
    for (int i = 0; i < 5; i++)
    {
        forumla[i] = (int)round((mass / massavgine) * avgine[i]);
    }
}
//fft
float fft_pep_list_to_dist(const int isolist[5], const int length, float* isodist)
{
    // Isolist as {hydrogen, carbon, nitrogen, oxygen, sulfur}
    int const complen = (int)(length / 2) + 1;

    fftw_complex* hft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* cft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* nft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* oft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* sft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    // Check for null pointers
    if (hft == NULL || cft == NULL || nft == NULL || oft == NULL || sft == NULL)
    {
        printf("Error: Could not allocate memory for fftw_complex objects\n");
        return 1;
    }
    setup_ft(1, hft, length, complen);
    setup_ft(6, cft, length, complen);
    setup_ft(7, nft, length, complen);
    setup_ft(8, oft, length, complen);
    setup_ft(16, sft, length, complen);

    fftw_complex* allft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    double* buffer = (double*)fftw_malloc(length * sizeof(double));
    // Check for null pointers
    if (allft == NULL || buffer == NULL)
    {
        printf("Error: Could not allocate memory for fftw_complex objects\n");
        return 1;
    }
    //ensure normalizing

    convolve_all(isolist, allft, cft, hft, nft, oft, sft, complen);
    fftw_free(hft);
    fftw_free(cft);
    fftw_free(nft);
    fftw_free(oft);
    fftw_free(sft);

    fftw_plan plan_irfft = fftw_plan_dft_c2r_1d(length, allft, buffer, FFTW_ESTIMATE);
    fftw_execute(plan_irfft);
    fftw_destroy_plan(plan_irfft);
    fftw_free(allft);

    const double max_val = normalize_isodist(buffer, length);

    //memcpy buffer to isodist
    for (int i = 0; i < length; i++)
    {
        isodist[i] = (float)buffer[i];
    }

    fftw_free(buffer);
    return (float)max_val;
}


float fft_pep_mass_to_dist(const float mass, float *isodist, const int isolen, const int offset)
{
    int formulalist[5];
    mass_to_formula_averaging(mass, formulalist);
    float max_val = fft_pep_list_to_dist(formulalist, isolen, isodist);

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


void formula_vector_from_nt_sequence(const char* sequence, int* vector)
{
    // Initialize the formulalist to zero but add the elements of water for the terminii
    vector[0] = 2; // Hydrogen
    vector[1] = 0; // Carbon
    vector[2] = 0; // Nitrogen
    vector[3] = 1; // Oxygen
    vector[4] = 0; // Sulfur

    int mode = 1;
    // If U is in the sequence, it's RNA
    for (int i = 0; i < strlen(sequence); i++)
    {
        if (sequence[i] == 'U')
        {
            mode = 0;
            break;
        }
    }

    // mode 0 is RNA, mode 1 is DNA
    char** names;
    const int (*vectors)[5];
    if (mode == 0)
    {
        names = rna_names;
        vectors = rna_vectors;
    }
    else
    {
        names = rna_names;
        vectors = dna_vectors;
    }

    for (int i = 0; i < strlen(sequence); i++)
    {
        const char nt = sequence[i];
        int nt_index = -1;
        for (int j = 0; j < 4; j++)
        {
            if (nt == names[j][0])
            {
                nt_index = j;
                break;
            }
        }
        if (nt_index == -1)
        {
            continue;
        }
        for (int j = 0; j < num_simp_elements; j++)
        {
            vector[j] += vectors[nt_index][j];
        }
    }
    if (mode == 0) { printf("Processed RNA sequence: "); }
    else { printf("Processed DNA sequence: "); }

    // Print forumla vector

    for (int i = 0; i < num_simp_elements; i++)
    {
        if (vector[i] == 0) { continue; }

    }
}


float isodist_from_nt_sequence(const char* sequence, float* isodist, const int isolen)
{
    int* formulalist = (int*)calloc(numelements, sizeof(int));
    // Check for null
    if (formulalist == NULL)
    {
        printf("Error: Could not allocate memory for formulalist\n");
        return 1;
    }
    formula_vector_from_nt_sequence(sequence, formulalist);
    float maxval = fft_pep_list_to_dist(formulalist, isolen, isodist);
    free(formulalist);
    return maxval;
}

int get_element_index(const char* symbol){
    for (int i = 0; i < 109; i++)
    {
        if (strcmp(elements[i], symbol) == 0)
        {
            return i;
        }
    }
    return -1;
}
float* formulaToVector(const char* formula)
{
    float* vector = malloc(109 * sizeof(float));
    if (!vector) {
        printf("Memory allocation failed!\n");
        exit(1);
    }
    memset(vector, 0, 109 * sizeof(float));
    int len = strlen(formula);
    for (int i = 0; i < len;){
        char element[3] = {0};
        int count = 0;
        if (i+1 < len && islower(formula[i+1])){
            element[0] = formula[i];
            element[1] = formula[i+1];
            element[2] = '\0';
            i += 2;
        } else{
            element[0] = formula[i];
            element[1] = '\0';
            i++;
        }
        while(i < len && isdigit(formula[i])){
            count = count * 10 + formula[i] - '0';
            i++;
        }
        if (count == 0){
            count = 1;
        }
        int index = get_element_index(element);
        if (index != -1){
            vector[index] += count;
        }
        else {
            printf("Error: Unknown element %s in formula %s\n", element, formula);
            exit(1);
        }
        vector[index]+=count;
    }
    return vector;
}

float isodist_from_averagine_mass(const float mass, float* isodist, const int isolen, const int offset)
{

    int formulalist[5];
    mass_to_formula_averaging(mass, formulalist);
    const float max_val = fft_pep_list_to_dist(formulalist, isolen, isodist);
    for (int i = isolen - offset - 1; i >= 0; i--)
    {
        isodist[i + offset] = isodist[i];
        if (i < offset) { isodist[i] = 0.0f; }
    }
    return max_val;
}




