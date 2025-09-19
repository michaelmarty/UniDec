#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include "fftw3.h"
#include "isogendep.h"
#include "isogenpep.h"
#include "isogenpep_model_16.h"
#include "isogenpep_model_64.h"
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
int pep_seq_to_nnvector(const char* seq, float* vector) {
    int length = strlen(seq);

    int aas = 0;

    int in_mod = 0;

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
            if (seq[i] == '['){ in_mod = 0; }
        }
    }
    return aas;
}


float nn_pep_seq_to_dist(const char* seq, float* isodist, int isolen, int offset){
    float* vector = (float *) calloc(20, sizeof(float));
    int aas = pep_seq_to_nnvector(seq, vector);

    if (aas > 300) {
        printf("Sequence contains too many amino acids (>300).");
        return -1.0f;
    }

    struct IsoGenWeights weights;
    float* nn_isodist;
    int nn_isolen;

    if (aas >= 1 && aas <= 50) {
        weights = SetupWeights(20, 16);
        weights = LoadWeights(weights, isogenpep_model_16_bin);
        nn_isolen = 16;
        nn_isodist = (float*)calloc(nn_isolen, sizeof(float));
    }
    else {
        weights = SetupWeights(20, 64);
        weights = LoadWeights(weights, isogenpep_model_64_bin);
        nn_isolen = 64;
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

//fft
void add_mod_to_fftlist(char* mod, int* fftlist) {
    int i = 0;
    while (mod[i] != '\0') {
        if (!isupper(mod[i])) {
            printf("Error while parsing modification formula:%s\n", mod);
            fflush(stdout);
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
                symbol_matched = 1;
                break;
            }
        }
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


int nn_pep_mass_to_isolen(const float mass) {
    if (mass < 1200) {
        return 8;
    }

    if (mass < 11000) {
        return 32;
    }
    if (mass > 11000 && mass < 55000) {
        return 64;
    }

    if (mass < 120000) {
        return 128;
    }
    return -1;
}


int fft_pep_mass_to_isolen(const float mass) {
    if (mass < 50000) {
        return 64;
    }
    if (mass < 120000) {
        return 128;
    }
    return 1024;
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


int get_pep_isolen_from_seq(const char* seq) {
    int aas = 0;

    int length = strlen(seq);
    int inmod = 0;

    for (int i = 0;i<length;i++) {
        if (seq[i] == '[') {
            inmod = 1;
        }
        else if (inmod == 1 && seq[i] == ']') {
            inmod = 0;
        }
        else {
            aas += 1;
        }
    }

    if (aas < 50){return 16;}
    if (aas < 300){return 64;}
    return 128;
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

    //int fft_isolen = get_pep_isolen_from_seq(sequence);
    int fft_isolen = 64;

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


//fft
float fft_pep_mass_to_dist(const float mass, float *isodist, const int isolen, const int offset)
{
    int* fftlist = (int*)calloc(11, sizeof(int));
    pep_mass_to_fftlist(mass, fftlist);

    int fft_isolen = fft_pep_mass_to_isolen(mass);

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


//nn
float nn_pep_mass_to_dist(const float mass, float* isodist, const int isolen, const int offset) {
    float* vector = (float*)calloc(5, sizeof(float));
    if (vector == NULL) {
        printf("Error: Could not allocate memory for vector\n");
    }

    mass_to_vector(mass, vector);

    int nn_isolen = nn_pep_mass_to_isolen(mass);

    if (nn_isolen == -1) {
        printf("Error: Mass outside of allowed NN mass range: %f\n", mass);
        return -1;
    }

    float* nn_isodist = (float*)calloc(nn_isolen, sizeof(float));

    struct IsoGenWeights weights = SetupWeights(5, nn_isolen);
    if (nn_isolen == 8){ weights = LoadWeights(weights, isogenmass_model_8_bin); }
    else if ( nn_isolen == 32 ){ weights = LoadWeights(weights, isogenmass_model_32_bin); }
    else if ( nn_isolen == 64 ){ weights = LoadWeights(weights, isogenmass_model_64_bin); }
    else { weights = LoadWeights(weights, isogenmass_model_128_bin); }

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




#define MAX_LINE_LENGTH 1024
#define INITIAL_CAPACITY 10

typedef struct {
    double mass;
    char* sequence;  // dynamically allocated string
} ProteinEntry;

typedef struct {
    ProteinEntry* entries;
    size_t size;
    size_t capacity;
} ProteinList;

void free_protein_list(ProteinList* list) {
    for (size_t i = 0; i < list->size; ++i) {
        free(list->entries[i].sequence);
    }
    free(list->entries);
    list->entries = NULL;
    list->size = 0;
    list->capacity = 0;
}

ProteinList read_masses_seqs_file(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file!\n");
        ProteinList empty = {NULL, 0, 0};
        return empty;
    }

    ProteinList list;
    list.size = 0;
    list.capacity = INITIAL_CAPACITY;
    list.entries = malloc(list.capacity * sizeof(ProteinEntry));
    if (!list.entries) {
        fclose(file);
        fprintf(stderr, "Memory allocation failed!\n");
        ProteinList empty = {NULL, 0, 0};
        return empty;
    }

    char line[MAX_LINE_LENGTH];
    while (fgets(line, sizeof(line), file)) {
        double mass;
        char seq_buf[MAX_LINE_LENGTH];

        if (sscanf(line, "%lf %s", &mass, seq_buf) != 2) {
            continue;  // skip malformed lines
        }

        if (list.size >= list.capacity) {
            list.capacity *= 2;
            ProteinEntry* temp = realloc(list.entries, list.capacity * sizeof(ProteinEntry));
            if (!temp) {
                fprintf(stderr, "Reallocation failed!\n");
                break;
            }
            list.entries = temp;
        }

        // Allocate and copy the sequence
        list.entries[list.size].mass = mass;
        list.entries[list.size].sequence = malloc(strlen(seq_buf) + 1);
        if (!list.entries[list.size].sequence) {
            fprintf(stderr, "String allocation failed!\n");
            break;
        }
        strcpy(list.entries[list.size].sequence, seq_buf);
        list.size++;
    }

    fclose(file);
    return list;
}



extern void run_file(const char* filename) {
    ProteinList list = read_masses_seqs_file(filename);

    int isolen = 64;

    for (int i = 0; i < list.size; i++) {
        float* isodist = (float*)calloc(isolen, sizeof(float));
        nn_pep_mass_to_dist(list.entries[i].mass, isodist, isolen, 0);
        fft_pep_mass_to_dist(list.entries[i].mass, isodist, isolen, 0);
        nn_pep_mass_to_dist(list.entries[i].mass, isodist, isolen, 0);
        fft_pep_seq_to_dist(list.entries[i].sequence, isodist, isolen, 0);
        nn_pep_seq_to_dist(list.entries[i].sequence, isodist, isolen, 0);
        free(isodist);
    }
    free_protein_list(&list);
    return;
}
