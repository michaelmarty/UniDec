#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fftw3.h"
#include "isogendep.h"



const char rnaOrder[] = "ACGU";
const int ISOLENGTH = 128;
const float rnaveragine_mass = 321.29163925;
const double rnaveragine_comp_numerical[] = {9.50, 10.70, 6.09, 3.70, 1.0};



float* fft_rna_list_to_dist(float* isolist)
{

    // for (int i = 0; i < 5; i++)
    // {
    //     printf("%f\n", isolist[i]);
    // }

    int complen = (ISOLENGTH / 2) + 1;

    fftw_complex* hft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* cft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* nft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* oft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* pft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    fftw_complex* allft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));

    if (hft == NULL || cft == NULL || nft == NULL || oft == NULL || pft == NULL || allft == NULL)
    {
        printf("Error: Memory allocation failed for FFTW arrays\n");
        return NULL;
    }


    setup_ft(1, hft, ISOLENGTH, complen); // Hydrogen
    setup_ft(6, cft, ISOLENGTH, complen); // Carbon
    setup_ft(7, nft, ISOLENGTH, complen); // Nitrogen
    setup_ft(8, oft, ISOLENGTH, complen); // Oxygen
    setup_ft(15, pft, ISOLENGTH, complen); // Phosphorus instead of sulfur


    double numc = isolist[0];
    double numh = isolist[1];
    double numn = isolist[2];
    double numo = isolist[3];
    double nump = isolist[4];

    for (int i = 0; i < complen; i++)
    {
        double real, imag;
        double temp_real, temp_imag;

        complex_power(cft[i], numc, &real, &imag);
        complex_power(hft[i], numh, &temp_real, &temp_imag);
        complex_multiplication(real, imag, temp_real, temp_imag, &real, &imag);

        complex_power(nft[i], numn, &temp_real, &temp_imag);
        complex_multiplication(real, imag, temp_real, temp_imag, &real, &imag);

        complex_power(oft[i], numo, &temp_real, &temp_imag);
        complex_multiplication(real, imag, temp_real, temp_imag, &real, &imag);

        complex_power(pft[i], nump, &temp_real, &temp_imag);
        complex_multiplication(real, imag, temp_real, temp_imag, &real, &imag);

        allft[i][0] = real;
        allft[i][1] = imag;
    }


    double* buffer = (double*)fftw_malloc(ISOLENGTH * sizeof(double));
    fftw_plan plan_irfft = fftw_plan_dft_c2r_1d(ISOLENGTH, allft, buffer, FFTW_ESTIMATE);
    fftw_execute(plan_irfft);
    fftw_destroy_plan(plan_irfft);
    fftw_free(allft);


    double max_val = normalize_isodist(buffer, ISOLENGTH);


    float* isodist = (float*)calloc(ISOLENGTH, sizeof(float));
    for (int i = 0; i < ISOLENGTH; i++)
    {
        isodist[i] = (float)buffer[i];
    }
    fftw_free(hft);
    fftw_free(cft);
    fftw_free(nft);
    fftw_free(oft);
    fftw_free(pft);
    fftw_free(buffer);

    return isodist;
}
float* rnaToVector(const char* rna)
//RNA -> Vector 4 length list of number of A C G U
{
    static float vector[4];
    memset(vector, 0, sizeof(vector));
    int len = strlen(rna);
    for (int i = 0; i < len; i++)
    {
        char *pos = strchr(rnaOrder, rna[i]);
        if (pos)
        {
            int index = pos - rnaOrder;
            vector[index] += 1.0f;
        }
    }
    return vector;
}

float* rna_mass_to_list(float initialMass)
//Mass -> List 5 length list of number of C H N O P
{
    float rnaMass = rnaveragine_mass;
    float valuePer = initialMass / rnaMass;

    float* isolist = (float*)calloc(5, sizeof(float));
    for (int i = 0; i < 5; i++)
    {

        isolist[i] = rnaveragine_comp_numerical[i] * valuePer;

    }
    return isolist;
}



float* fft_rna_seq_to_dist(char* sequence)
//String to isodist
{
    float* isolist = rnaToVector(sequence);
    float mass = rnaVectorToMass(isolist);
    float* rnaCountList = rna_mass_to_list(mass);
    return fft_rna_list_to_dist(rnaCountList);

}

float rnaVectorToMass(float* rnaVector)
//Totals mass from vector
{
    float totalMass = 0;
    float values[] = {329.21f, 305.18f, 345.21f, 306.17f};
    for (int i = 0; i < 4; i++)
    {
        totalMass += (rnaVector[i] * values[i]);
    }

    return totalMass;
}

float* fft_rna_mass_to_dist(float mass)
// Mass to dist
{
    float* masses = rna_mass_to_list(mass);
    return fft_rna_list_to_dist(masses);
}





