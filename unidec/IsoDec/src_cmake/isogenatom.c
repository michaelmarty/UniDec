#include <ctype.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "fftw3.h"
#include "isogendep.h"

float* get_formula_dist(char* sequence)
{
    int isolen = 128;

    int* formulalist = (int*)calloc(isolen, sizeof(int));


    formula_to_vector(sequence, formulalist);


    float* isodist = (float*)calloc(isolen, sizeof(float));
    if (isodist == NULL)
    {
        printf("Error: Could not allocate memory for isodist\n");
        free(formulalist);
        return NULL;
    }

    int len = (int)strlen(sequence);
    printf("Length of sequence: %d\n", len);
    fflush(stdout);
    isodist_from_formula_vector(formulalist, len, isodist, isolen);

    // Free temporary memory
    free(formulalist);

    return isodist;
}


float isodist_from_formula_vector(const int* formula, const int formulalen, float* isodist, const int isolen)
{


    printf("\n");
    fflush(stdout);
    // Create FFTW plans for each and fft objects for each
    int const complen = (int)(isolen / 2) + 1;

    fftw_complex* ftbuffer = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    // Check for null pointers
    if (ftbuffer == NULL)
    {
        printf("Error: Could not allocate memory for fftw_complex objects\n");
        return 1;
    }

    fftw_complex* allft = (fftw_complex*)fftw_malloc(complen * sizeof(fftw_complex));
    double* buffer = (double*)fftw_malloc(isolen * sizeof(double));
    // Check for null pointers
    if (allft == NULL || buffer == NULL)
    {
        printf("Error: Could not allocate memory for fftw_complex objects\n");
        return 1;
    }
    int firstelement = 1;
    for (int i = 0; i < formulalen; i++)
    {
        if (formula[i] > 0)

        if (formula[i] == 0) { continue; }
        if (i > numelements)
        {
            printf("Formula vector is too long\n");
            break;
        }
        setup_ft(i + 1, ftbuffer, isolen, complen);

        const int numc = formula[i];


        if (firstelement)
        {
            for (int j = 0; j < complen; j++)
            {
                double coutr = 0;
                double couti = 0;

                complex_power(ftbuffer[j], numc, &coutr, &couti);

                allft[j][0] = coutr;
                allft[j][1] = couti;
            }
            firstelement = 0;
        }
        else
        {
            for (int j = 0; j < complen; j++)
            {
                double coutr = 0;
                double couti = 0;

                complex_power(ftbuffer[j], numc, &coutr, &couti);

                double finalr = allft[j][0];
                double finali = allft[j][1];

                complex_multiplication(finalr, finali, coutr, couti, &finalr, &finali);

                allft[j][0] = finalr;
                allft[j][1] = finali;
            }
        }
    }

    fftw_free(ftbuffer);

    fftw_plan plan_irfft = fftw_plan_dft_c2r_1d(isolen, allft, buffer, FFTW_ESTIMATE);
    fftw_execute(plan_irfft);
    fftw_destroy_plan(plan_irfft);
    fftw_free(allft);

    const double max_val = normalize_isodist(buffer, isolen);

    fflush(stdout);
    //memcpy buffer to isodist
    for (int i = 0; i < isolen; i++)
    {
        isodist[i] = (float)buffer[i];
    }
    fftw_free(buffer);
    return (float)max_val;
}

void formula_to_vector(const char* formula, int *formulalist)
{
    // Initialize the formulalist to zero
    for (int i = 0; i < numelements; i++)
    {
        formulalist[i] = 0;
    }

    const int len = (int)strlen(formula);
    for (int i = 0; i < len;)
    {
        // Find the element symbol
        char element[3] = {0};
        element[0] = formula[i++];
        if (i < len && islower(formula[i]))
        {
            element[1] = formula[i++];
        }

        // Find the element index in the element_names array
        int element_index = -1;
        for (int j = 0; j < numelements; j++)
        {
            if (strcmp(element, element_names[j]) == 0)
            {
                element_index = j;
                break;
            }
        }

        // If the element is not found, continue to the next character
        if (element_index == -1)
        {
            continue;
        }

        // Find the number of atoms of this element
        int count = 0;
        while (i < len && isdigit(formula[i]))
        {
            count = count * 10 + (formula[i++] - '0');
        }
        if (count == 0)
        {
            count = 1;
        }

        // Update the formulalist
        formulalist[element_index] += count;
    }
    // printf("Processed formula: ");
    // // Print forumla vector
    // for (int i = 0; i < numelements; i++)
    // {
    //     if (formulalist[i] == 0) { continue; }
    //     printf("{%s %d} ", element_names[i], formulalist[i]);
    // }
    // printf("\n");
}

void fft_atom_sequence_to_dist(const char* sequence, float* isodist, int isolen)
{
    // Convert the sequence to a formula vector
    int* formulalist = (int*)calloc(109, sizeof(int));
    if (formulalist == NULL)
    {
        printf("Error: Could not allocate memory for formulalist\n");
        return;
    }

    formula_to_vector(sequence, formulalist);

    // Calculate the isodist from the formula vector
    isodist_from_formula_vector(formulalist, 109, isodist, isolen);

    // Free temporary memory
    free(formulalist);

}



