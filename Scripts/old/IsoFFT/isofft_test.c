//
// Created by marty on 12/31/2024.
//
#include "isofft.h"
#include <stdio.h>
// #include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

int main(const int argc, char *argv[])
{
    printf("Testing IsoFFT with ");
    float mass = 0.0f;
    char *input = NULL;
    int mode = 2;
    if (argc < 2) {
        mass = 20000.f;
        // input = "PdO3H5NClCO4H12C";
        input = "CRRKKKKW";
    }
    else {
        if (isdigit(argv[1][0])){
        mass = strtof(argv[1], NULL);
        mode = 0;
        }
        else {
            input = argv[1];
            mode = 2;
            // Check if there are any numbers in the input or any lower case letters
            for (int i = 0; i < strlen(input); i++) {
                if (isdigit(input[i]) || islower(input[i])) {
                    mode = 1;
                    break;
                }
            }
            // Check if length is less than 5
            if (strlen(input) < 5) {
                mode = 1;
            }
            if (mode==2) {
                mode = 3;
                // Check if all are AUTCG
                for (int i = 0; i < strlen(input); i++) {
                    if (input[i] != 'A' && input[i] != 'U' && input[i] != 'T' && input[i] != 'C' && input[i] != 'G') {
                        mode = 2;
                        break;
                    }
                }
            }
        }
    }

    if (mode==0) {
        printf("Mass: %f\n", mass);
        test_averagine(mass);
    }
    else if (mode==1) {
        printf("Formula: %s\n", input);
        test_formula(input);
    }
    else if (mode==2) {
        printf("Peptide Seq: %s\n", input);
        test_peptide_sequence(input);
    }
    else if (mode==3) {
        printf("NT Sequence: %s\n", input);
        test_nt_sequence(input);
    }
    else {
        printf("Invalid input\n");
    }
}