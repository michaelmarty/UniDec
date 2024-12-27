//
// Created by marty on 10/17/2024.
//

#include <stdio.h>
#include <cstring>
#include <stdlib.h>

#include "isogenmass.h"

void print_isodist(const float *isodist, const int size) {
    // Make a command line printout of the spectrum
    for (int i=10; i>=0; i--) {
        for (int j=0; j<size; j++) {
            const float val = isodist[j];
            if ((float) i < val*10) {
                printf("*");
            }
            else {
                printf(" ");
            }
        }
        printf("\n");
    }

    // Print floats
    // for (int i=0; i<size; i++) {
    //     printf("%.2f ", isodist[i]);
    // }
    printf("\n");
}

int main(const int argc, char *argv[])
{
    printf("IsoGen...\n");
    // Get arguments from command line
    if (argc < 2)
    {
        printf("Usage: %s -mass <mass float>\n", argv[0]);
        return 1;
    }

    if (argc == 3) {
        // Check if argv[1] is -mass
        if (strcmp(argv[1], "-mass") == 0) {
            // Run isogenmass
            float mass = 0.0f;
            mass = strtof(argv[2], nullptr);
            float *isodist = (float *) calloc(64, sizeof(float));
            if(isodist == nullptr) {
                printf("Error allocating memory for isodist\n");
                return 1;
            }
            printf("Running isogenmass on: %.2f Da\n", mass);
            isogenmass(mass, isodist);
            print_isodist(isodist, 64);
            free(isodist);
        }
    }
    return 0;
}
