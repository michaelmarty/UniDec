//#include <iostream>
// ReSharper disable CppDeprecatedEntity
#include "isodeclib.h"
#include <cstdio>
#include <cstring>


int main(const int argc, char *argv[])
{
    printf("Running Exe\n");
    // Get arguments from command line
    if (argc < 2)
    {
        printf("Usage: %s <inputfile> <optional weightfile>\n", argv[0]);
        return 1;
    }
    // Load the file name from the command line
    char *filename = argv[1];
    // Reformat input file to give output file name
    const auto outputfile = new char[strlen(filename) + 5];
    strcpy(outputfile, filename);
    strcat(outputfile, ".out");
    // Print output file
    printf("Output file: %s\n", outputfile);

    char *weightfile=nullptr; // File name for the modification file
    if (argc < 3) {
        // Set up default file paths
        printf("Using Default Weight File\n");
    }
    else {
        // Load the file name from the command line
        weightfile = argv[2];
        printf("Weight file: %s\n", weightfile);
    }

    run(filename, outputfile, weightfile);
    // Free memory
    delete[] outputfile;
}

