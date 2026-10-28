//#include <iostream>
// ReSharper disable CppDeprecatedEntity
#include "isodeclib.h"
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>

#define CHECK(x) do { if (!(x)) { std::fprintf(stderr, "failed line %d: %s\n", __LINE__, #x); return 1; } } while (0)

static int test_mass_group()
{
    double masses[] = {1000.};
    double theory[] = {1000., 2.};
    double observed[] = {999., 0., 1000., 2.};
    MassGroupPeak peaks[2] = {};
    for (int i = 0; i < 2; i++) {
        peaks[i].monoiso = peaks[i].avgmass = 1000.;
        peaks[i].mz = 501.;
        peaks[i].peakint = peaks[i].matchedintensity = 2.;
        peaks[i].charge = 2;
        peaks[i].scan = i;
        peaks[i].monoisos = masses;
        peaks[i].monoisos_count = 1;
        peaks[i].massdist = theory;
        peaks[i].massdist_count = 1;
        peaks[i].decon_centroids = observed;
        peaks[i].centroid_count = 2;
    }
    MassGroupResult* groups = nullptr;
    int count = -1;
    int* ids = nullptr;
    for (int repeat = 0; repeat < 100; repeat++) {
        CHECK(group_mass_peaks_batch(peaks, 2, 5., 3, 1.0033, .99, 100,
                                     repeat % 2, &groups, &count, &ids) == 0);
        CHECK(count == 1 && ids[0] == ids[1]);
        CHECK(groups[0].totalintensity == 4. && groups[0].totalpeaks == 2);
        CHECK(std::isfinite(groups[0].decon_centroids[0]));
        free_mass_groups(groups, count, ids);
    }
    CHECK(group_mass_peaks_batch(nullptr, 0, 5., 3, 1.0033, .7, 100, 0,
                                 &groups, &count, &ids) == 0);
    CHECK(groups == nullptr && ids == nullptr && count == 0);
    peaks[0].massdist_count = INT_MAX;
    CHECK(group_mass_peaks_batch(peaks, 1, 5., 3, 1.0033, .7, 100, 0,
                                 &groups, &count, &ids) == -1);
    CHECK(groups == nullptr && ids == nullptr && count == 0);
    peaks[0].massdist_count = 1;
    peaks[0].scan = INT_MIN;
    peaks[1].scan = INT_MAX;
    CHECK(group_mass_peaks_batch(peaks, 2, 5., 3, 1.0033, .7, 100, 0,
                                 &groups, &count, &ids) == 0);
    CHECK(count == 2);
    free_mass_groups(groups, count, ids);
    std::puts("native grouping ABI checks passed");
    return 0;
}

int main(const int argc, char *argv[])
{
    if (argc == 2 && std::strcmp(argv[1], "--test-mass-group") == 0)
        return test_mass_group();
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

    run(filename, outputfile, weightfile, "Pep");
    // Free memory
    delete[] outputfile;
}

