//
// Created by mm96978 on 7/7/2025.
//

#include <stdlib.h>
#include <stdio.h>
#include "udio.h"


//Reads in x y file.
void readfile(const char *infile, const int lengthmz, float *dataMZ, float *dataInt) {
    char x[500] = {0};
    char y[500] = {0};

    FILE * file_ptr = fopen(infile, "r");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }


    for (int i = 0; i < lengthmz; i++) {
        const int match = fscanf(file_ptr, "%s %s", x, y);
        if (match != 0) {
            // dataMZ[i] = atof(x);
            // dataInt[i] = atof(y);
            char * endptr, * endptr2;
            dataMZ[i] = strtof(x, &endptr);
            dataInt[i] = strtof(y, &endptr2);
            if (*endptr != '\0' || *endptr2 != '\0') {
                printf("Error converting string to float: %s %s\n", x, y);
                exit(2);
            }
        }
    }

    fclose(file_ptr);
}

//Reads in x y z file.
void readfile3(char *infile, int lengthmz, float *array1, float *array2, float *array3) {
    char x[500] = {0};
    char y[500] = {0};
    char z[500] = {0};

    FILE * file_ptr = fopen(infile, "r");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    for (int i = 0; i < lengthmz; i++) {
        const int match = fscanf(file_ptr, "%s %s %s", x, y, z);
        if (match != 0) {
            // array1[i] = atof(x);
            // array2[i] = atof(y);
            // array3[i] = atof(z);
            char * endptr, * endptr2, * endptr3;
            array1[i] = strtof(x, &endptr);
            array2[i] = strtof(y, &endptr2);
            array3[i] = strtof(z, &endptr3);
            if (*endptr != '\0' || *endptr2 != '\0' || *endptr3 != '\0') {
                printf("Error converting string to float: %s %s %s\n", x, y, z);
                exit(2);
            }
        }
    }
    fclose(file_ptr);
}

//Reads in x y z file.
void readfile3bin(const char *infile, const int lengthmz, float *array1, float *array2, float *array3) {
    const int l = lengthmz * 3;
    float * data = (float *) calloc(l, sizeof(float));

    FILE * file_ptr = fopen(infile, "rb");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    if (data) {
        fread(data, sizeof(float), l, file_ptr);
        for (int i = 0; i < lengthmz; i++) {
            array1[i] = data[i * 3];
            array2[i] = data[i * 3 + 1];
            array3[i] = data[i * 3 + 2];
            //printf("%f %f %f\n", array1[i], array2[i], array3[i]);
        }
        free(data);
    }
    fclose(file_ptr);
}

//Reads in list binary file.
void readfile1bin(const char *infile, const int lengthmz, float *data) {
    FILE *file_ptr = fopen(infile, "rb");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }
    fread(data, sizeof(float), lengthmz, file_ptr);
    fclose(file_ptr);
}

//Writes to binary file
void writefile1bin(const char *outstring, const int lengthmz, const float *data) {
    FILE *out_ptr = fopen(outstring, "wb");
    if (out_ptr == 0) {
        printf("Error Opening %s \n", outstring);
        exit(1);
    }
    fwrite(data, sizeof(float), lengthmz, out_ptr);
    fclose(out_ptr);
}


//Reads in single list of values
void readmfile(const char *infile, const int mfilelen, float *testmasses) {
    char input[500] = {0};

    FILE * file_ptr = fopen(infile, "r");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    for (int i = 0; i < mfilelen; i++) {
        const int match = fscanf(file_ptr, "%s", input);
        if (match != 0) {
            // testmasses[i] = atof(input);
            char * endptr;
            testmasses[i] = strtof(input, &endptr);
            if (*endptr != '\0') {
                printf("Error converting string to float: %s\n", input);
                exit(2);
            }
        }
    }

    fclose(file_ptr);
}

// count the number of lines we have in datafile
int getfilelength(const char *infile) {
    int l = 0;
    FILE * file_ptr = fopen(infile, "r");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    int read = -2;
    while (read != EOF) {
        char input[501];
        read = fscanf(file_ptr, "%500[^\n]%*c", input);
        if (read != 0 && read != EOF) {
            l += 1;
        }
    }

    fclose(file_ptr);
    return l;
}

// count the number of lines we have in datafile
int getfilelengthbin(const char *infile, const int size, const int width) {
    int l = 0;
    FILE * file_ptr = fopen(infile, "rb");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    fseek(file_ptr, 0, SEEK_END);
    l = ftell(file_ptr);
    l /= (size * width);

    fclose(file_ptr);
    return l;
}


//Reads in kernel file.
void readkernel(const char *infile, const int lengthmz, double *datax, double *datay) {
    char x[500] = {0};
    char y[500] = {0};

    FILE *file_ptr = fopen(infile, "r");

    if (file_ptr == 0) {
        printf("Could not open %s file\n", infile);
        exit(10);
    }

    for (int i = 0; i < lengthmz; i++) {
        const int match = fscanf(file_ptr, "%s %s", x, y);
        if (match != 0) {
            // datax[i] = atof(x);
            // datay[i] = atof(y);
            char * endptr, * endptr2;
            datax[i] = strtod(x, &endptr);
            datay[i] = strtod(y, &endptr2);
        }
    }

    fclose(file_ptr);
}