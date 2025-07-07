//
// Created by mm96978 on 7/7/2025.
//

#ifndef UDIO_H
#define UDIO_H

#endif //UDIO_H

void readfile(const char *infile, int lengthmz, float *dataMZ, float *dataInt);
void readfile3(char *infile, int lengthmz, float *array1, float *array2, float *array3);
void readfile3bin(const char *infile, int lengthmz, float *array1, float *array2, float *array3);
void readfile1bin(const char *infile, int lengthmz, float *data);
void writefile1bin(const char *outstring, int lengthmz, const float *data);
void readmfile(const char *infile, int mfilelen, float *testmasses);
int getfilelength(const char *infile);
int getfilelengthbin(const char *infile, int size, int width);