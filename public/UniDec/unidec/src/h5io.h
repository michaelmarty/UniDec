//
// Created by mm96978 on 7/8/2025.
//

#ifndef H5IO_H
#define H5IO_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "hdf5.h"
#include "hdf5_hl.h"
#include "udstruct.h"

int check_group(const hid_t file_id, char *dataname);
int check_group_noexit(hid_t file_id, char* dataname, const int silent);
void delete_group(const hid_t file_id, const char* dataname);

int mh5getfilelength(const hid_t file_id, char *dataname);
void mh5readfile2dcolumn(const hid_t file_id, char* dataname, float* outdata, const int col);
void mh5readfile2d(const hid_t file_id, char *dataname, const int lengthmz, float *dataMZ, float *dataInt);
void mh5readfile3d(const hid_t file_id, char *dataname, const int lengthmz, float *dataMZ, float *dataInt, float *data3);
void mh5readfile1d(const hid_t file_id, char *dataname, float *data);

void mh5writefile1d(const hid_t file_id, const char *dataname, const int length, const float *data1);
void mh5writefile2d(const hid_t file_id, const char *dataname, const int length, const float *data1, const float *data2);
void mh5writefile2d_grid(const hid_t file_id, const char *dataname, const int length1, const int length2, const float *data1);

int int_attr(hid_t file_id, const char *gname, char * attr_name, int value);
float float_attr(hid_t file_id, const char *gname, char * attr_name, float value);
void write_attr_int(hid_t file_id, const char *gname, const char *attr_name, int value);
void write_attr_float(hid_t file_id, const char *gname, const char *attr_name, float value);

Config mh5LoadConfig(Config config, hid_t file_id);

void makegroup(const hid_t file_id, const char *dataset);
void strjoin(const char *s1, const char *s2, char *newstring);

void set_needs_grids(const hid_t file_id);
void set_got_grids(const hid_t file_id);
int question_grids(const hid_t file_id);
int question_rawgrid(hid_t file_id);

#endif //H5IO_H
