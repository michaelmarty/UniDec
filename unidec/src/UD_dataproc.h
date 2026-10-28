
/*
* UD_dataproc.h
*
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//
#pragma once

#ifndef DATAPROC_HEADER
#define DATAPROC_HEADER

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "h5io.h"
#include "udio.h"

int pool1d(float **oldmz, float **oldint, const int oldlen, const float mzbins);
int chop1d(float **oldmz, float **oldint, const int oldlen, const float min, const float max);
int remove_duplicates(float **oldmz, float **oldint, int oldlen);
int remove_middle_zeros(float **oldmz, float **oldint, int oldlen);
void norm1d(float *dataInt, const int lengthmz);
void background_subtract(float *dataInt, const int bsub, const int lengthmz);
int compare(const void* a, const void* b);
int compare_highest_first(const void* b, const void* a);
int data_reduction(float **oldmz, float **oldint, int oldlen, const float redper);
void process_data(int argc, char *argv[], Config config);


#endif