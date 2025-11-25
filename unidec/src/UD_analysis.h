/*
* UD_analysis.h
*
*  Created on : 3 June 2017
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//



#ifndef ANALYSIS_HEADER
#define ANALYSIS_HEADER

#include <time.h>

#include "UD_dataproc.h"
#include "udio.h"
#include "UD_score.h"

void interpolate_merge(const float *massaxis, float *outint, const float *tempaxis, const float *tempint, int mlen, int templen);

Config get_global_min_max(int argc, char* argv[], Config config, const char* dtype);
void make_grid(int argc, char *argv[], Config config, const char *dtype, const char *out1, const char *out2, const char *out3);

int is_peak(const float *dataMZ, const float *dataInt, int lengthmz, float window, float thresh, int index);
int peak_detect(const float *dataMZ, const float *dataInt, int lengthmz, float window, float thresh, float *peakx, float *peaky);
void peak_norm(float *peaky, int plen, int peaknorm);

float extract_height(Config config, float peak, const float *xvals, const float *yvals, int length);
float extract_localmax(Config config, float peak, const float *xvals, const float *yvals, int length);
float extract_localmax_position(Config config, float peak, const float *xvals, const float *yvals, int length);
float extract_integral(Config config, float peak, const float *xvals, const float *yvals, int length, float thresh);
float extract_center_of_mass(Config config, float peak, const float *xvals, const float *yvals, int length, float thresh);
float extract_estimated_area(Config config, float peak, const float* xvals, const float* yvals, int length);
float extract_switch(Config config, float peak, const float* xvals, const float* yvals, int length);

void peak_extracts(Config config, const float *peakx, hid_t file_id, const char *dtype, int plen, int ultra);

void get_all_peaks(int argc, char *argv[], Config config);
void get_peaks(int argc, char *argv[], Config config, int ultra);


#endif