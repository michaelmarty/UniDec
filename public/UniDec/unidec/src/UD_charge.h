/*
* UD_charge.h
*
*  Created on : 19 Dec 2017
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//

#ifndef UD_CHARGE_H
#define UD_CHARGE_H

#include "udcore.h"
#include "udtools.h"
#include "udio.h"
#include "udstruct.h"

float extract_zcom(Config config, const float peak, const float *xvals, const float *yvals, const float *zvals, const int mlen, const int zlen, const float thresh);
float extract_zmax(Config config, const float peak, const float *xvals, const float *yvals, const float *zvals, const int mlen, const int zlen, const float thresh);
float charge_extract_switch(const Config config, const float peak, const float *xvals, const float *yvals, const float *zvals, const int mlen, const int zlen);
void charge_peak_extracts(int argc, char *argv[], Config config, const int ultra);
void z_slice(const float *ivals, const int mindex, float *zdat, const int zlen, const float thresh);
void z_slice_range(const float *ivals, const int mindex1, const int mindex2, float *zdat, const int zlen, const float thresh);

#endif