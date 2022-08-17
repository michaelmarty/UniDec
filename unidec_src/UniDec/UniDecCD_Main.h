/*
* UniDecCD_Main.h
*
*  Created on : 16 August 2022 (at home with covid :( )
* Author : Michael.Marty
*/

//
// 
// Copyright 2022 University of Arizona
//
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
//#include <fftw3.h>
//#include "UniDecIM.h"


void blur_it_CD(float * output, const float * input, const int* upinds, const int *loinds, const int length, const float floor)
{
	//#pragma omp parallel for private (i), schedule(dynamic)
	for (int i = 0; i < length; i++) {
		int lowindex = loinds[i];
		int uppindex = upinds[i];
		float i1 = input[i];
		float i2 = input[lowindex];
		float i3 = input[uppindex];
		if (floor > 0) {
			i1 = logf(i1 + floor);
			i2 = logf(i2 + floor);
			i3 = logf(i3 + floor);
			float newval = (i1 + i2 + i3) / 3;
			newval = expf(newval) - floor;
			if (newval < 0) {
				newval = 0;
			}
			output[i] = newval;
		}
		else {
			float ratio = fabs(floor);
			output[i] = (i1 + i2 * ratio + i3 * ratio) / 3;
		}		
	}
}


int run_unidec_CD(int argc, char* argv[], Config config) {

	time_t starttime, endtime;
	starttime = time(NULL);

	printf("Opening File: %s\n", config.infile);
	int lines=0;
	lines = getfilelength(config.infile);

	printf("Length of data: %d\n", lines);

	int size[3] = { 0,0,0 };

	int* zupind = NULL, * zloind = NULL;
	char* barr = NULL;

	float* mzdat = NULL, * mzext = NULL,
		* zdat = NULL, * zext = NULL,
		* dataInt = NULL,
		* peakshape = NULL,
		*blur = NULL, *newblur=NULL, *newblur2=NULL
		;
	
	//Reading In Data
	mzdat = calloc(lines, sizeof(float));
	zdat = calloc(lines, sizeof(float));
	dataInt = calloc(lines, sizeof(float));
	readfile3(config.infile, lines, mzdat, zdat, dataInt);

	//Determine the maximum intensity in the data
	float dmax = Max(dataInt, lines);
	//float betafactor = 1;
	//if (dmax > 1) { betafactor = dmax; }


	//Getting Dimension Sizes of Data
	size[0] = GetSize0(mzdat, lines);
	size[1] = GetSize1(zdat, lines);
	int totlen = size[0] * size[1] ;
	size[2] = totlen;
	printf("Dimensions of data: %d mz by %d z: %d total\n", size[0], size[1], totlen);
	if (totlen > 155E6) { printf("Warning: May exceed system memory capacity\n"); }

	//Extracting mz and z ranges
	mzext = calloc(size[0], sizeof(float));
	zext = calloc(size[1], sizeof(float));
	Extract(mzext, zext, mzdat, zdat, size);
	float mzranges[4] = { 0,0,0,0 };
	mzranges[0] = mzext[0];
	mzranges[1] = mzext[size[0] - 1];
	mzranges[2] = zext[0];
	mzranges[3] = zext[size[1] - 1];
	printf("MZ Range: %f to %f\n", mzranges[0], mzranges[1]);
	printf("Z Range: %f to %f\n", mzranges[2], mzranges[3]);

	//Make Peak Shape Kernel
	peakshape = calloc(lines, sizeof(float));
	GetPeaks(peakshape, size, mzext, zext, config.mzsig, config.csig, config.psfun);
	printf("Peak Shape Set\n");


	//Set up blur
	zupind = calloc(lines, sizeof(int));
	zloind = calloc(lines, sizeof(int));
	//#pragma omp parallel for private (i), schedule(dynamic)
	for (int i = 0; i < lines; i++) {
		float mz = mzdat[i];
		float z = zdat[i];
		int zind = (int) z - mzranges[2];
		//int mzind = (int)roundf((mz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (size[0] - 1));//i / size[1];
		float mass = calcMass(mz, z, config.adductmass);
		
		float uppermz = calcMz(mass, z+1, config.adductmass);

		float lowermz = 0;
		if ((z - 1) != 0) {
			lowermz = calcMz(mass, z - 1, config.adductmass);
		}

		int mzupind = (int) roundf((uppermz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (size[0] - 1));
		if (z == mzranges[3] || mzupind < 0 || mzupind > size[0]) {
			zupind[i] = i;
		}
		else
		{
			int newind = index2D(size[1], mzupind, zind + 1);
			if (newind < 0) { newind = 0; }
			if (newind >= lines) { newind = lines-1; }
			zupind[i] = newind;
		}

		int mzloind = (int) roundf((lowermz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (size[0] - 1));
		if (zind == 0 || mzloind < 0 || mzloind > size[0]) {
			zloind[i] = i;
		}
		else
		{
			int newind = index2D(size[1], mzloind, zind-1);
			if (newind < 0) { newind = 0; }
			if (newind >= lines) { newind = lines - 1; }
			zloind[i] = newind;
		}		
	}

	blur = calloc(lines, sizeof(float));
	newblur = calloc(lines, sizeof(float));
	newblur2 = calloc(lines, sizeof(float));
	barr = calloc(lines, sizeof(char));
	for (int i = 0; i < lines; i++) { barr[i] = 1; }
	rsize_t matsize = lines * sizeof(float);
	memcpy_s(blur, matsize, dataInt, matsize);

	printf("Iterating: \n");
	//Iterating
	for (int m = 0; m < config.numit; m++) {

		if (config.beta > 0) {
			softargmax(blur, size[0], size[1], config.beta);
		}

		if (config.psig > 0) {
			point_smoothing(blur, barr, size[0], size[1], abs((int)config.psig));
		}


		if (config.zsig!=0) {
			blur_it_CD(newblur, blur, zupind, zloind, lines, config.zsig);
			memcpy_s(blur, matsize, newblur, matsize);
			
		}
		
		
		fftconvolve2D(newblur, blur, peakshape, size);
		int i;
		//#pragma omp parallel for private (i), schedule(dynamic)
		for (i = 0; i < lines; i++) {
			if(newblur[i]!=0){
				newblur2[i] = dataInt[i] / newblur[i];
			}
			else {
				newblur2[i] = 0;
			}
		}
		fftconvolve2D(newblur, newblur2, peakshape, size);
		//#pragma omp parallel for private (i), schedule(dynamic)
		for (i = 0; i < lines; i++) {
			blur[i] = blur[i] * newblur[i];
		}
		if (config.numit < 10) { printf("Iteration: %d\n", m); }
		else { if (m % (config.numit / 10) == 0) { printf("Iteration: %d\n", m); } }
		//printf("Iteration: %d\n",m);
	}
	
	//Writing outputs

	//Reconvolve if necessary
	if (config.rawflag == 0) {
		memcpy_s(newblur, matsize, blur, matsize);
		fftconvolve2D(blur, newblur, peakshape, size);
	}
	
	//Normalize if necessary
	if (config.datanorm == 1) {
		float blurmax = getmax(lines, blur);
		if (dmax != 0) {
			normalize(lines, blur, blurmax / dmax);
		}
		printf("Maxes: %f %f\n",blurmax, dmax);
	}
	//Write the Output file
	ApplyCutoff1D(blur, 0, lines);
	
	// Write to binary
	//char* suffixfit = "";
	//write1D(config.outfile, suffixfit, blur, lines);

	// Write to text
	FILE* out_ptr = NULL;
	char outstring4[500];
	sprintf(outstring4, "%s_decon.txt", config.outfile);
	errno_t err = fopen_s(&out_ptr, outstring4, "w");
	if (err != 0) { printf("Error Opening %s %d\n", outstring4, err); exit(err); }
	for (int i = 0; i < lines; i++)
	{
		fprintf(out_ptr, "%f\n", blur[i]);
	}
	fclose(out_ptr);
	printf("Wrote Deconvolution Output to: %s\n", outstring4);


	free(mzdat);
	free(zdat);
	free(dataInt);

	free(blur);
	free(newblur);
	free(newblur2);
	free(peakshape);
	free(zloind);
	free(zupind);
	free(mzext);
	free(zext);
	
	endtime = time(NULL);
	printf("\nDone in %ds!\n", (int)difftime(endtime, starttime));
	return 0;
}
