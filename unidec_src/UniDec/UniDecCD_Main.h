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
	int i;
	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i < length; i++) {
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


void setup_blur_z(int * zupind, int * zloind, const float *mzdat, const float *zdat, const int lines, const float adductmass, const float mzranges[4], const int size[3])
{
	int i;
	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i < lines; i++) {
		float mz = mzdat[i];
		float z = zdat[i];
		int zind = (int)z - mzranges[2];
		float mass = calcMass(mz, z, adductmass);

		float uppermz = calcMz(mass, z + 1, adductmass);

		float lowermz = 0;
		if ((z - 1) != 0) {
			lowermz = calcMz(mass, z - 1, adductmass);
		}

		int mzupind = (int)roundf((uppermz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (size[0] - 1));
		if (z == mzranges[3] || mzupind < 0 || mzupind > size[0]) {
			zupind[i] = i;
		}
		else
		{
			int newind = index2D(size[1], mzupind, zind + 1);
			if (newind < 0) { newind = 0; }
			if (newind >= lines) { newind = lines - 1; }
			zupind[i] = newind;
		}

		int mzloind = (int)roundf((lowermz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (size[0] - 1));
		if (zind == 0 || mzloind < 0 || mzloind > size[0]) {
			zloind[i] = i;
		}
		else
		{
			int newind = index2D(size[1], mzloind, zind - 1);
			if (newind < 0) { newind = 0; }
			if (newind >= lines) { newind = lines - 1; }
			zloind[i] = newind;
		}
	}
}


void setup_blur_m(int* mupind, int* mloind, const float* mzdat, const float* zdat, const int lines, const float adductmass, const float mzranges[4], const int size[3], const float molig)
{
	int i;
	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i < lines; i++) {
		float mz = mzdat[i];
		float z = zdat[i];
		int zind = (int)z - mzranges[2];
		float mass = calcMass(mz, z, adductmass);

		float uppermz = calcMz(mass+molig, z, adductmass);
		float lowermz = calcMz(mass-molig, z, adductmass);

		int mzupind = (int)roundf((uppermz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (size[0] - 1));
		if (mzupind < 0 || mzupind > size[0]) {
			mupind[i] = i;
		}
		else
		{
			int newind = index2D(size[1], mzupind, zind);
			if (newind < 0) { newind = 0; }
			if (newind >= lines) { newind = lines - 1; }
			mupind[i] = newind;
		}

		int mzloind = (int)roundf((lowermz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (size[0] - 1));
		if (mzloind < 0 || mzloind > size[0]) {
			mloind[i] = i;
		}
		else
		{
			int newind = index2D(size[1], mzloind, zind);
			if (newind < 0) { newind = 0; }
			if (newind >= lines) { newind = lines - 1; }
			mloind[i] = newind;
		}
	}
}

void precompute_fft2D(float* list2, int* size, fftw_complex * out2) {
	fftw_complex* in2;
	fftw_plan p2;

	int i;
	int length = size[0] * size[1];
	in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i < length; i++)
	{
		in2[i][0] = list2[i];
		in2[i][1] = 0;
	}
	p2 = fftw_plan_dft_2d(size[0], size[1], in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p2);
	fftw_destroy_plan(p2);
	fftw_free(in2);
}

void fftconvolve2D_precomputed(float* corr, float* list1, fftw_complex * out2, int* size, fftw_plan p1, fftw_plan p3, fftw_complex * in1, fftw_complex * out1)
{
	int i;
	int length = size[0] * size[1];
	
	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i < length; i++)
	{
		in1[i][0] = list1[i];
		in1[i][1] = 0;
	}
	
	fftw_execute(p1);

	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i < length; i++)
	{
		float a = out1[i][0];
		float b = out1[i][1];
		float c = out2[i][0];
		float d = out2[i][1];
		in1[i][0] = a * c - b * d;
		in1[i][1] = b * c + a * d;
	}

	fftw_execute(p3);

	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i < length; i++)
	{
		corr[i] = out1[i][0] / length;
	}
}

void complex_conjugate(fftw_complex* in, fftw_complex* out, int length)
{
	int i;
	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i < length; i++)
	{
		out[i][0] = in[i][0];
		out[i][1] = -in[i][1];
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
	int* mupind = NULL, * mloind = NULL;
	char* barr = NULL;

	float* mzdat = NULL, * mzext = NULL,
		* zdat = NULL, * zext = NULL,
		* dataInt = NULL,
		* peakshape = NULL,
		*mkernel = NULL,
		*blur = NULL, *newblur=NULL, *newblur2=NULL, *oldblur=NULL
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

	// Set Up FFT stuff
	fftw_complex* in1, * out1;
	fftw_plan p1, p3;
	in1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lines);
	out1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lines);
	p1 = fftw_plan_dft_2d(size[0], size[1], in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
	p3 = fftw_plan_dft_2d(size[0], size[1], in1, out1, FFTW_BACKWARD, FFTW_ESTIMATE);

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

	// Set Up Peak Shapes
	peakshape = calloc(lines, sizeof(float));
	mkernel = calloc(lines, sizeof(float));
	fftw_complex* peakshape_FFT, *inverse_peakshape_FFT, *mkernel_FFT;
	peakshape_FFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lines);
	inverse_peakshape_FFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lines);
	mkernel_FFT = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * lines);

	//Make Peak Shape Kernel
	GetPeaks(peakshape, size, mzext, zext, config.mzsig, config.csig, config.psfun, config.zpsfun);
	// Precompute FFTs
	precompute_fft2D(peakshape, size, peakshape_FFT);
	complex_conjugate(peakshape_FFT, inverse_peakshape_FFT, lines);

	printf("Peak Shape Set\n");


	//Set up charge blur
	zupind = calloc(lines, sizeof(int));
	zloind = calloc(lines, sizeof(int));
	if (config.zsig != 0) { setup_blur_z(zupind, zloind, mzdat, zdat, lines, config.adductmass, mzranges, size); 
	printf("Charge Blur Set\n");
	}

	//Set up mass blur
	mupind = calloc(lines, sizeof(int));
	mloind = calloc(lines, sizeof(int));
	if (config.msig != 0) { setup_blur_m(mupind, mloind, mzdat, zdat, lines, config.adductmass, mzranges, size, config.molig); 
	printf("Mass Blur Set: %f %f\n", config.molig, config.msig);	
	}


	//Set up iteration
	blur = calloc(lines, sizeof(float));
	newblur = calloc(lines, sizeof(float));
	newblur2 = calloc(lines, sizeof(float));
	oldblur = calloc(lines, sizeof(float));
	barr = calloc(lines, sizeof(char));
	for (int i = 0; i < lines; i++) { barr[i] = 1; }
	rsize_t matsize = lines * sizeof(float);
	memcpy(blur, dataInt, matsize);
	memcpy(oldblur, blur, matsize);

	printf("Iterating: \n");
	//Iterating
	float conv=0;
	int off = 0;
	for (int m = 0; m < config.numit; m++) {
		// Apply softmax
		if (config.beta > 0) {
			softargmax(blur, size[0], size[1], config.beta);
		}
		// Apply point smoothing
		if (config.psig > 0) {
			point_smoothing(blur, barr, size[0], size[1], abs((int)config.psig));
		}
		// Apply charge smoothing
		if (config.zsig!=0) {
			blur_it_CD(newblur, blur, zupind, zloind, lines, config.zsig);
			memcpy(blur, newblur, matsize);
		}
		// Apply mass smoothing
		if (config.msig != 0) {
			blur_it_CD(newblur, blur, mupind, mloind, lines, config.msig);
			memcpy(blur, newblur, matsize);
		}
		
		// Richardson Lucy
		fftconvolve2D_precomputed(newblur, blur, peakshape_FFT, size, p1, p3, in1, out1);
		int i;
		#pragma omp parallel for private (i), schedule(dynamic)
		for (i = 0; i < lines; i++) {
			if(newblur[i]!=0){
				newblur2[i] = dataInt[i] / newblur[i];
			}
			else {
				newblur2[i] = 0;
			}
		}
		fftconvolve2D_precomputed(newblur, newblur2, inverse_peakshape_FFT, size, p1, p3, in1, out1);
		#pragma omp parallel for private (i), schedule(dynamic)
		for (i = 0; i < lines; i++) {
			blur[i] = blur[i] * newblur[i];
		}
		
		//Determine the metrics for conversion. Only do this every 10% to speed up. Stop loop if converged
		if ((config.numit < 10 || m % 10 == 0 || m % 10 == 1 || m>0.9 * config.numit)) {
			float diff = 0;
			float tot = 0;
			for (int i = 0; i < lines; i++)
			{
				diff += pow((blur[i] - oldblur[i]), 2);
				tot += blur[i];
			}
			if (tot != 0) { conv = (diff / tot); }
			else {
				if (conv == 12345678) { printf("m/z vs. charge grid is zero. Iteration: %d\n", m); break; }
				else { conv = 12345678; }
			}

			if (conv < 0.000001) {
				if (off == 1 && config.numit > 0) {
					printf("Converged in %d iterations.\n", m);
					break;
				}
				off = 1;
			}
			memcpy(oldblur, blur, matsize);
		}
		
	}
	printf("Completed Iterations\n");
	//Writing outputs

	//Outputting Fit Reconvolved Data as newblur2
	memcpy(newblur, blur, matsize);
	fftconvolve2D_precomputed(newblur2, newblur, peakshape_FFT, size, p1, p3, in1, out1);
	//Normalize if necessary
	if (config.datanorm == 1) {
		float blurmax = getmax(lines, newblur2);
		if (dmax != 0) {
			normalize(lines, newblur2, blurmax / dmax);
		}
	}
	ApplyCutoff1D(newblur2, 0, lines);
	//Write the fitdat file
	char* suffixfit = "fitdat";
	write1D(config.outfile, suffixfit, newblur2, lines);


	//Writing Main Output

	//Reconvolve if necessary
	if (config.rawflag == 0) {
		GetPeaks(mkernel, size, mzext, zext, config.mzsig, 0, config.psfun, config.zpsfun);
		precompute_fft2D(mkernel, size, mkernel_FFT);
		
		memcpy(newblur, blur, matsize);
		fftconvolve2D_precomputed(blur, newblur, mkernel_FFT, size, p1, p3, in1, out1);
		printf("Reconvolved with m/z dimension\n");
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
	out_ptr = fopen(outstring4, "w");
	if (out_ptr == 0) { printf("Error Opening %s\n", outstring4); exit(1); }
	for (int i = 0; i < lines; i++)
	{
		fprintf(out_ptr, "%f\n", blur[i]);
	}
	fclose(out_ptr);
	printf("Wrote Deconvolution Output to: %s\n", outstring4);

	//Free memory
	free(mzdat);
	free(zdat);
	free(dataInt);

	free(blur);
	free(newblur);
	free(newblur2);
	free(oldblur);

	free(peakshape);
	fftw_free(peakshape_FFT);
	free(mkernel);
	fftw_free(mkernel_FFT);
	fftw_free(inverse_peakshape_FFT);

	free(zloind);
	free(zupind);
	free(mloind);
	free(mupind);
	free(mzext);
	free(zext);

	fftw_destroy_plan(p1);
	fftw_destroy_plan(p3);
	fftw_free(in1); fftw_free(out1);
	
	endtime = time(NULL);
	printf("\nDone in %ds!\n", (int)difftime(endtime, starttime));
	return 0;
}
