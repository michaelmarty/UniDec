/*
* UniDec_Main.h
*
*  Created on : 29 December 2016
* Author : Michael.Marty
*/

//
// Copyright 2015 University of Oxford
// Copyright 2016 University of Arizona
//
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "UniDec.h"
#include "UD_H5_IO.h"
#include "UD_score.h"


Decon MainDeconvolution(const Config config, const Input inp, const int silent)
{
	Decon decon = SetupDecon();
	char* barr = NULL;

	int mlength, zlength, numclose,
		* mind = NULL,
		* zind = NULL,
		* closemind = NULL,
		* closezind = NULL,
		* closeind = NULL,
		* starttab = NULL,
		* endtab = NULL;
	float
		* mdist = NULL,
		* zdist = NULL,
		* mzdist = NULL,
		* rmzdist = NULL,
		* oldblur = NULL,
		* closeval = NULL,
		* closearray = NULL;

	//...................................................................
	//
	//     Sets the mzdist with the peak shape
	//
	//....................................................................
	barr = calloc(config.lengthmz * config.numz, sizeof(char));
	memcpy(barr, inp.barr, config.lengthmz * config.numz * sizeof(char));

	//Sets a threshold for m/z values to check. Things that are far away in m/z space don't need to be considered in the iterations.
	float threshold = config.psthresh * fabs(config.mzsig) * config.peakshapeinflate;
	if (silent == 0) { printf("Threshold: %f\t", threshold); }
	//Create a list of start and end values to box in arrays based on the above threshold
	starttab = calloc(config.lengthmz, sizeof(int));
	endtab = calloc(config.lengthmz, sizeof(int));
	int maxlength = 1;
	if (config.mzsig != 0) {
		//Gets maxlength and sets start and endtab
		maxlength = SetStartsEnds(config, &inp, starttab, endtab, threshold);

		//Changes dimensions of the peak shape function. 1D for speedy and 2D otherwise
		int pslen = config.lengthmz;
		if (config.speedyflag == 0) { pslen = config.lengthmz * maxlength; }
		mzdist = calloc(pslen, sizeof(float));
		rmzdist = calloc(pslen, sizeof(float));
		memset(mzdist, 0, pslen * sizeof(float));
		memset(rmzdist, 0, pslen * sizeof(float));

		int makereverse = 0;
		if (config.mzsig < 0 || config.beta < 0) { makereverse = 1; }

		//Calculates the distance between mz values as a 2D or 3D matrix
		if (config.speedyflag == 0)
		{
			MakePeakShape2D(config.lengthmz, maxlength, starttab, endtab, inp.dataMZ, fabs(config.mzsig) * config.peakshapeinflate, config.psfun, config.speedyflag, mzdist, rmzdist, makereverse);
		}
		else
		{
			//Calculates peak shape as a 1D list centered at the first element for circular convolutions
			MakePeakShape1D(inp.dataMZ, threshold, config.lengthmz, config.speedyflag, fabs(config.mzsig) * config.peakshapeinflate, config.psfun, mzdist, rmzdist, makereverse);
		}
		if (silent == 0) { printf("mzdist set: %f\t maxlength: %d\n", mzdist[0], maxlength); }

	}
	else
	{
		mzdist = calloc(0, sizeof(float));
		maxlength = 0;
	}

	//....................................................
	//
	//    Setting up the neighborhood blur
	//
	//......................................................

	//sets some parameters regarding the neighborhood blur function
	if (config.zsig >= 0 && config.msig >= 0) {
		zlength = 1 + 2 * (int)config.zsig;
		mlength = 1 + 2 * (int)config.msig;
	}
	else {
		if (config.zsig != 0) { zlength = 1 + 2 * (int)(3 * fabs(config.zsig) + 0.5); }
		else { zlength = 1; }
		if (config.msig != 0) { mlength = 1 + 2 * (int)(3 * fabs(config.msig) + 0.5); }
		else { mlength = 1; }
	}
	numclose = mlength * zlength;

	//Sets up the blur function in oligomer mass and charge
	mind = calloc(mlength, sizeof(int));
	mdist = calloc(mlength, sizeof(float));

	for (int i = 0; i < mlength; i++)
	{
		mind[i] = i - (mlength - 1) / 2;
		if (config.msig != 0) { mdist[i] = exp(-(pow((i - (mlength - 1) / 2.), 2)) / (2.0 * config.msig * config.msig)); }
		else { mdist[i] = 1; }
	}

	zind = calloc(zlength, sizeof(int));
	zdist = calloc(zlength, sizeof(float));
	for (int i = 0; i < zlength; i++)
	{
		zind[i] = i - (zlength - 1) / 2;
		if (config.zsig != 0) { zdist[i] = exp(-(pow((i - (zlength - 1) / 2.), 2)) / (2.0 * config.zsig * config.zsig)); }
		else { zdist[i] = 1; }
		//printf("%f\n", zdist[i]);
	}

	//Initializing memory
	closemind = calloc(numclose, sizeof(int));
	closezind = calloc(numclose, sizeof(int));
	closeval = calloc(numclose, sizeof(float));
	closeind = calloc(numclose * config.lengthmz * config.numz, sizeof(int));
	closearray = calloc(numclose * config.lengthmz * config.numz, sizeof(float));

	//Determines the indexes of things that are close as well as the values used in the neighborhood convolution
	for (int k = 0; k < numclose; k++)
	{
		closemind[k] = mind[k % mlength];
		closezind[k] = zind[(int)k / mlength];
		closeval[k] = zdist[(int)k / mlength] * mdist[k % mlength];
	}
	simp_norm_sum(mlength, mdist);
	simp_norm_sum(zlength, zdist);
	simp_norm_sum(numclose, closeval);

	//Set up blur
	//MakeBlur(config.lengthmz, config.numz, numclose, barr, closezind, closemind, inp.mtab, config.molig, config.adductmass, inp.nztab, inp.dataMZ, closeind, threshold, config);
	MakeSparseBlur(numclose, barr, closezind, closemind, inp.mtab, inp.nztab, inp.dataMZ, closeind, closeval, closearray, config);

	int badness = 1;
	for (int i = 0; i < config.lengthmz * config.numz; i++)
	{
		if (barr[i] == 1) { badness = 0;}
	}
	if (badness == 1) { printf("ERROR: Setup is bad. No points are allowed\n"); exit(10); }

	if (silent == 0) { printf("Charges blurred: %d  Masses blurred: %d\n", zlength, mlength); }
	
	//IntPrint(closeind, numclose * config.lengthmz * config.numz);

	//Determine the maximum intensity in the data
	float dmax = Max(inp.dataInt, config.lengthmz);
	float betafactor = 1;
	if (dmax > 1) { betafactor=dmax; }

	//...................................................
	//
	//  Setting up and running the iteration
	//
	//.........................................................


	//Applies the intensity threshold to kill peaks
	if (config.intthresh != -1) { KillB(inp.dataInt, barr, config.intthresh, config.lengthmz, config.numz, config.isolength, inp.isotopepos, inp.isotopeval); }
	

	//Creates an intial probability matrix, decon.blur, of 1 for each element
	decon.blur = calloc(config.lengthmz * config.numz, sizeof(float));
	decon.newblur = calloc(config.lengthmz * config.numz, sizeof(float));
	oldblur = calloc(config.lengthmz * config.numz, sizeof(float));

	if (config.baselineflag == 1) {
		printf("Auto Baseline Mode On: %d\n", config.aggressiveflag);
		decon.baseline = calloc(config.lengthmz, sizeof(float));
		decon.noise = calloc(config.lengthmz, sizeof(float));
	}


	//#pragma omp parallel for private (i,j), schedule(auto)
	for (int i = 0; i < config.lengthmz; i++)
	{
		float val = inp.dataInt[i] / ((float)(config.numz + 2));
		if (config.baselineflag == 1) {

			decon.baseline[i] = val;
			decon.noise[i] = val;
		}

		for (int j = 0; j < config.numz; j++)
		{
			if (barr[index2D(config.numz, i, j)] == 1) {
				if (config.isotopemode == 0) {
					decon.blur[index2D(config.numz, i, j)] = val;
				}
				else { decon.blur[index2D(config.numz, i, j)] = 1; }
			}
			else
			{
				decon.blur[index2D(config.numz, i, j)] = 0;
			}
		}
	}
	memcpy(oldblur, decon.blur, sizeof(float) * config.lengthmz * config.numz);
	memcpy(decon.newblur, decon.blur, sizeof(float) * config.lengthmz * config.numz);



	float* dataInt2 = NULL;
	dataInt2 = calloc(config.lengthmz, sizeof(float));
	memcpy(dataInt2, inp.dataInt, sizeof(float) * config.lengthmz);
	if (config.baselineflag == 1)
	{
		if (config.mzsig != 0)
		{
			deconvolve_baseline(config.lengthmz, inp.dataMZ, inp.dataInt, decon.baseline, fabs(config.mzsig));
			if (config.aggressiveflag == 2)
			{
				for (int i = 0; i < 10; i++)
				{
					deconvolve_baseline(config.lengthmz, inp.dataMZ, inp.dataInt, decon.baseline, fabs(config.mzsig));
				}
				for (int i = 0; i < config.lengthmz; i++)
				{
					if (decon.baseline[i] > 0) {
						dataInt2[i] -= decon.baseline[i];
					}
				}
				//memcpy(decon.baseline, dataInt2, sizeof(float)*config.lengthmz);
			}
		}
		else
		{
			printf("Ignoring baseline subtraction because peak width is 0\n");
		}

	}


	//Run the iteration
	float blurmax = 0;
	decon.conv = 0;
	int off = 0;
	if (silent == 0) { printf("Iterating..."); }

	for (int iterations = 0; iterations < abs(config.numit); iterations++)
	{
		decon.iterations = iterations;
		if (config.beta > 0 && iterations > 0)
		{

			softargmax(decon.blur, config.lengthmz, config.numz, config.beta/betafactor);
			//printf("Beta %f\n", beta);
		}
		else if (config.beta < 0 && iterations >0)
		{
			softargmax_transposed(decon.blur, config.lengthmz, config.numz, fabs(config.beta/betafactor), barr, maxlength, config.isolength, inp.isotopepos, inp.isotopeval, config.speedyflag, starttab, endtab, rmzdist, config.mzsig);
		}

		if (config.psig >= 1 && iterations > 0)
		{
			point_smoothing(decon.blur, barr, config.lengthmz, config.numz, abs((int)config.psig));
			//printf("Point Smoothed %f\n", config.psig);
		}
		else if (config.psig < 0 && iterations >0)
		{
			point_smoothing_peak_width(config.lengthmz, config.numz, maxlength, starttab, endtab, mzdist, decon.blur, config.speedyflag, barr);
		}

		
		//Run Blurs
		if (config.zsig >= 0 && config.msig >= 0) {
			blur_it_mean(config.lengthmz, config.numz, numclose, closeind, decon.newblur, decon.blur, barr, closearray, config.zerolog);
		}
		else if (config.zsig > 0 && config.msig < 0)
		{
			blur_it_hybrid1(config.lengthmz, config.numz, zlength, mlength, closeind, closemind, closezind, mdist, zdist, decon.newblur, decon.blur, barr, closearray, config.zerolog);
		}
		else if (config.zsig < 0 && config.msig > 0)
		{
			blur_it_hybrid2(config.lengthmz, config.numz, zlength, mlength, closeind, closemind, closezind, mdist, zdist, decon.newblur, decon.blur, barr, closearray, config.zerolog);
		}
		else {
			blur_it(config.lengthmz, config.numz, numclose, closeind, closearray, decon.newblur, decon.blur, barr);
		}

		//Run Richardson-Lucy Deconvolution
		deconvolve_iteration_speedy(config.lengthmz, config.numz, maxlength,
			decon.newblur, decon.blur, barr, config.aggressiveflag, dataInt2,
			config.isolength, inp.isotopepos, inp.isotopeval, starttab, endtab, mzdist, rmzdist, config.speedyflag,
			config.baselineflag, decon.baseline, decon.noise, config.mzsig, inp.dataMZ, config.filterwidth, config.psig);
		
		//Determine the metrics for conversion. Only do this every 10% to speed up.
		if ((config.numit < 10 || iterations % 10 == 0 || iterations % 10 == 1 || iterations>0.9 * config.numit)) {
			float diff = 0;
			float tot = 0;
			for (int i = 0; i < config.lengthmz * config.numz; i++)
			{
				if (barr[i] == 1)
				{
					diff += pow((decon.blur[i] - oldblur[i]), 2);
					tot += decon.blur[i];
				}
			}
			if (tot != 0) { decon.conv = (diff / tot); }
			else { decon.conv = 12345678; printf("m/z vs. charge grid is zero. Iteration: %d\n", iterations); }

			//printf("Iteration: %d Convergence: %f\n", iterations, decon.conv);
			if (decon.conv < 0.000001) {
				if (off == 1 && config.numit > 0) {
					printf("Converged in %d iterations.\n\n", iterations);
					break;
				}
				off = 1;
			}
			memcpy(oldblur, decon.blur, config.lengthmz * config.numz * sizeof(float));
		}
		
	}

	free(dataInt2);
	free(oldblur);

	//................................................................
	//
	//     Setting up the outputs
	//
	//...............................................................

	

	//Reset the peak shape if it was inflated
	if (config.peakshapeinflate != 1 && config.mzsig != 0) {
		if (config.speedyflag == 0)
		{
			MakePeakShape2D(config.lengthmz, maxlength, starttab, endtab, inp.dataMZ, fabs(config.mzsig), config.psfun, config.speedyflag, mzdist, rmzdist, 0);
		}
		else
		{
			MakePeakShape1D(inp.dataMZ, threshold, config.lengthmz, config.speedyflag, fabs(config.mzsig), config.psfun, mzdist, rmzdist, 0);
		}
		printf("mzdist reset: %f\n", config.mzsig);
	}


	//Determine the maximum intensity in the blur matrix
	blurmax = Max(decon.blur, config.lengthmz * config.numz);
	float cutoff = 0;
	if (blurmax != 0) { cutoff = 0.000001; }

	//Apply The Cutoff
	ApplyCutoff1D(decon.blur, blurmax * cutoff, config.lengthmz * config.numz);
	
	
	//Calculate the fit data and error.
	decon.fitdat = calloc(config.lengthmz, sizeof(float));
	decon.error = errfunspeedy(config, decon, barr, inp.dataInt, maxlength, inp.isotopepos, inp.isotopeval, starttab, endtab, mzdist, &decon.rsquared);

	//Fix issues with fitdat and consecutive zero data points
	//TODO: It might be possible to build this in to convolve_simp so that this isn't necessary but it would require a 1D barr.
	if (config.intthresh != -1)
	{
		#pragma omp parallel for schedule(auto)
		for (int i = 0; i < config.lengthmz-1; i++)
		{
			if (inp.dataInt[i] == 0 && inp.dataInt[i + 1] == 0)
			{
				decon.fitdat[i] = 0;
				decon.fitdat[i + 1] = 0;
			}
		}
	}
	

	// Charge scaling (orbimode)
	if (config.orbimode == 1)
	{
		printf("Rescaling charge states and normalizing ");
		charge_scaling(decon.blur, inp.nztab, config.lengthmz, config.numz);
		//simp_norm(config.lengthmz * config.numz, decon.blur);
		printf("Done\n");
	}

	//Change Monoisotopic to Average if necessary
	if (config.isotopemode == 2)
	{
		monotopic_to_average(config.lengthmz, config.numz, decon.blur, barr, config.isolength, inp.isotopepos, inp.isotopeval);
	}

	//newblur is repurposed as the convolution of blur by the mz peaks shape
	float newblurmax = blurmax;
	if ((config.rawflag == 0 || config.rawflag == 2)) {
		if (config.mzsig != 0) {
			newblurmax = Reconvolve(config.lengthmz, config.numz, maxlength, starttab, endtab, mzdist, decon.blur, decon.newblur, config.speedyflag, barr);
		}
		else
		{
			memcpy(decon.newblur, decon.blur, config.lengthmz * config.numz * sizeof(float));
		}
	}


	//.......................................................
	//
	//  Mass space outputs
	//
	//..........................................................

	//Determine the maximum and minimum allowed masses.
	float massmax = config.masslb;
	float massmin = config.massub;
	if (config.fixedmassaxis == 0) {
		for (int i = 0; i < config.lengthmz; i++)
		{
			for (int j = 0; j < config.numz; j++)
			{
				if (decon.newblur[index2D(config.numz, i, j)] * barr[index2D(config.numz, i, j)] > newblurmax * cutoff)
				{
					float testmax = inp.mtab[index2D(config.numz, i, j)] + threshold * inp.nztab[j]+config.massbins;
					float testmin = inp.mtab[index2D(config.numz, i, j)] - threshold * inp.nztab[j];

					//To prevent really wierd decimals
					testmin = round(testmin / config.massbins) * config.massbins;
					testmax = round(testmax / config.massbins) * config.massbins;

					if (testmax > massmax){	massmax = testmax;}
					if (testmin < massmin){	massmin = testmin;}
				}
			}
		}
		printf("Massmin: %f  ", massmin);
		printf("Massmax: %f  ", massmax);
	}
	else { massmax = config.massub; massmin = config.masslb; }

	//Checks to make sure the mass axis is good and makes a dummy axis if not
	decon.mlen = (int)(massmax - massmin) / config.massbins;
	if (decon.mlen < 1) {
		printf("Bad mass axis length: %d\n", decon.mlen);
		massmax = config.massub;
		massmin = config.masslb;
		decon.mlen = (int)(massmax - massmin) / config.massbins;

		//Declare the memory
		decon.massaxis = calloc(decon.mlen, sizeof(float));
		decon.massaxisval = calloc(decon.mlen, sizeof(float));
		decon.massgrid = calloc(decon.mlen * config.numz, sizeof(float));
		memset(decon.massaxisval, 0, decon.mlen * sizeof(float));
		memset(decon.massgrid, 0, decon.mlen * config.numz * sizeof(float));

		//Create the mass axis
		for (int i = 0; i < decon.mlen; i++)
		{
			decon.massaxis[i] = massmin + i * config.massbins;
		}
		decon.uniscore = 0;
		printf("ERROR: No masses detected.\n");
	}
	else {

		//Declare the memory
		decon.massaxis = calloc(decon.mlen, sizeof(float));
		decon.massaxisval = calloc(decon.mlen, sizeof(float));
		decon.massgrid = calloc(decon.mlen * config.numz, sizeof(float));
		memset(decon.massaxisval, 0, decon.mlen * sizeof(float));
		memset(decon.massgrid, 0, decon.mlen * config.numz * sizeof(float));
		if (silent == 0) { printf("Mass axis length: %d\n", decon.mlen); }

		
		//Create the mass axis
		for (int i = 0; i < decon.mlen; i++)
		{
			decon.massaxis[i] = massmin + i * config.massbins;
		}

		//Determine the mass intensities from m/z grid
		if (config.poolflag == 0) {
			if (config.rawflag == 1 || config.rawflag == 3) {
				IntegrateTransform(config.lengthmz, config.numz, inp.mtab, massmax, massmin, decon.mlen, decon.massaxis, decon.massaxisval, decon.blur, decon.massgrid);
			}
			if (config.rawflag == 0 || config.rawflag == 2) {
				IntegrateTransform(config.lengthmz, config.numz, inp.mtab, massmax, massmin, decon.mlen, decon.massaxis, decon.massaxisval, decon.newblur, decon.massgrid);
			}
		}
		else if (config.poolflag == 1) {
			if (config.rawflag == 1 || config.rawflag == 3) {
				InterpolateTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
					inp.dataMZ, decon.massgrid, decon.massaxisval, decon.blur);
			}
			if (config.rawflag == 0 || config.rawflag == 2) {
				InterpolateTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
					inp.dataMZ, decon.massgrid, decon.massaxisval, decon.newblur);
			}
		}
		else if (config.poolflag == 2) {
			if (config.rawflag == 1 || config.rawflag == 3) {
				SmartTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
					inp.dataMZ, decon.massgrid, decon.massaxisval, decon.blur);
			}
			if (config.rawflag == 0 || config.rawflag == 2) {
				SmartTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
					inp.dataMZ, decon.massgrid, decon.massaxisval, decon.newblur);
			}
		}
		else {
			printf("Invalid poolflag %d\n", config.poolflag);
			exit(1987);
		}

	//.....................................
	// Scores
	// .......................................

	//Note this will not execute if the mass axis is bad
	float scorethreshold = 0;
	decon.uniscore = score(config, &decon, inp, scorethreshold);

	}

	//Free Memory
	free(mzdist);
	free(rmzdist);
	free(closeval);
	free(closearray);
	free(closemind);
	free(closezind);
	free(endtab);
	free(starttab);

	free(mdist);
	free(mind);
	free(zind);
	free(zdist);
	free(barr);
	free(closeind);
	return decon;
}

void RunAutotune(Config *config, Input *inp, Decon *decon) {
	printf("Starting Autotune...\n");
	float start_peakwindow = config->peakwin;
	float start_peakthresh = config->peakthresh;
	config->peakwin = 3 * config->massbins;
	config->peakthresh = 0.01;

	int start_numit = config->numit;
	config->numit = 10;
	*decon = MainDeconvolution(*config, *inp, 1);
	float bestscore = decon->uniscore;

	float start_mzsig = config->mzsig;
	float start_zsig = config->zsig;
	float start_beta = config->beta;
	float start_psig = config->psig;

	//float mz_sigs[5] = { 0, 0.1, 1, 10, 100};
	float mz_sigs[2] = { 0, config->mzsig };
	//float z_sigs[2] = { -1, 1 };
	float z_sigs[1] = { config->zsig };
	float betas[3] = { 0, 50, 500 };
	float psigs[3] = { 0, 1, 10 };

	//int n_mzsigs = 5;
	int n_mzsigs = 2;
	//int n_zsigs = 2;
	int n_zsigs = 1;
	int n_betas = 3;
	int n_psigs = 3;

	for (int i = 0; i < n_mzsigs; i++)
	{
		for (int j = 0; j < n_zsigs; j++)
		{
			for (int k = 0; k < n_betas; k++)
			{
				for (int h = 0; h < n_psigs; h++)
				{
					config->mzsig = mz_sigs[i];
					config->zsig = z_sigs[j];
					config->beta = betas[k];
					config->psig = psigs[h];

					*decon = MainDeconvolution(*config, *inp, 1);
					float newscore = decon->uniscore;
					if (newscore > bestscore) {
						start_mzsig = config->mzsig;
						start_zsig = config->zsig;
						start_beta = config->beta;
						start_psig = config->psig;
						bestscore = newscore;
					}
				}
			}
		}
	}

	config->numit = start_numit;
	config->mzsig = start_mzsig;
	config->zsig = start_zsig;
	config->beta = start_beta;
	config->psig = start_psig;
	config->peakthresh = start_peakthresh;
	config->peakwin = start_peakwindow;
	printf("Best mzsig: %f zsig: %f beta: %f psig: %f Score:%f\n", start_mzsig, start_zsig, start_beta, start_psig, bestscore);
	*decon = MainDeconvolution(*config, *inp, 0);
}

int run_unidec(int argc, char *argv[], Config config) {
	//Initialize
	clock_t starttime;
	starttime = clock();

	Input inp = SetupInputs();

	bool autotune = 0;
	if (argc>2)
	{
		if (strcmp(argv[2], "-test") == 0)
		{
			float testmass = atof(argv[3]);
			printf("Testing Mass: %f\n", testmass);
			test_isotopes(testmass, inp.isoparams);
			exit(0);
		}

		if (strcmp(argv[2], "-autotune") == 0){autotune = 1;}
	}

	//..................................
	//
	// File Inputs
	//
	//.................................

	if (argc >= 2)
	{
		ReadInputs(argc, argv, &config, &inp);
	}
	else{ exit(88); }


	//.............................................................
	//
	//          Set the Mass Values and Values to Exclude
	//
	//...............................................................

	//Fills inp.mtab by multiplying each mz by each z
	inp.mtab = calloc(config.lengthmz * config.numz, sizeof(float));
	#pragma omp parallel for schedule(auto)
	for (int i = 0; i<config.lengthmz; i++)
	{
		for (int j = 0; j<config.numz; j++)
		{
			inp.mtab[index2D(config.numz, i, j)] = (inp.dataMZ[i] * inp.nztab[j] - config.adductmass*inp.nztab[j]);
		}
	}

	//Allocates the memory for the boolean array of whether to test a value
	inp.barr = calloc(config.lengthmz * config.numz, sizeof(char));
	//Tells the algorithm to ignore data that are equal to zero
	ignorezeros(inp.barr, inp.dataInt, config.lengthmz, config.numz);

	//Sets limits based on mass range and any test masses
	SetLimits(config, &inp);

	//Manual Assignments
	if (config.manualflag == 1)
	{
		ManualAssign(inp.dataMZ, inp.barr, inp.nztab, config);
	}

	//Setup Isotope Distributions
	if (config.isotopemode > 0)
	{
		setup_and_make_isotopes(&config, &inp);
	}

	//................................................................
	//
	// Deconvolution
	//
	//...................................................................

	//Setup the Deconvolution
	Decon decon=SetupDecon();

	//Autotuning
	if (autotune == 1) {
		RunAutotune(&config, &inp, &decon);
	}
	else{ 
		//Run the main Deconvolution		
		decon = MainDeconvolution(config, inp, 0); }

	//................................................................
	//
	//  Wrapping up
	//
	//...................................................................

	//Write Everything
	WriteDecon(config, &decon, &inp);

	//Writes a file with a number of key parameters such as error and number of significant parameters.
	clock_t end = clock();
	float totaltime = (float)(end - starttime) / CLOCKS_PER_SEC;
	if (config.filetype == 0) {
		FILE* out_ptr = NULL;
		char outstring3[500];
		sprintf(outstring3, "%s_error.txt", config.outfile);
		out_ptr = fopen(outstring3, "w");
		fprintf(out_ptr, "error = %f\n", decon.error);
		fprintf(out_ptr, "time = %f\n", totaltime);
		fprintf(out_ptr, "iterations = %d\n", decon.iterations);
		fprintf(out_ptr, "uniscore = %f\n", decon.uniscore);
		fprintf(out_ptr, "mzsig = %f\n", config.mzsig);
		fprintf(out_ptr, "zzsig = %f\n", config.zsig);
		fprintf(out_ptr, "beta = %f\n", config.beta);
		fprintf(out_ptr, "psig = %f\n", config.psig);
		fclose(out_ptr);
		//printf("Stats and Error written to: %s\n", outstring3);
	}
	else {
		write_attr_float(config.file_id, config.dataset, "error", decon.error);
		write_attr_int(config.file_id, config.dataset, "iterations", decon.iterations);
		write_attr_float(config.file_id, config.dataset, "time", totaltime);
		write_attr_int(config.file_id, config.dataset, "length_mz", config.lengthmz);
		write_attr_int(config.file_id, config.dataset, "length_mass", decon.mlen);
		write_attr_float(config.file_id, config.dataset, "uniscore", decon.uniscore);
		write_attr_float(config.file_id, config.dataset, "rsquared", decon.rsquared);
		write_attr_float(config.file_id, config.dataset, "mzsig", config.mzsig);
		write_attr_float(config.file_id, config.dataset, "zsig", config.zsig);
		write_attr_float(config.file_id, config.dataset, "psig", config.psig);
		write_attr_float(config.file_id, config.dataset, "beta", config.beta);
		set_needs_grids(config.file_id);
		H5Fclose(config.file_id);
	}

	//Free memory
	FreeInputs(inp);
	FreeDecon(decon);

	//printf("Error in the Fit: %f\n", decon.error);

	//Final Check that iterations worked and a reporter of the time consumed
	printf("Finished with %d iterations in %f seconds!\n\n", decon.iterations, totaltime);
	return 0;
}