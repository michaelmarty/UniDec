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
	double
		* mdist = NULL,
		* zdist = NULL,
		* mzdist = NULL,
		* rmzdist = NULL,
		* oldblur = NULL,
		* closeval = NULL;

	//...................................................................
	//
	//     Sets the mzdist with the peak shape
	//
	//....................................................................
	barr = calloc(config.lengthmz * config.numz, sizeof(char));
	memcpy(barr, inp.barr, config.lengthmz * config.numz * sizeof(char));

	//Sets a threshold for m/z values to check. Things that are far away in m/z space don't need to be considered in the iterations.
	double threshold = config.psthresh * fabs(config.mzsig) * config.peakshapeinflate;
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
		mzdist = calloc(pslen, sizeof(double));
		rmzdist = calloc(pslen, sizeof(double));
		memset(mzdist, 0, pslen * sizeof(double));
		memset(rmzdist, 0, pslen * sizeof(double));

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
		if (silent == 0) { printf("mzdist set: %f\t", mzdist[0]); }
	}
	else
	{
		mzdist = calloc(0, sizeof(double));
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
	mdist = calloc(mlength, sizeof(double));

	for (int i = 0; i < mlength; i++)
	{
		mind[i] = i - (mlength - 1) / 2;
		if (config.msig != 0) { mdist[i] = exp(-(pow((i - (mlength - 1) / 2.), 2)) / (2.0 * config.msig * config.msig)); }
		else { mdist[i] = 1; }
	}

	zind = calloc(zlength, sizeof(int));
	zdist = calloc(zlength, sizeof(double));
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
	closeval = calloc(numclose, sizeof(double));
	closeind = calloc(numclose * config.lengthmz * config.numz, sizeof(double));

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
	MakeBlur(config.lengthmz, config.numz, numclose, barr, closezind, closemind, inp.mtab, config.molig, config.adductmass, inp.nztab, inp.dataMZ, closeind);

	if (silent == 0) { printf("Charges blurred: %d  Oligomers blurred: %d\n", zlength, mlength); }


	//...................................................
	//
	//  Setting up and running the iteration
	//
	//.........................................................

	//Creates an intial probability matrix, decon.blur, of 1 for each element
	decon.blur = calloc(config.lengthmz * config.numz, sizeof(double));
	decon.newblur = calloc(config.lengthmz * config.numz, sizeof(double));
	oldblur = calloc(config.lengthmz * config.numz, sizeof(double));

	if (config.baselineflag == 1) {
		printf("Auto Baseline Mode On: %d\n", config.aggressiveflag);
		decon.baseline = calloc(config.lengthmz, sizeof(double));
		decon.noise = calloc(config.lengthmz, sizeof(double));
	}


	//#pragma omp parallel for private (i,j), schedule(auto)
	for (int i = 0; i < config.lengthmz; i++)
	{
		double val = inp.dataInt[i] / ((float)(config.numz + 2));
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
	memcpy(oldblur, decon.blur, sizeof(double) * config.lengthmz * config.numz);
	memcpy(decon.newblur, decon.blur, sizeof(double) * config.lengthmz * config.numz);

	if (config.intthresh != 0) { KillB(inp.dataInt, barr, config.intthresh, config.lengthmz, config.numz, config.isolength, inp.isotopepos, inp.isotopeval); }

	double* dataInt2 = NULL;
	dataInt2 = calloc(config.lengthmz, sizeof(double));
	memcpy(dataInt2, inp.dataInt, sizeof(double) * config.lengthmz);
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
				//memcpy(decon.baseline, dataInt2, sizeof(double)*config.lengthmz);
			}
		}
		else
		{
			printf("Ignoring baseline subtraction because peak width is 0\n");
		}

	}


	//Run the iteration
	double blurmax = 0;
	decon.conv = 0;
	int off = 0;
	if (silent == 0) { printf("Iterating..."); }

	for (int iterations = 0; iterations < abs(config.numit); iterations++)
	{
		decon.iterations = iterations;
		if (config.beta > 0 && iterations > 0)
		{

			softargmax(decon.blur, config.lengthmz, config.numz, config.beta);
			//printf("Beta %f\n", beta);
		}
		else if (config.beta < 0 && iterations >0)
		{
			softargmax_transposed(decon.blur, config.lengthmz, config.numz, fabs(config.beta), barr, maxlength, config.isolength, inp.isotopepos, inp.isotopeval, config.speedyflag, starttab, endtab, rmzdist, config.mzsig);
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
			blur_it_mean(config.lengthmz, config.numz, numclose, closeind, decon.newblur, decon.blur, barr, config.zerolog);
		}
		else if (config.zsig > 0 && config.msig < 0)
		{
			blur_it_hybrid1(config.lengthmz, config.numz, zlength, mlength, closeind, closemind, closezind, mdist, zdist, decon.newblur, decon.blur, barr, config.zerolog);
		}
		else if (config.zsig < 0 && config.msig > 0)
		{
			blur_it_hybrid2(config.lengthmz, config.numz, zlength, mlength, closeind, closemind, closezind, mdist, zdist, decon.newblur, decon.blur, barr, config.zerolog);
		}
		else {
			blur_it(config.lengthmz, config.numz, numclose, closeind, closeval, decon.newblur, decon.blur, barr);
		}

		//Run Richardson-Lucy Deconvolution
		deconvolve_iteration_speedy(config.lengthmz, config.numz, maxlength,
			decon.newblur, decon.blur, barr, config.aggressiveflag, dataInt2,
			config.isolength, inp.isotopepos, inp.isotopeval, starttab, endtab, mzdist, rmzdist, config.speedyflag,
			config.baselineflag, decon.baseline, decon.noise, config.mzsig, inp.dataMZ, config.filterwidth, config.psig);

		//Determine the metrics for conversion. Only do this every 10% to speed up.
		if ((config.numit < 10 || iterations % 10 == 0 || iterations % 10 == 1 || iterations>0.9 * config.numit)) {
			double diff = 0;
			double tot = 0;
			for (int i = 0; i < config.lengthmz * config.numz; i++)
			{
				if (barr[i] == 1)
				{
					diff += pow((decon.blur[i] - oldblur[i]), 2);
					tot += decon.blur[i];
				}
			}
			if (tot != 0) { decon.conv = (diff / tot); }
			else { decon.conv = 12345678; }

			//printf("Iteration: %d Convergence: %f\n", iterations, decon.conv);
			if (decon.conv < 0.000001) {
				if (off == 1 && config.numit > 0) {
					printf("Converged in %d iterations.\n\n", iterations);
					break;
				}
				off = 1;
			}
			memcpy(oldblur, decon.blur, config.lengthmz * config.numz * sizeof(double));
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
	double cutoff = 0;
	if (blurmax != 0) { cutoff = 0.0000001 / blurmax; }

	//Apply The Cutoff
	ApplyCutoff1D(decon.blur, blurmax * cutoff, config.lengthmz * config.numz);

	//Calculate the fit data and error.
	decon.fitdat = calloc(config.lengthmz, sizeof(double));
	decon.error = errfunspeedy(inp.dataInt, decon.fitdat, decon.blur, config.lengthmz, config.numz, maxlength,
		config.isolength, inp.isotopepos, inp.isotopeval, starttab, endtab, mzdist, config.speedyflag);
	if (config.baselineflag == 1)
	{
#pragma omp parallel for schedule(auto)
		for (int i = 0; i < config.lengthmz; i++) {
			decon.fitdat[i] += decon.baseline[i];// +decon.noise[i];
			//decon.fitdat[i] = decon.noise[i]+0.1;
		}
	}
	ApplyCutoff1D(decon.fitdat, 0, config.lengthmz);

	// Charge scaling (orbimode)
	if (config.orbimode == 1)
	{
		printf("Rescaling charge states and normalizing ");
		charge_scaling(decon.blur, inp.nztab, config.lengthmz, config.numz);
		simp_norm(config.lengthmz * config.numz, decon.blur);
		printf("Done\n");
	}

	//Change Monoisotopic to Average if necessary
	if (config.isotopemode == 2)
	{
		monotopic_to_average(config.lengthmz, config.numz, decon.blur, barr, config.isolength, inp.isotopepos, inp.isotopeval);
	}

	//newblur is repurposed as the convolution of blur by the mz peaks shape
	double newblurmax = blurmax;
	if ((config.rawflag == 0 || config.rawflag == 2)) {
		if (config.mzsig != 0) {
			newblurmax = Reconvolve(config.lengthmz, config.numz, maxlength, starttab, endtab, mzdist, decon.blur, decon.newblur, config.speedyflag, barr);
		}
		else
		{
			memcpy(decon.newblur, decon.blur, config.lengthmz * config.numz * sizeof(double));
		}
	}


	//.......................................................
	//
	//  Mass space outputs
	//
	//..........................................................

	//Determine the maximum and minimum allowed masses.
	double massmax = config.masslb;
	double massmin = config.massub;
	if (config.fixedmassaxis == 0) {
		for (int i = 0; i < config.lengthmz; i++)
		{
			for (int j = 0; j < config.numz; j++)
			{
				if (decon.newblur[index2D(config.numz, i, j)] * barr[index2D(config.numz, i, j)] > newblurmax * cutoff)
				{
					if (inp.mtab[index2D(config.numz, i, j)] + threshold * inp.nztab[j] > massmax)
					{
						massmax = inp.mtab[index2D(config.numz, i, j)] + threshold * inp.nztab[j];
					}
					if (inp.mtab[index2D(config.numz, i, j)] - threshold * inp.nztab[j] < massmin)
					{
						massmin = inp.mtab[index2D(config.numz, i, j)] - threshold * inp.nztab[j];
					}
				}
			}
		}
		printf("Massmax: %f  ", massmax);
		printf("Massmin: %f  ", massmin);
	}
	else { massmax = config.massub; massmin = config.masslb; }

	//Performs an interpolation to get the zero charge mass values starting at config.masslb and going to config.massub in steps of config.massbins
	decon.mlen = (int)(massmax - massmin) / config.massbins;
	if (decon.mlen < 1) {
		printf("Bad mass axis length: %d\n", decon.mlen);
		massmax = config.massub;
		massmin = config.masslb;
		decon.mlen = (int)(massmax - massmin) / config.massbins;
	}
	else {

		//Declare the memory
		decon.massaxis = calloc(decon.mlen, sizeof(double));
		decon.massaxisval = calloc(decon.mlen, sizeof(double));
		decon.massgrid = calloc(decon.mlen * config.numz, sizeof(double));
		memset(decon.massaxisval, 0, decon.mlen * sizeof(double));
		memset(decon.massgrid, 0, decon.mlen * config.numz * sizeof(double));
		if (silent == 0) { printf("Mass axis length: %d\n", decon.mlen); }

		//To prevent really wierd decimals
		//Only do this if the axis isn't already fixed
		if (config.fixedmassaxis == 0) { massmin = round(massmin / config.massbins) * config.massbins; }

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
		else {
			if (config.rawflag == 1 || config.rawflag == 3) {
				InterpolateTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
					inp.dataMZ, inp.dataInt, decon.massgrid, decon.massaxisval, decon.blur);
			}
			if (config.rawflag == 0 || config.rawflag == 2) {
				InterpolateTransform(decon.mlen, config.numz, config.lengthmz, inp.nztab, decon.massaxis, config.adductmass,
					inp.dataMZ, inp.dataInt, decon.massgrid, decon.massaxisval, decon.newblur);
			}
		}

	}

	//.....................................
	// Scores
	// .......................................

	double scorethreshold = 0;
	decon.uniscore = score(config, decon.mlen, inp.dataMZ, inp.dataInt, decon.newblur, decon.massaxis, decon.massaxisval, decon.massgrid, inp.nztab, scorethreshold);

	//Free Memory
	free(mzdist);
	free(rmzdist);
	free(closeval);
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


int run_unidec(int argc, char *argv[], Config config) {
	//Initialize
	clock_t starttime;
	starttime = clock();

	FILE *out_ptr = NULL;
	hid_t file_id;

	Input inp = SetupInputs();

	char dataset[1024];
	char outdat[1024];
	bool autotune = 0;

	//Convert aggressiveflag to baselineflag
	if (config.aggressiveflag == 1 || config.aggressiveflag == 2) { config.baselineflag = 1; }
	else { config.baselineflag = 0; }

	//Experimental correction. Not sure why this is necessary.
	if (config.psig < 0) { config.mzsig /= 3; }

	//Set default parameters. These can be overwritten by config file.
	double isoparams[10] = { 1.00840852e+00, 1.25318718e-03, 2.37226341e+00, 8.19178000e-04,
		-4.37741951e-01, 6.64992972e-04, 9.94230511e-01, 4.64975237e-01, 1.00529041e-02, 5.81240792e-01 };

	if (argc>2)
	{
		if (strcmp(argv[2], "-test") == 0)
		{
			double testmass = atof(argv[3]);
			printf("Testing Mass: %f\n", testmass);
			test_isotopes(testmass, isoparams);
			exit(0);
		}

		if (strcmp(argv[2], "-autotune") == 0)
		{
			autotune = 1;
		}
	}

	if (config.metamode != -2)
	{
		strcpy(dataset, "/ms_dataset");
		char strval[1024];
		sprintf(strval, "/%d", config.metamode);
		strcat(dataset, strval);
		printf("HDF5 Data Set: %s\n", dataset);
	}
	else
	{
		strcpy(dataset, "/ms_data");
	}

	//..................................
	//
	// File Inputs
	//
	//.................................

	if (argc >= 2)
	{
		if (config.filetype == 1) {
			file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
			strjoin(dataset, "/processed_data", outdat);
			config.lengthmz = mh5getfilelength(file_id, outdat);
			inp.dataMZ = calloc(config.lengthmz, sizeof(double));
			inp.dataInt = calloc(config.lengthmz, sizeof(double));
			mh5readfile2d(file_id, outdat, config.lengthmz, inp.dataMZ, inp.dataInt);
			printf("Length of Data: %d \n", config.lengthmz);

			//Check the length of the mfile and then read it in.
			if (config.mflag == 1)
			{
				config.mfilelen = mh5getfilelength(file_id, "/config/masslist");
				printf("Length of mfile: %d \n", config.mfilelen);
				inp.testmasses = malloc(sizeof(double)*config.mfilelen);
				mh5readfile1d(file_id, "/config/masslist", inp.testmasses);
			}
			else {
				inp.testmasses = malloc(sizeof(double)*config.mfilelen);
			}
		}
		else {

			//Calculate the length of the data file automatically
			config.lengthmz = getfilelength(config.infile);
			inp.dataMZ = calloc(config.lengthmz, sizeof(double));
			inp.dataInt = calloc(config.lengthmz, sizeof(double));

			readfile(config.infile, config.lengthmz, inp.dataMZ, inp.dataInt);//load up the data array
			printf("Length of Data: %d \n", config.lengthmz);

			//Check the length of the mfile and then read it in.

			if (config.mflag == 1)
			{
				config.mfilelen = getfilelength(config.mfile);
				printf("Length of mfile: %d \n", config.mfilelen);
			}
			inp.testmasses = malloc(sizeof(double)*config.mfilelen);
			if (config.mflag == 1)
			{
				readmfile(config.mfile, config.mfilelen, inp.testmasses);//read in mass tab
			}
		}
	}

	//This for loop creates a list of charge values
	inp.nztab = calloc(config.numz, sizeof(int));
	for (int i = 0; i<config.numz; i++) { inp.nztab[i] = i + config.startz; }
	//printf("nzstart %d\n",inp.nztab[0]);

	//Test to make sure no charge state is zero
	for (int j = 0; j < config.numz; j++)
	{
		if (inp.nztab[j] == 0) { printf("Error: Charge state cannot be 0"); exit(100); }
	}
	//Test to make sure no two data points has the same x value
	for (int i = 0; i < config.lengthmz - 1; i++)
	{
		if (inp.dataMZ[i] == inp.dataMZ[i + 1]) { printf("Error: Two data points are identical"); exit(104); }
	}

	//.............................................................
	//
	//          Set the Mass Values and Values to Exclude
	//
	//...............................................................

	//Fills inp.mtab by multiplying each mz by each z
	inp.mtab = calloc(config.lengthmz * config.numz, sizeof(double));
	#pragma omp parallel for schedule(auto)
	for (int i = 0; i<config.lengthmz; i++)
	{
		for (int j = 0; j<config.numz; j++)
		{
			inp.mtab[index2D(config.numz, i, j)] = (inp.dataMZ[i] * inp.nztab[j] - config.adductmass*inp.nztab[j]);
		}
	}

	inp.barr = calloc(config.lengthmz * config.numz, sizeof(char));
	//Tells the algorithm to ignore data that are equal to zero
	ignorezeros(inp.barr, inp.dataInt, config.lengthmz, config.numz);

	//Sets limits based on mass range and any test masses
	SetLimits(config, &inp);

	//Manual Assignments
	if (config.manualflag == 1)
	{
		ManualAssign(config.manualfile, config.lengthmz, config.numz, inp.dataMZ, inp.barr, inp.nztab, file_id, config);
	}

	//Setup Isotope Distributions
	if (config.isotopemode > 0)
	{
		printf("Isotope Mode: %d\n", config.isotopemode);
		config.isolength = setup_isotopes(isoparams, inp.isotopepos, inp.isotopeval, inp.mtab, inp.nztab, inp.barr, inp.dataMZ, config.lengthmz, config.numz);

		inp.isotopepos = calloc(config.isolength * config.lengthmz * config.numz, sizeof(int));
		inp.isotopeval = calloc(config.isolength * config.lengthmz * config.numz, sizeof(float));

		make_isotopes(isoparams, inp.isotopepos, inp.isotopeval, inp.mtab, inp.nztab, inp.barr, inp.dataMZ, config.lengthmz, config.numz);

		printf("Isotopes set up, Length: %d\n", config.isolength);
	}

	Decon decon=SetupDecon();

	if (autotune == 1) {
		printf("Starting Autotune...\n");
		double start_peakwindow = config.peakwin;
		double start_peakthresh = config.peakthresh;
		config.peakwin = 3 * config.massbins;
		config.peakthresh = 0.01;

		int start_numit = config.numit;
		config.numit = 10;
		decon = MainDeconvolution(config, inp, 1);
		double bestscore = decon.uniscore;

		double start_mzsig = config.mzsig;
		double start_zsig = config.zsig;
		double start_beta = config.beta;
		double start_psig = config.psig;

		//double mz_sigs[5] = { 0, 0.1, 1, 10, 100};
		double mz_sigs[2] = { 0, config.mzsig };
		//double z_sigs[2] = { -1, 1 };
		double z_sigs[1] = { config.zsig };
		double betas[3] = { 0, 50, 500 };
		double psigs[3] = { 0, 1, 10 };

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
						config.mzsig = mz_sigs[i];
						config.zsig = z_sigs[j];
						config.beta = betas[k];
						config.psig = psigs[h];

						decon = MainDeconvolution(config, inp, 1);
						double newscore = decon.uniscore;
						if (newscore > bestscore) {
							start_mzsig = config.mzsig;
							start_zsig = config.zsig;
							start_beta = config.beta;
							start_psig = config.psig;
							bestscore = newscore;
						}
					}
				}
			}
		}

		config.numit = start_numit;
		config.mzsig = start_mzsig;
		config.zsig = start_zsig;
		config.beta = start_beta;
		config.psig = start_psig;
		config.peakthresh = start_peakthresh;
		config.peakwin = start_peakwindow;
		printf("Best mzsig: %f zsig: %f beta: %f psig: %f Score:%f\n", start_mzsig, start_zsig, start_beta, start_psig, bestscore);
		decon = MainDeconvolution(config, inp, 0);
	}
	else{ decon = MainDeconvolution(config, inp, 0); }

	//................................................................
	//
	//  Wrapping up
	//
	//...................................................................

	//Write Everything
	WriteDecon(config, &decon, &inp, file_id, dataset);

	//Writes a file with a number of key parameters such as error and number of significant parameters.
	clock_t end = clock();
	double totaltime = (double)(end - starttime) / CLOCKS_PER_SEC;
	if (config.filetype == 0) {
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
		printf("Stats and Error written to: %s\n", outstring3);
	}
	else {
		write_attr_double(file_id, dataset, "error", decon.error);
		write_attr_int(file_id, dataset, "iterations", decon.iterations);
		write_attr_double(file_id, dataset, "time", totaltime);
		write_attr_int(file_id, dataset, "length_mz", config.lengthmz);
		write_attr_int(file_id, dataset, "length_mass", decon.mlen);
		write_attr_double(file_id, dataset, "uniscore", decon.uniscore);
		write_attr_double(file_id, dataset, "mzsig", config.mzsig);
		write_attr_double(file_id, dataset, "zsig", config.zsig);
		write_attr_double(file_id, dataset, "psig", config.psig);
		write_attr_double(file_id, dataset, "beta", config.beta);
		set_needs_grids(file_id);
		H5Fclose(file_id);
	}

	//Free memory
	FreeInputs(inp);
	FreeDecon(decon);

	printf("\nError: %f\n", decon.error);

	//Final Check that iterations worked and a reporter of the time consumed
	printf("\nFinished with %d iterations and %f convergence in %f seconds!\n\n", decon.iterations, decon.conv, totaltime);
	return 0;
}