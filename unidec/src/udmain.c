//
// Created by mm96978 on 7/8/2025.
//
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

#include "udmain.h"


Decon ExitToBlank(const Config config, Decon decon)
{
	const int l = config.lengthmz * config.numz;
	decon.blur = calloc(l, sizeof(float));
	decon.newblur = calloc(l, sizeof(float));
	return decon;
}

Decon MainDeconvolution(const Config config, const Input inp, const int silent, const int verbose)
{
	//Insert time
	time_t starttime;
	starttime = clock();

	Decon decon = SetupDecon();
	char* barr = NULL;

	int mlength, zlength, numclose,
		* mind = NULL,
		* zind = NULL,
		* closemind = NULL,
		* closezind = NULL,
		* closeind = NULL;
	float * mdist = NULL,
		* zdist = NULL,
		* oldblur = NULL,
		* closeval = NULL,
		* closearray = NULL;

	//...................................................................
	//
	//     Sets the mzdist with the peak shape
	//
	//....................................................................
	int ln = config.lengthmz * config.numz;
	size_t sizelnc = (size_t)ln * sizeof(char);
	size_t sizelnf = (size_t)ln * sizeof(float);
	barr = calloc(ln, sizeof(char));
	memcpy(barr, inp.barr, sizelnc);

	int maxlength = SetUpPeakShape(config, inp, &decon, 0, 1);



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
		if (config.zsig != 0) { zlength = 1 + 2 * (int)(3 * fabsf(config.zsig) + 0.5); }
		else { zlength = 1; }
		if (config.msig != 0) { mlength = 1 + 2 * (int)(3 * fabsf(config.msig) + 0.5); }
		else { mlength = 1; }
	}
	numclose = mlength * zlength;

	//Sets up the blur function in oligomer mass and charge
	mind = calloc(mlength, sizeof(int));
	mdist = calloc(mlength, sizeof(float));

	for (int i = 0; i < mlength; i++)
	{
		mind[i] = i - (mlength - 1) / 2;
		if (config.msig != 0) { mdist[i] = expf(-(powf(((float) i - ((float) mlength - 1) / 2.f), 2.f)) / (2.0f * config.msig * config.msig)); }
		else { mdist[i] = 1; }
	}

	zind = calloc(zlength, sizeof(int));
	zdist = calloc(zlength, sizeof(float));
	for (int i = 0; i < zlength; i++)
	{
		zind[i] = i - (zlength - 1) / 2;
		if (config.zsig != 0) { zdist[i] = expf(-(powf(((float) i - ((float) zlength - 1) / 2.f), 2.f)) / (2.0f * config.zsig * config.zsig)); }
		else { zdist[i] = 1; }
		//printf("%f\n", zdist[i]);
	}

	//Initializing memory
	closemind = calloc(numclose, sizeof(int));
	closezind = calloc(numclose, sizeof(int));
	closeval = calloc(numclose, sizeof(float));
	int newlen = numclose * config.lengthmz * config.numz;
	closeind = calloc(newlen, sizeof(int));
	closearray = calloc(newlen, sizeof(float));

	//Determines the indexes of things that are close as well as the values used in the neighborhood convolution
	for (int k = 0; k < numclose; k++)
	{
		closemind[k] = mind[k % mlength];
		closezind[k] = zind[(int)k / mlength];
		closeval[k] = zdist[(int)k / mlength] * mdist[k % mlength];
	}
	NormSum(mdist, mlength);
	NormSum(zdist, zlength);
	NormSum(closeval, numclose);

	int badness1 = 1;
	for (int i = 0; i < config.lengthmz * config.numz; i++)
	{
		if (barr[i] == 1) { badness1 = 0; }
	}
	if (badness1 == 1) { printf("ERROR: No points are allowed. Setup is bad\n"); exit(10); }

	//Set up blur
	MakeSparseBlur(numclose, barr, closezind, closemind, closeind, closeval, closearray, config, &inp);

	if (silent == 0) { printf("Charges blurred: %d  Masses blurred: %d\n", zlength, mlength); }

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
	if (config.intthresh != -1) { KillB(inp.dataInt, barr, config.intthresh, config.lengthmz, config.numz); }


	//Creates an intial probability matrix, decon.blur, of 1 for each element
	decon.blur = calloc(ln, sizeof(float));
	decon.newblur = calloc(ln, sizeof(float));
	oldblur = calloc(ln, sizeof(float));

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
			deconvolve_baseline(config.lengthmz, inp.dataMZ, inp.dataInt, decon.baseline, fabsf(config.mzsig));
			if (config.aggressiveflag == 2)
			{
				for (int i = 0; i < 10; i++)
				{
					deconvolve_baseline(config.lengthmz, inp.dataMZ, inp.dataInt, decon.baseline, fabsf(config.mzsig));
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

	// Final check that there are are actually possible data in the set
	int badness = 1;
	for (int i = 0; i < config.lengthmz * config.numz; i++)
	{
		if (barr[i] == 1) { badness = 0; }
	}
	if (badness == 1) {
		if (silent == 0) { printf("ERROR: Setup is bad. No points are allowed.\n Check that either mass smoothing, charge smoothing, manual assignment, or isotope mode are on.\n"); };
		decon = ExitToBlank(config, decon);
		free(closeval);
		free(closearray);
		free(closemind);
		free(closezind);

		free(mdist);
		free(mind);
		free(zind);
		free(zdist);
		free(barr);
		free(closeind);
		return(decon);
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

		if (config.psig >= 1 && iterations > 0)
		{
			point_smoothing(decon.blur, barr, config.lengthmz, config.numz, abs((int)config.psig));
			//printf("Point Smoothed %f\n", config.psig);
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
			decon.starttab, decon.endtab, decon.mzdist, decon.rmzdist, config.speedyflag,
			config.baselineflag, decon.baseline, decon.noise, config.mzsig, inp.dataMZ, config.filterwidth, config.psig);

		//Determine the metrics for conversion. Only do this every 10% to speed up.
		if ((config.numit < 10 || iterations % 10 == 0 || iterations % 10 == 1 || iterations>0.9 * config.numit)) {
			float diff = 0;
			float tot = 0;
			for (int i = 0; i < config.lengthmz * config.numz; i++)
			{
				if (barr[i] == 1)
				{
					diff += powf((decon.blur[i] - oldblur[i]), 2);
					tot += decon.blur[i];
				}
			}
			if (tot != 0) { decon.conv = (diff / tot); }
			else {
				if (decon.conv == 12345678) { printf("m/z vs. charge grid is zero. Iteration: %d\n", iterations); break; }
				else{decon.conv = 12345678;}
			}

			//printf("Iteration: %d Convergence: %f\n", iterations, decon.conv);
			if (decon.conv < 0.000001) {
				if (off == 1 && config.numit > 0) {
					if (silent == 0) { printf("Converged in %d iterations.\n", iterations); }
					break;
				}
				off = 1;
			}
			memcpy(oldblur, decon.blur, sizelnf);
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
			MakePeakShape2D(config, &decon, maxlength, inp.dataMZ, 0, 0);
		}
		else
		{
			MakePeakShape1D(config, &decon, inp.dataMZ, 0, 0);
		}
		printf("mzdist reset: %f\n", config.mzsig);
	}


	//Determine the maximum intensity in the blur matrix
	blurmax = Max(decon.blur, config.lengthmz * config.numz);
	float cutoff = 0;
	if (blurmax != 0) { cutoff = 0.000001f; }

	//Apply The Cutoff
	ApplyCutoff(decon.blur, blurmax * cutoff, config.lengthmz * config.numz);


	//Calculate the fit data and error.
	decon.fitdat = calloc(config.lengthmz, sizeof(float));
	decon.error = errfunspeedy(config, decon, barr, inp.dataInt, maxlength, decon.starttab, decon.endtab, decon.mzdist, &decon.rsquared);

	//Fix issues with fitdat and consecutive zero data points
	//TODO: It might be possible to build this in to convolve_simp so that this isn't necessary but it would require a 1D barr.
	if (config.intthresh != -1)
	{
		int i;
		#pragma omp parallel for private(i)
		for (i = 0; i < config.lengthmz-1; i++)
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

	//newblur is repurposed as the convolution of blur by the mz peaks shape
	float newblurmax = blurmax;
	if ((config.rawflag == 0 || config.rawflag == 2)) {
		if (config.mzsig != 0) {
			newblurmax = Reconvolve(config, maxlength, &decon, barr);
		}
		else
		{
			memcpy(decon.newblur, decon.blur, sizelnf);
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
				if (decon.newblur[index2D(config.numz, i, j)] * (float) barr[index2D(config.numz, i, j)] > newblurmax * cutoff)
				{
					float testmax = inp.mtab[index2D(config.numz, i, j)] + config.psmzthresh * (float) inp.nztab[j]+config.massbins;
					float testmin = inp.mtab[index2D(config.numz, i, j)] - config.psmzthresh * (float) inp.nztab[j];

					//To prevent really wierd decimals
					testmin = roundf(testmin / config.massbins) * config.massbins;
					testmax = roundf(testmax / config.massbins) * config.massbins;

					if (testmax > massmax){	massmax = testmax;}
					if (testmin < massmin){	massmin = testmin;}
				}
			}
		}
		if (silent == 0){printf("Massmin: %f  ", massmin); printf("Massmax: %f  ", massmax);}
	}
	else { massmax = config.massub; massmin = config.masslb; }

	//Checks to make sure the mass axis is good and makes a dummy axis if not
	decon.mlen = (int)((massmax - massmin) / config.massbins);
	int mn = decon.mlen * config.numz;
	size_t sizemn = (size_t)mn * sizeof(float);
	if (decon.mlen < 1) {
		printf("ERROR: No masses detected. Length: %d\n", decon.mlen);
		massmax = config.massub;
		massmin = config.masslb;
		decon.mlen = (int)((massmax - massmin) / config.massbins);

		//Declare the memory
		decon.massaxis = calloc(decon.mlen, sizeof(float));
		decon.massaxisval = calloc(decon.mlen, sizeof(float));
		decon.massgrid = calloc(mn, sizeof(float));
		memset(decon.massaxisval, 0, decon.mlen * sizeof(float));
		memset(decon.massgrid, 0, sizemn);

		//Create the mass axis
		for (int i = 0; i < decon.mlen; i++)
		{
			decon.massaxis[i] = massmin + (float) i * config.massbins;
		}
		decon.uniscore = 0;
	}
	else {

		//Declare the memory
		decon.massaxis = calloc(decon.mlen, sizeof(float));
		decon.massaxisval = calloc(decon.mlen, sizeof(float));
		decon.massgrid = calloc(mn, sizeof(float));
		memset(decon.massaxisval, 0, decon.mlen * sizeof(float));
		memset(decon.massgrid, 0, sizemn);
		if (silent == 0) { printf("Mass axis length: %d\n", decon.mlen); }


		//Create the mass axis
		for (int i = 0; i < decon.mlen; i++)
		{
			decon.massaxis[i] = massmin + (float) i * config.massbins;
		}

		//Determine the mass intensities from m/z grid
		if (config.poolflag == 0) {
			IntegrateTransform(config, &decon, inp.mtab, massmax, massmin);
		}
		else if (config.poolflag == 1) {
			InterpolateTransform(config, &decon, &inp);
		}
		else if (config.poolflag == 2) {
			SmartTransform(config, &decon, &inp);
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
	decon.uniscore = score(config, &decon, inp, scorethreshold, silent);

	}
	// Print out time
	printf("Deconvolution Time: %f\n", (float)(clock() - starttime) / CLOCKS_PER_SEC);

	//Free Memory
	free(closeval);
	free(closearray);
	free(closemind);
	free(closezind);

	free(mdist);
	free(mind);
	free(zind);
	free(zdist);
	free(barr);
	free(closeind);
	return decon;
}

void RunAutotune(Config *config, const Input *inp, Decon *decon) {
	printf("Starting Autotune...\n");
	const float start_peakwindow = config->peakwin;
	const float start_peakthresh = config->peakthresh;
	config->peakwin = 3 * config->massbins;
	config->peakthresh = 0.01f;

	const int start_numit = config->numit;
	config->numit = 10;
	*decon = MainDeconvolution(*config, *inp, 1, 0);
	float bestscore = decon->uniscore;

	float start_mzsig = config->mzsig;
	float start_zsig = config->zsig;
	float start_beta = config->beta;
	float start_psig = config->psig;

	//float mz_sigs[5] = { 0, 0.1, 1, 10, 100};
	const float mz_sigs[2] = { 0, config->mzsig };
	//float z_sigs[2] = { -1, 1 };
	const float z_sigs[1] = { config->zsig };
	const float betas[3] = { 0, 50, 500 };
	const float psigs[3] = { 0, 1, 10 };

	//int n_mzsigs = 5;
	const int n_mzsigs = 2;
	//int n_zsigs = 2;
	const int n_zsigs = 1;
	const int n_betas = 3;
	const int n_psigs = 3;

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

					*decon = MainDeconvolution(*config, *inp, 1, 0);
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
	*decon = MainDeconvolution(*config, *inp, 0, 0);
}

int run_unidec_core(Config config, Input inp, Decon *decon, const int verbose, const int autotune) {
	//.............................................................
	//
	//          Set the Mass Values to Exclude
	//
	//...............................................................

	// Check inputs
	if (verbose ==1 ) {printf("Checking inputs: %f %d\n", inp.dataMZ[0], config.lengthmz);}

	CalcMasses(&config, &inp);

	//Allocates the memory for the boolean array of whether to test a value
	int newlen = config.lengthmz * config.numz;
	inp.barr = calloc(newlen, sizeof(char));
	//Tells the algorithm to ignore data that are equal to zero
	ignorezeros(inp.barr, inp.dataInt, config.lengthmz, config.numz);
	if (verbose == 1) { printf("Ignored Zeros. Length: %d\n", config.lengthmz);  }

	//Sets limits based on mass range and any test masses
	SetLimits(config, &inp);
	if (verbose == 1) { printf("Set limits\n"); }

	//Manual Assignments
	if (config.manualflag == 1)
	{
		ManualAssign(inp.dataMZ, inp.barr, inp.nztab, config);
		if (verbose == 1) { printf("Setup Manual Assign\n"); }
	}

	//................................................................
	//
	// Deconvolution
	//
	//...................................................................

	//Autotuning
	if (autotune == 1) {
		RunAutotune(&config, &inp, decon);
	}
	else{
		//Run the main Deconvolution
		*decon = MainDeconvolution(config, inp, config.silent, verbose); }

	// Run DoubleDec if it's checked
	if (config.doubledec) { // Use Python to confirm kernel file is selected
		printf("Running DoubleDec now\n");
		DoubleDecon(&config, decon);
	}
	if (verbose == 1) { printf("Decon Done\n"); }
	return 0;
}


int run_unidec(int argc, char *argv[], Config config) {
	//Initialize
	clock_t starttime;
	starttime = clock();

	bool autotune = 0;
	int verbose = 0;
	if (argc>2)
	{
		if (strcmp(argv[2], "-autotune") == 0){autotune = 1;}
		if (strcmp(argv[2], "-verbose") == 0) { verbose = 1; }
	}

	//..................................
	//
	// File Inputs
	//
	//.................................

	Input inp = SetupInputs();

	if (argc >= 2)
	{
		ReadInputs(argc, argv, &config, &inp);
	}
	else{ exit(88); }
	if (verbose == 1) { printf("Read Inputs\n"); }

	//Setup the Deconvolution
	Decon decon=SetupDecon();

	// ........................................................
	//
	//  Run the Deconvolution
	//
	//........................................................

	run_unidec_core(config, inp, &decon, verbose, autotune);

	//................................................................
	//
	//  Wrapping up
	//
	//...................................................................

	//Write Everything
	WriteDecon(config, &decon, &inp);
	if (verbose == 1) { printf("Write Done\n"); }

	//Writes a file with a number of key parameters such as error and number of significant parameters.
	clock_t end = clock();
	float totaltime = (float)(end - starttime) / CLOCKS_PER_SEC;
	if (config.filetype == 0) {
		FILE* out_ptr = NULL;
		char outstring3[1024];
		sprintf(outstring3, "%s_error.txt", config.outfile);
		out_ptr = fopen(outstring3, "w");
		if (out_ptr == 0) { printf("Error Opening %s\n", outstring3); exit(1); }
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
		//H5Fclose(config.file_id);
	}

	//Free memory
	FreeInputs(inp);
	FreeDecon(decon);

	//printf("Error in the Fit: %f\n", decon.error);

	//Final Check that iterations worked and a reporter of the time consumed
	if (config.silent == 0) { printf("Finished with %d iterations in %f seconds!\n\n", decon.iterations, totaltime); }
	return 0;
}