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

#ifdef _OPENMP
#include <omp.h>
#endif

typedef struct IterationTimings {
	double softmax;
	double point_smoothing;
	double suppression;
	double blur;
	double richardson_lucy;
	double convergence;
} IterationTimings;

static double IterationWallTime(void) {
#ifdef _OPENMP
	return omp_get_wtime();
#else
	struct timespec now;
	timespec_get(&now, TIME_UTC);
	return (double)now.tv_sec + (double)now.tv_nsec * 1e-9;
#endif
}


Decon ExitToBlank(const Config config, Decon decon)
{
	const int l = config.lengthmz * config.numz;
	decon.blur = calloc(l, sizeof(float));
	decon.newblur = calloc(l, sizeof(float));
	return decon;
}

void SetupInputs(const Config config, Input *inp, const int verbose){
	//.............................................................
	//
	//          Set the Mass Values to Exclude
	//
	//...............................................................

	// Check inputs
	if (verbose ==1 ) {printf("Checking inputs: %f %d\n", inp->dataMZ[0], config.lengthmz);}

	CalcMasses(&config, inp);

	//Allocates the memory for the boolean array of whether to test a value
	inp->barr = calloc(config.lengthmz * config.numz, sizeof(char));
	//Tells the algorithm to ignore data that are equal to zero
	ignorezeros(inp->barr, inp->dataInt, config.lengthmz, config.numz);
	if (verbose == 1) { printf("Ignored Zeros. Length: %d\n", config.lengthmz);  }

	//Sets limits based on mass range and any test masses
	SetLimits(config, inp);
	if (verbose == 1) { printf("Set limits\n"); }

	//Manual Assignments
	if (config.manualflag == 1)
	{
		ManualAssign(inp->dataMZ, inp->barr, inp->nztab, config);
		if (verbose == 1) { printf("Setup Manual Assign\n"); }
	}
}


int SetupDeconvolution(const Config config, const Input inp, Decon *decon, IntraDecon *intra, const int silent)
{
	// Set up memory and copy in binary array
	intra->ln = config.lengthmz * config.numz;
	intra->barr = calloc(intra->ln, sizeof(char));
	memcpy(intra->barr, inp.barr, (size_t)intra->ln * sizeof(char));

	//...................................................................
	//
	//     Sets the mzdist with the peak shape
	//
	//....................................................................

	SetUpPeakShape(config, inp, decon, silent, silent == 0);

	//....................................................
	//
	//    Setting up the neighborhood blur
	//
	//......................................................

	SetUpBlur(config, inp, intra, silent);

	//...................................................
	//
	//  Setting up and running the iteration
	//
	//.........................................................

	//Applies the intensity threshold to kill peaks
	if (config.intthresh != -1) { KillB(inp.dataInt, intra->barr, config.intthresh, config.lengthmz, config.numz); }

	//Creates an intial probability matrix, decon->blur, of 1 for each element
	decon->blur = calloc(intra->ln, sizeof(float));
	decon->newblur = calloc(intra->ln, sizeof(float));
	intra->oldblur = calloc(intra->ln, sizeof(float));
	intra->rl_deltas = malloc((size_t)config.lengthmz * sizeof(float));
	intra->rl_denom = malloc((size_t)config.lengthmz * sizeof(float));
	const int needs_grid_scratch = config.psig >= 1 || config.suppression_satellite > 0 ||
		config.suppression_harmonic > 0 || config.suppression_topn > 0 || config.suppression_topx > 0;
	if (needs_grid_scratch) {
		intra->smoothing_scratch = malloc((size_t)intra->ln * sizeof(float));
	}
	if (config.psig >= 1) {
		intra->smoothing_sums = malloc((size_t)config.numz * sizeof(float));
	}
	if (config.zsig >= 0 && config.msig >= 0 && intra->numclose > 1) {
		intra->log_blur = malloc((size_t)intra->ln * sizeof(float));
	}
	if (intra->rl_deltas == NULL || intra->rl_denom == NULL ||
		(needs_grid_scratch && intra->smoothing_scratch == NULL) ||
		(config.psig >= 1 && intra->smoothing_sums == NULL) ||
		(config.zsig >= 0 && config.msig >= 0 && intra->numclose > 1 && intra->log_blur == NULL)) {
		fprintf(stderr, "Error allocating iteration scratch arrays.\n");
		exit(11);
	}

	if (config.baselineflag == 1) {
		printf("Auto Baseline Mode On: %d\n", config.aggressiveflag);
		decon->baseline = calloc(config.lengthmz, sizeof(float));
		decon->noise = calloc(config.lengthmz, sizeof(float));
	}

	//#pragma omp parallel for private (i,j), schedule(auto)
	for (int i = 0; i < config.lengthmz; i++)
	{
		float val = inp.dataInt[i] / ((float)(config.numz + 2));
		if (config.baselineflag == 1) {
			decon->baseline[i] = val;
			decon->noise[i] = val;
		}

		for (int j = 0; j < config.numz; j++)
		{
			if (intra->barr[index2D(config.numz, i, j)] == 1) {
				decon->blur[index2D(config.numz, i, j)] = val;
			}
			else
			{
				decon->blur[index2D(config.numz, i, j)] = 0;
			}
		}
	}
	memcpy(intra->oldblur, decon->blur, sizeof(float) * config.lengthmz * config.numz);
	memcpy(decon->newblur, decon->blur, sizeof(float) * config.lengthmz * config.numz);


	// ............
	//
	// Run Baseline Subtraction if needed, otherwise just copy data over
	//
	// ..................

	intra->dataInt2 = calloc(config.lengthmz, sizeof(float));
	memcpy(intra->dataInt2, inp.dataInt, sizeof(float) * config.lengthmz);

	if (config.baselineflag == 1)
	{
		if (config.mzsig != 0)
		{
			deconvolve_baseline(config.lengthmz, inp.dataMZ, inp.dataInt, decon->baseline, fabsf(config.mzsig));
			if (config.aggressiveflag == 2)
			{
				for (int i = 0; i < 10; i++)
				{
					deconvolve_baseline(config.lengthmz, inp.dataMZ, inp.dataInt, decon->baseline, fabsf(config.mzsig));
				}
				for (int i = 0; i < config.lengthmz; i++)
				{
					if (decon->baseline[i] > 0) {
						intra->dataInt2[i] -= decon->baseline[i];
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

	//Determine the maximum intensity in the data for the betafactor
	float dmax = Max(intra->dataInt2, config.lengthmz);
	if (dmax > 1) { intra->betafactor=dmax; }

	// Final check that there are actually possible data in the set
	int badness = 1;
	for (int i = 0; i < config.lengthmz * config.numz; i++)
	{
		if (intra->barr[i] == 1) { badness = 0; break;}
	}
	return badness;
}

int CheckConvergence(const Config *config, Decon * decon, IntraDecon *intra, int iterations, const int silent) {
	float diff = 0;
	float tot = 0;
	#pragma omp parallel for reduction(+:diff,tot) schedule(static)
	for (int i = 0; i < intra->ln; i++)
	{
		if (intra->barr[i] == 1)
		{
			const float delta = decon->blur[i] - intra->oldblur[i];
			diff += delta * delta;
			tot += decon->blur[i];
		}
	}
	if (tot != 0) { decon->conv = (diff / tot); }
	else {
		if (decon->conv == 12345678) { printf("m/z vs. charge grid is zero. Iteration: %d\n", iterations); return 1; }
	}

	//printf("Iteration: %d Convergence: %f\n", iterations, decon.conv);
	if (decon->conv < 0.000001) {
		if (intra->off == 1 && config->numit > 0) {
			if (silent == 0) { printf("Converged in %d iterations.\n", iterations); }
			return 1;
		}
		intra->off = 1;
	}
	memcpy(intra->oldblur, decon->blur, (size_t) intra->ln * sizeof(float));
	return 0;
}

void RunIteration(const Config *config, Decon *decon, IntraDecon * intra, const Input *inp, int iterations,
	IterationTimings *timings) {
	decon->iterations = iterations;
	// if (config.minratio > 0){
	// 	adjust_ratios(config, barr, intra.numclose, closeind, decon.blur);
	// }

	// Softmax
	if (config->beta > 0 && iterations > 0)
	{
		const double start = timings != NULL ? IterationWallTime() : 0;
		softargmax(decon->blur, config->lengthmz, config->numz, config->beta/intra->betafactor);
		if (timings != NULL) { timings->softmax += IterationWallTime() - start; }
		//printf("Beta %f\n", beta);
		// softmax_peakwidth(config, decon, decon.blur, barr, config.beta/betafactor);
	}

	// Point Smoothing
	if (config->psig >= 1 && iterations > 0)
	{
		const double start = timings != NULL ? IterationWallTime() : 0;
		point_smoothing(decon->blur, intra->smoothing_scratch, intra->smoothing_sums, intra->barr,
			config->lengthmz, config->numz, abs((int)config->psig));
		if (timings != NULL) { timings->point_smoothing += IterationWallTime() - start; }
		//printf("Point Smoothed %f\n", config.psig);
	}

	// Suppression options
	const double suppression_start = timings != NULL ? IterationWallTime() : 0;
	if (iterations > config->suppression_startit &&
		(config->suppression_satellite > 0 || config->suppression_harmonic > 0 ||
		 config->suppression_topn > 0 || config->suppression_topx > 0)) {
		apply_suppressions(decon->blur, intra->smoothing_scratch, config->lengthmz, config->numz,
			config->suppression_satellite, config->suppression_harmonic, inp->nztab,
			config->suppression_topn, config->suppression_topx, config->suppression_percent);
	}
	if (timings != NULL) { timings->suppression += IterationWallTime() - suppression_start; }


	//Run Blurs
	const double blur_start = timings != NULL ? IterationWallTime() : 0;
	if (config->zsig >= 0 && config->msig >= 0) {
		blur_it_mean(intra, decon->newblur, decon->blur, config->zerolog);
	}
	else if (config->zsig > 0 && config->msig < 0)
	{
		blur_it_hybrid1(intra, config->lengthmz, config->numz, decon->newblur, decon->blur, config->zerolog);
	}
	else if (config->zsig < 0 && config->msig > 0)
	{
		blur_it_hybrid2(intra, config->lengthmz, config->numz, decon->newblur, decon->blur, config->zerolog);
	}
	else {
		blur_it(intra, decon->newblur, decon->blur);
	}
	if (timings != NULL) { timings->blur += IterationWallTime() - blur_start; }

	//Run Richardson-Lucy Deconvolution
	const double rl_start = timings != NULL ? IterationWallTime() : 0;
	deconvolve_iteration_speedy(config, decon, intra, inp->dataMZ);
	if (timings != NULL) { timings->richardson_lucy += IterationWallTime() - rl_start; }
}

void SetupOutputs(const Config config, Decon * decon, const IntraDecon intra, const Input inp, const int silent) {
		//................................................................
	//
	//     Setting up the outputs
	//
	//...............................................................

	//Reset the peak shape if it was inflated
	if (config.peakshapeinflate != 1 && config.mzsig != 0) {
		if (config.speedyflag == 0)
		{
			MakePeakShape2D(config, decon, &inp, 0, 0);
		}
		else
		{
			MakePeakShape1D(config, decon, inp.dataMZ, 0, 0, silent);
		}
		if (silent == 0) { printf("mzdist reset: %f\n", config.mzsig); }
	}


	//Determine the maximum intensity in the blur matrix
	float blurmax = Max(decon->blur, config.lengthmz * config.numz);
	float cutoff = 0;
	if (blurmax != 0) { cutoff = 0.000001f; }

	//Apply The Cutoff
	ApplyCutoff(decon->blur, blurmax * cutoff, config.lengthmz * config.numz);


	//Calculate the fit data and error.
	decon->fitdat = calloc(config.lengthmz, sizeof(float));
	decon->error = errfunspeedy(config, *decon, inp.dataInt, &decon->rsquared);

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
				decon->fitdat[i] = 0;
				decon->fitdat[i + 1] = 0;
			}
		}
	}


	// Charge scaling (orbimode)
	if (config.orbimode == 1)
	{
		printf("Rescaling charge states and normalizing ");
		charge_scaling(decon->blur, inp.nztab, config.lengthmz, config.numz);
		//simp_norm(config.lengthmz * config.numz, decon.blur);
		printf("Done\n");
	}

	//newblur is repurposed as the convolution of blur by the mz peaks shape
	float newblurmax = blurmax;
	if ((config.rawflag == 0 || config.rawflag == 2)) {
		if (config.mzsig != 0) {
			newblurmax = Reconvolve(config, decon, intra.barr);
		}
		else
		{
			memcpy(decon->newblur, decon->blur, (size_t)intra.ln * sizeof(float));
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
				if (decon->newblur[index2D(config.numz, i, j)] * (float) intra.barr[index2D(config.numz, i, j)] > newblurmax * cutoff)
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
	decon->mlen = (int)((massmax - massmin) / config.massbins);
	const int keep_mass_grid = config.rawflag == 0 || config.rawflag == 1;
	if (decon->mlen < 1) {
		printf("ERROR: No masses detected. Length: %d\n", decon->mlen);
		massmax = config.massub;
		massmin = config.masslb;
		decon->mlen = (int)((massmax - massmin) / config.massbins);

		//Declare the memory
		decon->massaxis = calloc(decon->mlen, sizeof(float));
		decon->massaxisval = calloc(decon->mlen, sizeof(float));
		if (keep_mass_grid) {
			decon->massgrid = calloc((size_t)decon->mlen * config.numz, sizeof(float));
		}

		//Create the mass axis
		for (int i = 0; i < decon->mlen; i++)
		{
			decon->massaxis[i] = massmin + (float) i * config.massbins;
		}
		decon->uniscore = 0;
	}
	else {

		//Declare the memory
		decon->massaxis = calloc(decon->mlen, sizeof(float));
		decon->massaxisval = calloc(decon->mlen, sizeof(float));
		if (keep_mass_grid) {
			decon->massgrid = calloc((size_t)decon->mlen * config.numz, sizeof(float));
		}
		if (silent == 0) { printf("Mass axis length: %d\n", decon->mlen); }


		//Create the mass axis
		for (int i = 0; i < decon->mlen; i++)
		{
			decon->massaxis[i] = massmin + (float) i * config.massbins;
		}

		//Determine the mass intensities from m/z grid
		if (config.poolflag == 0) {
			IntegrateTransform(config, decon, inp.mtab, massmax, massmin);
		}
		else if (config.poolflag == 1) {
			InterpolateTransform(config, decon, &inp);
		}
		else if (config.poolflag == 2) {
			SmartTransform(config, decon, &inp);
		}
		else {
			printf("Invalid poolflag %d\n", config.poolflag);
			exit(1987);
		}

		//.....................................
		// Scores
		// .......................................

		if (silent == 0 && keep_mass_grid) {
			//Note this will not execute if the mass axis is bad
			float scorethreshold = 0;
			decon->uniscore = score(config, decon, inp, scorethreshold, silent);
		}

		// Run DoubleDec if it's checked
		if (config.doubledec) { // Use Python to confirm kernel file is selected
			printf("Running DoubleDec now\n");
			DoubleDecon(&config, decon);
		}
	}
}


Decon MainDeconvolution(const Config config, const Input inp, const int silent)
{
	//Insert time
	time_t starttime;
	starttime = clock();

	Decon decon = InitDecon();
	IntraDecon intra = InitIntraDecon();

	int badness = SetupDeconvolution(config, inp, &decon, &intra, silent);

	// If badness found, abort
	if (badness == 1) {
		if (silent == 0) { printf("ERROR: Setup is bad. No points are allowed.\n Check that either mass smoothing, charge smoothing, manual assignment, or isotope mode are on.\n"); };
		decon = ExitToBlank(config, decon);
		FreeIntraDecon(intra);
		return(decon);
	}

	//Run the iteration

	// if (config.minratio > 0 && config.variablepw > 0){
	// 	check_ratios(config, inp, barr, intra->numclose, closeind, decon.blur);
	// }

	if (silent == 0) { printf("Iterating..."); }

	const int max_iterations = abs(config.numit);
	const int final_convergence_start = 9 * max_iterations / 10;
	IterationTimings timings = {0};
	IterationTimings *timing_ptr = silent == 0 ? &timings : NULL;
	for (int iterations = 0; iterations < max_iterations; iterations++)
	{
		// Run the iteration
		RunIteration(&config, &decon, &intra, &inp, iterations, timing_ptr);

		//Determine the metrics for conversion. Only do this every 10% to speed up.
		if (max_iterations < 10 || iterations == 1 || iterations % 10 == 0 ||
			iterations >= final_convergence_start) {
			const double convergence_start = timing_ptr != NULL ? IterationWallTime() : 0;
			int converged = CheckConvergence(&config, &decon, &intra, iterations, silent);
			if (timing_ptr != NULL) { timings.convergence += IterationWallTime() - convergence_start; }
			if (converged == 1) { break; }
		}
	}
	if (timing_ptr != NULL) {
		const double timed_total = timings.softmax + timings.point_smoothing + timings.suppression +
			timings.blur + timings.richardson_lucy + timings.convergence;
		printf("\nIteration Stage Times (s): Softmax %.6f  PointSmooth %.6f  Suppression %.6f  "
			"Blur %.6f  RichardsonLucy %.6f  Convergence %.6f  Total %.6f\n",
			timings.softmax, timings.point_smoothing, timings.suppression, timings.blur,
			timings.richardson_lucy, timings.convergence, timed_total);
	}

	SetupOutputs(config, &decon, intra, inp, silent);

	// Print out time
	if (silent == 0) { printf("Deconvolution Time: %f\n", (float)(clock() - starttime) / CLOCKS_PER_SEC); }

	//Free Memory
	FreeIntraDecon(intra);
	return decon;
}



int run_unidec_core(Config config, Input inp, Decon *decon, const int verbose) {
	// Setup the inputs
	SetupInputs(config, &inp, verbose);

	//Run the main Deconvolution
	*decon = MainDeconvolution(config, inp, config.silent);

	if (verbose == 1) { printf("Decon Done\n"); }
	return 0;
}


int run_unidec(int argc, char *argv[], Config config) {
	//Initialize
	clock_t starttime;
	starttime = clock();

	int verbose = 0;
	if (argc>2)
	{
		if (strcmp(argv[2], "-verbose") == 0) { verbose = 1; }
	}

	// File Inputs
	Input inp = InitInputs();
	ReadInputs(&config, &inp);

	if (verbose == 1) { printf("Read Inputs\n"); }

	// Setup the Deconvolution
	Decon decon=InitDecon();

	//  Run the Deconvolution
	run_unidec_core(config, inp, &decon, verbose);

	//Write Everything
	WriteDecon(config, &decon, &inp);
	if (verbose == 1) { printf("Write Done\n"); }

	clock_t end = clock();
	float totaltime = (float)(end - starttime) / CLOCKS_PER_SEC;
	WriteGlobalOutputs(config, decon, totaltime);

	//Free memory
	FreeInputs(inp);
	FreeDecon(decon);

	//printf("Error in the Fit: %f\n", decon.error);

	//Final Check that iterations worked and a reporter of the time consumed
	if (config.silent == 0) { printf("Finished with %d iterations in %f seconds!\n\n", decon.iterations, totaltime); }
	return 0;
}
