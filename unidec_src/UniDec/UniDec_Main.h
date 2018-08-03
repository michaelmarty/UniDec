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


int run_unidec(int argc, char *argv[], Config config) {
	//Initialize
	clock_t starttime;
	starttime = clock();

	FILE *out_ptr = NULL;
	hid_t file_id;

	unsigned int i, j, k, iterations;
	char *barr = NULL;

	int mlength, zlength, numclose,
		*nztab = NULL,
		*mind = NULL,
		*zind = NULL,
		*closemind = NULL,
		*closezind = NULL,
		*closeind=NULL,
		*starttab = NULL,
		*endtab = NULL,
		*isotopepos = NULL;
	double
		*mdist = NULL,
		*zdist = NULL,
		*mzdist = NULL,
		*mtab = NULL,
		*blur = NULL,
		*newblur = NULL,
		*oldblur = NULL,
		*fitdat = NULL,
		*massaxis = NULL,
		*massaxisval = NULL,
		*massgrid = NULL,
		*closeval = NULL,
		*dataMZ = NULL,
		*testmasses = NULL,
		*dataInt = NULL,
		*baseline=NULL,
		*noise=NULL;

	float *isotopeval = NULL;
	int lengthmz;
	int mfilelen = 0;
	char dataset[1024];
	char outdat[1024];


	//Set default parameters. These can be overwritten by config file.
	double killmass = 0;
	double isoparams[10] = { 1.00840852e+00, 1.25318718e-03, 2.37226341e+00, 8.19178000e-04,
		-4.37741951e-01, 6.64992972e-04, 9.94230511e-01, 4.64975237e-01, 1.00529041e-02, 5.81240792e-01 };

	if (argc>2)
	{
		if (strcmp(argv[2], "-kill") == 0)
		{
			killmass = atof(argv[3]);
			printf("Killing Peak: %f\n", killmass);
		}
		if (strcmp(argv[2], "-test") == 0)
		{
			killmass = atof(argv[3]);
			printf("Testing Mass: %f\n", killmass);
			test_isotopes(killmass, isoparams);
			exit(0);
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
			lengthmz = mh5getfilelength(file_id, outdat);
			dataMZ = calloc(lengthmz, sizeof(double));
			dataInt = calloc(lengthmz, sizeof(double));
			mh5readfile2d(file_id, outdat, lengthmz, dataMZ, dataInt);
			printf("Length of Data: %d \n", lengthmz);

			//Check the length of the mfile and then read it in.
			if (config.mflag == 1)
			{
				mfilelen = mh5getfilelength(file_id, "/config/masslist");
				printf("Length of mfile: %d \n", mfilelen);
				testmasses = malloc(sizeof(double)*mfilelen);
				mh5readfile1d(file_id, "/config/masslist", testmasses);
			}
			else {
				testmasses = malloc(sizeof(double)*mfilelen);
			}
		}
		else {

			//Calculate the length of the data file automatically
			lengthmz = getfilelength(config.infile);
			dataMZ = calloc(lengthmz, sizeof(double));
			dataInt = calloc(lengthmz, sizeof(double));

			readfile(config.infile, lengthmz, dataMZ, dataInt);//load up the data array

			printf("Length of Data: %d \n", lengthmz);

			//Check the length of the mfile and then read it in.

			if (config.mflag == 1)
			{
				mfilelen = getfilelength(config.mfile);
				printf("Length of mfile: %d \n", mfilelen);
			}
			testmasses = malloc(sizeof(double)*mfilelen);
			if (config.mflag == 1)
			{
				readmfile(config.mfile, mfilelen, testmasses);//read in mass tab
			}
		}
	}


	//This for loop creates a list of charge values
	nztab = calloc(config.numz, sizeof(int));
	for (i = 0; i<config.numz; i++) { nztab[i] = i + config.startz; }
	//printf("nzstart %d\n",nztab[0]);

	//.............................................................
	//
	//          Set the Mass and Intensity values
	//
	//...............................................................

	//Initializes mtab and barr matrices to the same dimensions of length mz by number of charges.
	//mztab = calloc(lengthmz*config.numz, sizeof(double));
	mtab = calloc(lengthmz*config.numz, sizeof(double));
	barr = calloc(lengthmz*config.numz, sizeof(char));

	//Fills mztab from dataMZ, mtab by multiplying each mz by each z
#pragma omp parallel for private (i,j), schedule(dynamic)
	for (i = 0; i<lengthmz; i++)
	{
		for (j = 0; j<config.numz; j++)
		{
			//mztab[index2D(config.numz, i, j)] = dataMZ[i];
			mtab[index2D(config.numz, i, j)] = (dataMZ[i] * nztab[j] - config.adductmass*nztab[j]);
		}
	}

	//Determines the indexes of each test mass from mfile in m/z space
	int *testmasspos = malloc(sizeof(double)*mfilelen*config.numz);
	if (config.mflag == 1 && config.limitflag == 1) {
		for (i = 0; i<mfilelen; i++)
		{
			for (j = 0; j<config.numz; j++)
			{
				double mztest = (testmasses[i] + config.adductmass*nztab[j]) / (double)nztab[j];
				testmasspos[index2D(config.numz, i, j)] = nearfast(dataMZ, mztest, lengthmz);
			}
		}
	}

	//If there is a mass file read, it will only allow masses close to those masses within some config.mtabsig window.
	if (config.mflag == 1 && config.limitflag == 0)
	{
		TestMassListWindowed(lengthmz, config.numz, barr, mtab, config.nativezub, config.nativezlb, config.massub, config.masslb, nztab, testmasses, mfilelen, config.mtabsig);
	}
	//If there is a mass file read and the mass table window (config.mtabsig) is 0, it will only write intensities at the m/z values closest to the m/z values read in from the mfile.
	else if (config.mflag == 1 && config.limitflag == 1)
	{
		TestMassListLimit(lengthmz, config.numz, barr, mtab, config.nativezub, config.nativezlb, config.massub, config.masslb, nztab, testmasspos, mfilelen);
	}
	//Normally, write the intensity values if the values fall within the mass upperbound and lower bound
	else
	{
		TestMass(lengthmz, config.numz, barr, mtab, config.nativezub, config.nativezlb, config.massub, config.masslb, nztab);
	}

	if (config.manualflag == 1)
	{
		ManualAssign(config.manualfile, lengthmz, config.numz, dataMZ, barr, nztab, file_id, config);
	}

	if (config.killmass != 0 && config.mzsig!=0)
	{
		KillMass(config.killmass, lengthmz, config.numz, barr, nztab, config.adductmass, dataMZ, config.psfun, config.mzsig);
	}

	//Get max intensity
	double maxint = 0;
	for (i = 0; i<lengthmz; i++)
	{
		if (dataInt[i]>maxint) { maxint = dataInt[i]; }
	}


	//...................................................................
	//
	//     Sets the mzdist with the peak shape
	//
	//....................................................................


	//Sets a threshold for m/z values to check. Things that are far away in m/z space don't need to be considered in the iterations.
	double threshold = config.psthresh*config.mzsig*config.peakshapeinflate;
	printf("Threshold: %f\n", threshold);
	//Create a list of start and end values to box in arrays based on the above threshold
	starttab = calloc(lengthmz, sizeof(int));
	endtab = calloc(lengthmz, sizeof(int));
	int maxlength = 1;
	if(config.mzsig!=0){
		for (i = 0; i<lengthmz; i++)
		{
			double point = dataMZ[i] - threshold;
			int start, end;
			if (point < dataMZ[0] && config.speedyflag == 0) { start = (int)((point - dataMZ[0]) / (dataMZ[1] - dataMZ[0])); }
			else {
				start = nearfast(dataMZ, point, lengthmz);
			}

			starttab[i] = start;
			point = dataMZ[i] + threshold;
			if (point > dataMZ[lengthmz - 1] && config.speedyflag == 0) { end = lengthmz - 1 + (int)((point - dataMZ[lengthmz - 1]) / (dataMZ[lengthmz - 1] - dataMZ[lengthmz - 2])); }
			else {
				end = nearfast(dataMZ, point, lengthmz);
			}
			endtab[i] = end;
			if (end - start>maxlength) { maxlength = end - start; }
		}
		printf("maxlength %d\n", maxlength);


		//Changes dimensions of the peak shape function. 1D for speedy and 2D otherwise
		int pslen = lengthmz;
		if (config.speedyflag == 0) { pslen = lengthmz*maxlength; }
		mzdist = calloc(pslen, sizeof(double));
		memset(mzdist, 0, pslen * sizeof(double));


		//Calculates the distance between mz values as a 2D or 3D matrix
		if (config.speedyflag == 0)
		{
			MakePeakShape2D(lengthmz, maxlength, starttab, endtab, dataMZ, config.mzsig*config.peakshapeinflate, config.psfun, config.speedyflag, mzdist);
		}
		else
		{
			//Calculates peak shape as a 1D list centered at the first element for circular convolutions
			MakePeakShape1D(dataMZ, threshold, lengthmz, config.speedyflag, config.mzsig*config.peakshapeinflate, config.psfun, mzdist);
		}
		printf("mzdist set: %f\n", mzdist[0]);
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
	int chargewidth = 0;
	if (config.zsig > 1) { chargewidth = config.zsig; config.zsig = 1; }
	if (config.msig > 1) { chargewidth = config.msig; config.msig = 1; }
	if (config.zsig >= 0) {
		zlength = 1 + 2 * (int)config.zsig;
		mlength = 1 + 2 * (int)config.msig;
	}
	else {
		zlength = 1 + 2 * (int)(3 * fabs(config.zsig) + 0.5);
		mlength = 1 + 2 * (int)(3 * fabs(config.msig) + 0.5);
	}
	numclose = mlength*zlength;

	//Sets up the blur function in oligomer mass and charge
	mind = calloc(mlength, sizeof(int));
	mdist = calloc(mlength, sizeof(double));
	for (i = 0; i<mlength; i++)
	{
		mind[i] = i - (mlength - 1) / 2;
		mdist[i] = exp(-(pow((i - (mlength - 1) / 2.), 2)) / (2.0 * config.msig*config.msig));
		//mdist[i]=secderndis((mlength-1)/2.,config.msig,i);
	}

	zind = calloc(zlength, sizeof(int));
	zdist = calloc(zlength, sizeof(double));
	for (i = 0; i<zlength; i++)
	{
		zind[i] = i - (zlength - 1) / 2;
		zdist[i] = exp(-(pow((i - (zlength - 1) / 2.), 2)) / (2.0 * config.zsig*config.zsig));
		//zdist[i]=secderndis((zlength-1)/2.,zsig,i);
	}

	//Initializing memory
	closemind = calloc(numclose, sizeof(int));
	closezind = calloc(numclose, sizeof(int));
	closeval = calloc(numclose, sizeof(double));
	closeind = calloc(numclose*lengthmz*config.numz, sizeof(double));

	//Determines the indexes of things that are close as well as the values used in the neighborhood convolution
	double norm = 0;
	for (k = 0; k<numclose; k++)
	{
		norm += zdist[(int)k / mlength] * mdist[k%mlength];
	}

	for (k = 0; k<numclose; k++)
	{
		closemind[k] = mind[k%mlength];
		closezind[k] = zind[(int)k / mlength];
		closeval[k] = zdist[(int)k / mlength] * mdist[k%mlength] / norm;
	}

	//Set up blur
	MakeBlur(lengthmz, config.numz, numclose, barr, closezind, closemind, mtab, config.molig, config.adductmass, nztab, dataMZ, closeind);
	

	printf("Number of Oligomers blurred: %d\n", mlength);
	printf("Number of Charges blurred: %d\n", zlength);

	int isolength = 0;
	if (config.isotopemode == 1)
	{
		isolength = setup_isotopes(isoparams, isotopepos, isotopeval, mtab, nztab, barr, dataMZ, lengthmz, config.numz);

		isotopepos = calloc(isolength*lengthmz*config.numz, sizeof(int));
		isotopeval = calloc(isolength*lengthmz*config.numz, sizeof(float));

		make_isotopes(isoparams, isotopepos, isotopeval, mtab, nztab, barr, dataMZ, lengthmz, config.numz);

		printf("Isotopes set up, Length: %d\n", isolength);
	}

	//...................................................
	//
	//  Setting up and running the iteration
	//
	//.........................................................


	//Convert aggressiveflag to baselineflag
	if (config.aggressiveflag == 1 || config.aggressiveflag==2) { config.baselineflag = 1; }
	else { config.baselineflag = 0; }


	//Creates an intial probability matrix, blur, of 1 for each element
	blur = calloc(lengthmz*config.numz, sizeof(double));
	newblur = calloc(lengthmz*config.numz, sizeof(double));
	oldblur = calloc(lengthmz*config.numz, sizeof(double));

	if (config.baselineflag == 1) {
		printf("Auto Baseline Mode On: %d\n",config.aggressiveflag);
		baseline = calloc(lengthmz, sizeof(double));
		noise = calloc(lengthmz, sizeof(double));
	}

	//#pragma omp parallel for private (i,j), schedule(dynamic)
	for (i = 0; i<lengthmz; i++)
	{
		double val = dataInt[i] / ((float)(config.numz+2));
		if (config.baselineflag==1){
			
			baseline[i] = val;
			noise[i] = val;
		}

		for (j = 0; j<config.numz; j++)
		{
			if (barr[index2D(config.numz, i, j)] == 1) {
				if (config.isotopemode == 0) {
					blur[index2D(config.numz, i, j)] = val ;
				}
				else { blur[index2D(config.numz, i, j)] = 1; }
			}
			else
			{
				blur[index2D(config.numz, i, j)] = 0;
			}
		}
	}
	memcpy(oldblur, blur, sizeof(double)*lengthmz*config.numz);
	memcpy(newblur, blur, sizeof(double)*lengthmz*config.numz);

	if (config.intthresh != 0) { KillB(dataInt, barr, config.intthresh, lengthmz, config.numz, isolength, isotopepos, isotopeval); }

	
	double *dataInt2 = NULL;
	dataInt2 = calloc(lengthmz, sizeof(double));
	memcpy(dataInt2, dataInt, sizeof(double)*lengthmz);
	if (config.baselineflag == 1)
	{
		if(config.mzsig!=0)
		{
			deconvolve_baseline(lengthmz, dataMZ, dataInt, baseline, config.mzsig);
			if (config.aggressiveflag == 2)
			{
				for (int i = 0; i < 10; i++)
				{
					deconvolve_baseline(lengthmz, dataMZ, dataInt, baseline, config.mzsig);
				}
				for (int i = 0; i < lengthmz; i++)
				{
					if (baseline[i]>0) {
						dataInt2[i] -= baseline[i];
					}
				}
				//memcpy(baseline, dataInt2, sizeof(double)*lengthmz);
			}
		}
		else
		{
			printf("Ignoring baseline subtraction because peak width is 0\n");
		}
		
	}
	

	//Run the iteration
	double blurmax = 0, conv = 0;
	int off = 0;
	printf("Iterating\n\n");
	for (iterations = 0; iterations<abs(config.numit); iterations++)
	{
		if (chargewidth > 0 && iterations>20)
		{
			charge_smoothing(blur, lengthmz, config.numz, chargewidth - 1);
			//printf("Charge Smoothed %d\n", chargewidth - 1);
		}

		if (config.zsig >= 0) {
			blur_it_mean(lengthmz,
				config.numz,
				numclose,
				closeind,
				newblur,
				blur,
				barr,
				config.zerolog);
		}
		else {
			blur_it(lengthmz,
				config.numz,
				numclose,
				closeind,
				closeval,
				newblur,
				blur,
				barr);
		}
		
				
		deconvolve_iteration_speedy(lengthmz, config.numz, maxlength,
			newblur, blur, barr, config.aggressiveflag, dataInt2,
			isolength, isotopepos, isotopeval, starttab, endtab, mzdist, config.speedyflag,
			config.baselineflag, baseline,noise,config.mzsig,dataMZ, config.filterwidth);
		
		

		//Determine the metrics for conversion. Only do this every 10% to speed up.
		if ((config.numit<10 || iterations % 10 == 0 || iterations % 10 == 1 || iterations>0.9*config.numit)) {
			double diff = 0;
			double tot = 0;
			for (i = 0; i<lengthmz*config.numz; i++)
			{
				if (barr[i] == 1)
				{
					diff += pow((blur[i] - oldblur[i]), 2);
					tot += blur[i];
				}
			}
			conv = (diff / tot);

			printf("Iteration: %d Convergence: %f\n", iterations, conv);
			if (conv<0.000001) {
				if (off == 1 && config.numit>0) {
					printf("\nConverged in %d iterations.\n\n", iterations);
					break;
				}
				off = 1;
			}
			memcpy(oldblur, blur, lengthmz*config.numz * sizeof(double));
		}
	}
	free(dataInt2);
	//................................................................
	//
	//     Writing and reporting the outputs
	//
	//...............................................................

	//Reset the peak shape if it was inflated
	if (config.peakshapeinflate != 1 && config.mzsig!=0) {
		if (config.speedyflag == 0)
		{
			MakePeakShape2D(lengthmz, maxlength, starttab, endtab, dataMZ, config.mzsig, config.psfun, config.speedyflag, mzdist);
		}

		else
		{
			MakePeakShape1D(dataMZ, threshold, lengthmz, config.speedyflag, config.mzsig, config.psfun, mzdist);
		}
		printf("mzdist reset: %f\n", config.mzsig);
	}


	//Determine the maximum intensity in the blur matrix
	blurmax = Max(blur, lengthmz * config.numz);
	config.cutoff = 0.0000001 / blurmax;
	//Apply The Cutoff
	ApplyCutoff1D(blur, blurmax*config.cutoff, lengthmz*config.numz);


	//Calculate the fit data and error.
	fitdat = calloc(lengthmz, sizeof(double));
	double error;
	error = errfunspeedy(dataInt, fitdat,blur,lengthmz, config.numz, maxlength, maxint,
		isolength, isotopepos, isotopeval, starttab, endtab, mzdist, config.speedyflag);

	
	if (config.baselineflag == 1)
	{
		#pragma omp parallel for private(i), schedule(dynamic)
		for (int i = 0; i<lengthmz; i++) {
			fitdat[i] += baseline[i];// +noise[i];
			//fitdat[i] = noise[i]+0.1;
		}
	}
	

	ApplyCutoff1D(fitdat, 0, lengthmz);

	//Write the fit data to a file.
	if (config.rawflag >= 0) {
		if (config.filetype == 0) {
			char outstring2[500];
			sprintf(outstring2, "%s_fitdat.bin", config.outfile);
			out_ptr = fopen(outstring2, "wb");
			fwrite(fitdat, sizeof(double), lengthmz, out_ptr);
			fclose(out_ptr);
			printf("Fit data written to: %s\n", outstring2);
		}
		else {
			strjoin(dataset, "/fit_data", outdat);
			mh5writefile1d(file_id, outdat, lengthmz, fitdat);
		}
	}

	//Write the baseline to a file.
	if (config.rawflag >= 0 && config.baselineflag==1) {
		if (config.filetype == 0) {
			char outstring2[500];
			sprintf(outstring2, "%s_baseline.bin", config.outfile);
			out_ptr = fopen(outstring2, "wb");
			fwrite(baseline, sizeof(double), lengthmz, out_ptr);
			fclose(out_ptr);
			printf("Background written to: %s\n", outstring2);
		}
		else {
			strjoin(dataset, "/baseline", outdat);
			mh5writefile1d(file_id, outdat, lengthmz, baseline);
		}
	}

	// Charge scaling (orbimode)
	if (config.orbimode == 1)
	{
		printf("Rescaling charge states and normalizing ");
		charge_scaling(blur, nztab, lengthmz, config.numz);
		simp_norm(lengthmz*config.numz,blur);
		printf("Done\n");
	}

	//newblur is repurposed as the convolution of blur by the mz peaks shape
	double newblurmax = blurmax;
	if ((config.rawflag == 0 || config.rawflag == 2)) {
		if (config.mzsig != 0) {
			newblurmax = Reconvolve(lengthmz, config.numz, maxlength, starttab, endtab, mzdist, blur, newblur, config.speedyflag);
		}
		else
		{
			memcpy(newblur, blur, lengthmz*config.numz * sizeof(double));
		}
	}

	//Writes the convolved m/z grid in binary format
	if (config.rawflag == 0 || config.rawflag == 1)
	{
		if (config.filetype == 0) {
			char outstring9[500];
			sprintf(outstring9, "%s_grid.bin", config.outfile);
			out_ptr = fopen(outstring9, "wb");
			if (config.rawflag == 0) { fwrite(newblur, sizeof(double), lengthmz*config.numz, out_ptr); }
			if (config.rawflag == 1) { fwrite(blur, sizeof(double), lengthmz*config.numz, out_ptr); }
			fclose(out_ptr);
			printf("m/z grid written to: %s\n", outstring9);
		}
		else {
			strjoin(dataset, "/mz_grid", outdat);
			if (config.rawflag == 0) { mh5writefile1d(file_id, outdat, lengthmz*config.numz, newblur); }
			if (config.rawflag == 1) { mh5writefile1d(file_id, outdat, lengthmz*config.numz, blur); }

			strjoin(dataset, "/mz_grid", outdat);
			double *chargedat = NULL;
			chargedat = calloc(config.numz, sizeof(double));
			double *chargeaxis = NULL;
			chargeaxis = calloc(config.numz, sizeof(double));
			
			for (int j = 0; j < config.numz; j++) {
				double val = 0;
				chargeaxis[j] = (double)nztab[j];
				for (int i = 0; i<lengthmz; i++) {
				
					val += newblur[index2D(config.numz, i, j)];
				}
				chargedat[j] = val;
			}
			strjoin(dataset, "/charge_data", outdat);
			mh5writefile2d(file_id, outdat, config.numz, chargeaxis, chargedat);
			free(chargedat);
			free(chargeaxis);
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
		for (i = 0; i<lengthmz; i++)
		{
			for (j = 0; j<config.numz; j++)
			{
				if (newblur[index2D(config.numz, i, j)] * barr[index2D(config.numz, i, j)]>newblurmax*config.cutoff)
				{
					if (mtab[index2D(config.numz, i, j)] + threshold*nztab[j]>massmax)
					{
						massmax = mtab[index2D(config.numz, i, j)] + threshold*nztab[j];
					}
					if (mtab[index2D(config.numz, i, j)] - threshold*nztab[j]<massmin)
					{
						massmin = mtab[index2D(config.numz, i, j)] - threshold*nztab[j];
					}
				}
			}
		}
		printf("Massmax: %f  ", massmax);
		printf("Massmin: %f  ", massmin);
	}
	else { massmax = config.massub; massmin = config.masslb; }

	//Performs an interpolation to get the zero charge mass values starting at config.masslb and going to config.massub in steps of config.massbins
	int maaxle = (int)(massmax - massmin) / config.massbins;
	if (maaxle<1) {
		printf("Bad mass axis length: %d\n", maaxle);
		massmax = config.massub;
		massmin = config.masslb;
		maaxle = (int)(massmax - massmin) / config.massbins;
	}
	else {

		//Declare the memory

		massaxis = calloc(maaxle, sizeof(double));
		massaxisval = calloc(maaxle, sizeof(double));
		massgrid = calloc(maaxle*config.numz, sizeof(double));
		memset(massaxisval, 0, maaxle * sizeof(double));
		memset(massgrid, 0, maaxle*config.numz * sizeof(double));
		printf("Mass axis length: %d\n", maaxle);

		//To prevent really wierd decimals
		//Only do this if the axis isn't already fixed
		if (config.fixedmassaxis == 0) { massmin = round(massmin / config.massbins)*config.massbins; }

		//Create the mass axis
		for (i = 0; i < maaxle; i++)
		{
			massaxis[i] = massmin + i*config.massbins;
		}

		//Determine the mass intensities from m/z grid
		if (config.poolflag == 0) {
			if (config.rawflag == 1 || config.rawflag == 3) {
				IntegrateTransform(lengthmz, config.numz, mtab, massmax, massmin, maaxle, massaxis, massaxisval, blur, massgrid);
			}
			if (config.rawflag == 0 || config.rawflag == 2) {
				IntegrateTransform(lengthmz, config.numz, mtab, massmax, massmin, maaxle, massaxis, massaxisval, newblur, massgrid);
			}
		}
		else {
			if (config.rawflag == 1 || config.rawflag == 3) {
				InterpolateTransform(maaxle, config.numz, lengthmz, nztab, massaxis, config.adductmass, dataMZ, dataInt, massgrid, massaxisval, blur);
			}
			if (config.rawflag == 0 || config.rawflag == 2) {
				InterpolateTransform(maaxle, config.numz, lengthmz, nztab, massaxis, config.adductmass, dataMZ, dataInt, massgrid, massaxisval, newblur);
			}
		}


		//Writes the convolved mass grid in binary format
		if (config.rawflag == 0 || config.rawflag == 1) {
			if (config.filetype == 0) {
				char outstring10[500];
				sprintf(outstring10, "%s_massgrid.bin", config.outfile);
				out_ptr = fopen(outstring10, "wb");
				fwrite(massgrid, sizeof(double), maaxle*config.numz, out_ptr);
				fclose(out_ptr);
				printf("Mass Grid Written: %s\n", outstring10);
			}
			else {
				strjoin(dataset, "/mass_grid", outdat);
				mh5writefile1d(file_id, outdat, maaxle*config.numz, massgrid);
			}
		}

		//Writes the mass values convolved with the peak shape
		if (config.rawflag == 0 || config.rawflag == 1 || config.rawflag == 2 || config.rawflag == 3) {
			if (config.filetype == 0) {
				char outstring4[500];
				sprintf(outstring4, "%s_mass.txt", config.outfile);
				out_ptr = fopen(outstring4, "w");
				for (i = 0; i < maaxle; i++)
				{
					fprintf(out_ptr, "%f %f\n", massaxis[i], massaxisval[i]);
				}
				fclose(out_ptr);
				printf("Masses written: %s\n", outstring4);
			}
			else {
				strjoin(dataset, "/mass_data", outdat);
				mh5writefile2d(file_id, outdat, maaxle, massaxis, massaxisval);
			}
		}
	}


	//................................................................
	//
	//  Wrapping up
	//
	//...................................................................

	//Writes a file with a number of key parameters such as error and number of significant parameters.
	clock_t end = clock();
	double totaltime = (double)(end - starttime) / CLOCKS_PER_SEC;
	if (config.filetype == 0) {
		char outstring3[500];
		sprintf(outstring3, "%s_error.txt", config.outfile);
		out_ptr = fopen(outstring3, "w");
		fprintf(out_ptr, "error = %f\n", error);
		fprintf(out_ptr, "blurmax = %f\n", blurmax);
		fprintf(out_ptr, "newblurmax = %f\n", newblurmax);
		fprintf(out_ptr, "cutoff = %f\n", config.cutoff);
		fprintf(out_ptr, "time = %f\n", totaltime);
		fprintf(out_ptr, "interations = %d\n", iterations);
		fclose(out_ptr);
		printf("Stats and Error written to: %s\n", outstring3);
	}
	else {
		write_attr_double(file_id, dataset, "error", error);
		write_attr_int(file_id, dataset, "interations", iterations);
		write_attr_double(file_id, dataset, "time", totaltime);
		write_attr_double(file_id, dataset, "blurmax", blurmax);
		write_attr_double(file_id, dataset, "newblurmax", newblurmax);
		write_attr_double(file_id, dataset, "cutoff", config.cutoff);
		write_attr_int(file_id, dataset, "length_mz", lengthmz);
		write_attr_int(file_id, dataset, "length_mass", maaxle);
		set_needs_grids(file_id);
		H5Fclose(file_id);
	}

	//Free memory
	free(mtab);
	free(mzdist);;
	free(closeval);
	free(closemind);
	free(closezind);
	free(fitdat);
	free(blur);
	free(endtab);
	free(dataMZ);
	free(massaxisval);
	free(dataInt);
	free(starttab);
	free(massaxis);
	free(newblur);
	free(oldblur);
	free(mdist);
	free(mind);
	free(zind);
	free(zdist);
	free(nztab);
	free(barr);
	free(massgrid);
	free(isotopepos);
	free(isotopeval);
	free(testmasses);
	free(testmasspos);
	free(baseline);
	free(noise);
	free(closeind);

	printf("\nError: %f\n", error);

	//Final Check that iterations worked and a reporter of the time consumed
	printf("\nFinished with %d iterations and %f convergence in %f seconds!\n\n", iterations, conv, totaltime);
	return 0;
}