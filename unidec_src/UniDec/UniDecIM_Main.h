/*
* UniDecIM_Main.h
*
*  Created on : 23 December 2016
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
#include <fftw3.h>
#include "UniDecIM.h"


int run_unidec_IM(int argc, char *argv[], Config config) {

	time_t starttime, endtime;
	starttime = time(NULL);

	//printf("config.tcal1 %f config.tcal2 %f config.edc %f gmass %f\n", config.tcal1, config.tcal2, config.edc, config.hmass);
	//printf("dt calc %f\n", calcDtTwaveLog(98000., 10, 5000., config.tcal1, config.tcal2, config.hmass, config.edc));a

	printf("Opening File: %s\n", config.infile);
	int lines;
	lines = getfilelength(config.infile);

	printf("config.length of data: %d\n", lines);

	int i, j, k, l, m, tmlen, manlength;

	int *ztab = NULL, *closemind = NULL,
		*closezind = NULL,
		*closecind = NULL,
		*barr = NULL, *closetab = NULL;
	int size[4];

	double *mzdat = NULL,
		*dtdat = NULL,
		*dataInt = NULL,
		*mzext = NULL,
		*dtext = NULL,
		*masstab = NULL,
		*ccstab = NULL,
		*blur = NULL,
		*newblur = NULL,
		*peakshape = NULL,
		*denom = NULL,
		*deltas = NULL,
		*testmasses = NULL,
		*testmassvals = NULL,
		*testmassCCSavg = NULL,
		*testmassCCSstd = NULL,
		*testmassZavg = NULL,
		*testmassZstd = NULL
		;

	//Mass Limit File
	if (config.mflag == 1) {
		tmlen = getfilelength(config.mfile);
		testmasses = calloc(tmlen, sizeof(double));
		testmassvals = calloc(tmlen, sizeof(double));
		testmassCCSavg = calloc(tmlen, sizeof(double));
		testmassCCSstd = calloc(tmlen, sizeof(double));
		testmassZavg = calloc(tmlen, sizeof(double));
		testmassZstd = calloc(tmlen, sizeof(double));
		readmfile(config.mfile, tmlen, testmasses);
		printf("Read mass file of length: %d\n", tmlen);
	}

	//CCS Constants
	double e = 1.60217657E-19;
	double kb = 1.3806488E-23;
	double n = 2.6867774E25;
	double po = 760;
	double tempo = 273.15;
	double pi = 3.14159265359;
	config.temp = config.temp + tempo;
	double ccsconst = (sqrt(18 * pi) / 16)*(e / sqrt(kb*config.temp)) / n*(config.volt / pow(config.len, 2))*(po / config.press)*(config.temp / tempo)*1E20;
	if (config.twaveflag>0) { printf("Ridin' the T-Wave!\n"); }

	//Reading In Data
	mzdat = calloc(lines, sizeof(double));
	dtdat = calloc(lines, sizeof(double));
	dataInt = calloc(lines, sizeof(double));
	readfile3(config.infile, lines, mzdat, dtdat, dataInt);

	//Charge States
	int numz = config.endz - config.startz + 1;
	ztab = calloc(numz, sizeof(int));
	for (i = 0; i<numz; i++)
	{
		ztab[i] = config.startz + i;
	}

	//Getting Dimension Sizes of Data
	size[0] = GetSize0(mzdat, lines);
	size[1] = GetSize1(dtdat, lines);
	size[2] = numz;
	int totlen = size[0] * size[1] * size[2];
	printf("Dimensions of data: %d mz by %d dt by %d z: %d lines: %d total\n", size[0], size[1], size[2], size[0] * size[1], totlen);
	if (totlen>155E6) { printf("Warning: May exceed system memory capacity"); }

	//Extracting mz and dt ranges
	mzext = calloc(size[0], sizeof(double));
	dtext = calloc(size[1], sizeof(double));
	Extract(mzext, dtext, mzdat, dtdat, size);
	double mzranges[4];
	mzranges[0] = mzext[0];
	mzranges[1] = mzext[size[0] - 1];
	mzranges[2] = dtext[0];
	mzranges[3] = dtext[size[1] - 1];
	printf("MZ Range: %f config.to %f\n", mzranges[0], mzranges[1]);
	printf("DT Range: %f config.to %f\n", mzranges[2], mzranges[3]);

	peakshape = calloc(lines, sizeof(double));
	GetPeaks(peakshape, size, mzext, dtext, config.mzsig, config.dtsig, config.psfun);
	printf("Peak Shape Set\n");

	//Filling the mass table
	ccstab = calloc(totlen, sizeof(double));
	masstab = calloc(size[0] * size[2], sizeof(double));
	barr = calloc(totlen, sizeof(int));
#pragma omp parallel for private (i,j), schedule(dynamic)
	for (i = 0; i<size[0]; i++)
	{
		for (j = 0; j<size[2]; j++)
		{
			masstab[index2D(size[2], i, j)] = calcMass(mzext[i], ztab[j], config.adductmass);
		}
	}
	//Fill CCS table and simultaneously set limit array


	double tempmass, tempccs, testmassclose, zlimit, climit;
#pragma omp parallel for private (i,j,k,tempmass,tempccs,testmassclose,zlimit,climit), schedule(dynamic)
	for (i = 0; i<size[0]; i++)
	{
		for (j = 0; j<size[1]; j++)
		{
			for (k = 0; k<size[2]; k++)
			{
				tempmass = masstab[index2D(size[2], i, k)];
				if (config.twaveflag == 0)
				{
					tempccs = calcCCS(masstab[index2D(size[2], i, k)], ztab[k], dtext[j], ccsconst, config.hmass, config.to);
				}
				else if (config.twaveflag == 1)
				{
					tempccs = calcCCSTwaveLog(masstab[index2D(size[2], i, k)], ztab[k], dtext[j], config.tcal1, config.tcal2, config.hmass, config.edc);
				}
				else if (config.twaveflag == 2)
				{
					tempccs = calcCCSTwaveLinear(masstab[index2D(size[2], i, k)], ztab[k], dtext[j], config.tcal1, config.tcal2, config.hmass, config.edc);
				}
				else if (config.twaveflag == 3)
				{
					tempccs = calcCCSTwavePower(masstab[index2D(size[2], i, k)], ztab[k], dtext[j], config.tcal1, config.tcal2, config.hmass, config.edc);
				}
				else
				{
					TwaveError(config.twaveflag);
				}
				ccstab[index3D(size[1], size[2], i, j, k)] = tempccs;
				zlimit = nativecharge(tempmass, 0);
				climit = nativeCCS(tempmass, 0, config.hmass);
				if (tempmass<config.massub&&tempmass>config.masslb&&tempccs<config.ccsub&&tempccs>config.ccslb&&ztab[k]<zlimit + config.nativezub&&ztab[k]>zlimit + config.nativezlb&&tempccs<climit + config.nativeccsub&&tempccs>climit + config.nativeccslb)
				{
					if (config.mflag == 1)
					{
						testmassclose = nearfastpoint(testmasses, tempmass, tmlen);
						if (fabs(testmassclose - tempmass) <= config.mtabsig) { barr[index3D(size[1], size[2], i, j, k)] = 1; }
						else { barr[index3D(size[1], size[2], i, j, k)] = 0; }
					}
					else
					{
						barr[index3D(size[1], size[2], i, j, k)] = 1;
					}
				}
				else
				{
					barr[index3D(size[1], size[2], i, j, k)] = 0;
				}
			}
		}
	}

	if (config.manualflag == 1)
	{
		ManualAssign_IM(config.manualfile, size, mzdat, dtdat, ztab, barr);
	}


	//Working out the Blur
	int mlength = 2 * config.msig + 1;
	int zlength = 2 * config.zsig + 1;
	int clength = 2 * config.csig + 1;
	if (config.csig < 0) { clength = 1; }
	int numclose = mlength*zlength*clength;
	size[3] = numclose;
	closemind = calloc(numclose, sizeof(int));
	closezind = calloc(numclose, sizeof(int));
	closecind = calloc(numclose, sizeof(int));
	closetab = calloc(numclose*totlen, sizeof(int));
	int index;
#pragma omp parallel for private (i,j,k,index), schedule(dynamic)
	for (i = 0; i<mlength; i++)
	{
		for (j = 0; j<zlength; j++)
		{
			for (k = 0; k<clength; k++)
			{
				index = index3D(zlength, clength, i, j, k);
				closemind[index] = config.msig - i;
				closezind[index] = config.zsig - j;
				closecind[index] = config.csig - k;
				//printf("m z c: %d %d %d %d\n",index,config.msig-i,config.zsig-j,config.csig-k);
			}
		}
	}


	printf("Number of Blurs: %d\n", numclose);

	int indm, indc, indz;
	double  point2;
#pragma omp parallel for private (indz,indc,indm,index,point2), schedule(dynamic)
	for (int i = 0; i<size[0]; i++)
	{
		for (int j = 0; j<size[1]; j++)
		{
			for (int k = 0; k<size[2]; k++)
			{
				for (int l = 0; l<size[3]; l++)
				{
					double ccspt = ccstab[index3D(size[1], size[2], i, j, k)];
					ccspt = ccspt + closecind[l] * config.ccsbins;
					indz = k + closezind[l];
					int badflag = 0;
					if (indz<0 || indz >= numz) { badflag = 1; }
					else
					{
						double point = calcMz(masstab[index2D(size[2], i, k)] + closemind[l] * config.molig, ztab[k] + closezind[l], config.adductmass);
						if (point<mzranges[0] || point>mzranges[1]) { badflag = 1; }
						else
						{
							indm = nearfast(mzext, point, size[0]);
							if (config.twaveflag == 0) { point2 = calcDt(point, ztab[k] + closezind[l], ccspt, ccsconst, config.hmass, config.to); }
							else if (config.twaveflag == 1) { point2 = calcDtTwaveLog(point, ztab[k] + closezind[l], ccspt, config.tcal1, config.tcal2, config.hmass, config.edc); }
							else if (config.twaveflag == 2) { point2 = calcDtTwaveLinear(point, ztab[k] + closezind[l], ccspt, config.tcal1, config.tcal2, config.hmass, config.edc); }
							else if (config.twaveflag == 3) { point2 = calcDtTwavePower(point, ztab[k] + closezind[l], ccspt, config.tcal1, config.tcal2, config.hmass, config.edc); }
							else { TwaveError(config.twaveflag); }
							if ((point2<mzranges[2] || point2>mzranges[3]) && config.csig >= 0) { badflag = 1; }
							else
							{
								indc = nearfast(dtext, point2, size[1]);
							}
						}
					}
					if (badflag == 1) { index = -1; }
					else { index = index3D(size[1], size[2], indm, indc, indz); }
					int index2 = index4D(size, i, j, k, l);
					closetab[index2] = index;

				}
			}
		}
	}
	printf("Finished Blur\n");
	//Setting Up the Iteration
	blur = calloc(totlen, sizeof(double));
	newblur = calloc(totlen, sizeof(double));
	//memset(barr,1,totlen);
#pragma omp parallel for private (i,j,k), schedule(dynamic)
	for (i = 0; i<size[0]; i++)
	{
		for (j = 0; j<size[1]; j++)
		{
			for (k = 0; k<size[2]; k++)
			{
				blur[index3D(size[1], size[2], i, j, k)] = dataInt[index2D(size[1], i, j)];
			}
		}
	}

	if (config.intthresh > 0) { KillB_IM(dataInt, barr, size, config.intthresh); }

	denom = calloc(lines, sizeof(double));
	deltas = calloc(lines, sizeof(double));
	printf("Iterating: \n");

	//Iterating
	for (m = 0; m<config.numit; m++) {

		blur_it_IM(size, blur, newblur, closetab, barr, config.csig);

		sumdeltas(size, deltas, newblur);
		fftconvolve2D(denom, deltas, peakshape, size);
		bayes(size, denom, blur, newblur, dataInt, barr);
		if (config.numit<10) { printf("Iteration: %d\n", m); }
		else { if (m % (config.numit / 10) == 0) { printf("Iteration: %d\n", m); } }
		//printf("Iteration: %d\n",m);
	}

	//Writing outputs

	//Write the Fit config.to the Experimental Data
	sumdeltas(size, deltas, blur);
	fftconvolve2D(denom, deltas, peakshape, size);

	double denommax = getmax(lines, denom);
	double datamax = getmax(lines, dataInt);
	normalize(lines, dataInt, datamax);
	normalize(lines, denom, denommax);
	ApplyCutoff1D(denom, 0, lines);
	//printf("maxes: %f %f\n",denommax,datamax);
	char *suffixfit = "fitdat";
	write1D(config.outfile, suffixfit, denom, lines);
	double error = errfun(lines, dataInt, denom);

	//Convolve the deltas with the peak shape and put that inconfig.to newblur
	printf("Convolving...");
	convolve3D(newblur, blur, peakshape, size);
	printf("Done\n");

	//Get Max and Min mass and ccs values
	double ranges[4];
	ranges[0] = config.masslb; ranges[1] = config.massub; ranges[2] = config.ccslb; ranges[3] = config.ccsub;
	if (config.fixedmassaxis == 0) { getranges(size, newblur, masstab, ccstab, ranges, barr); }
	printf("Mass Range: %f to %f  Mass Bins: %f Da\n", ranges[0], ranges[1], config.massbins);
	printf("CCS Range: %f to %f CCS Bins: %f\n", ranges[2], ranges[3], config.ccsbins);

	int maaxle = 1 + (int)((ranges[1] - ranges[0]) / config.massbins);
	int ccaxle = 1 + (int)((ranges[3] - ranges[2]) / config.ccsbins);
	int newsize[3];
	newsize[0] = maaxle;
	newsize[1] = ccaxle;
	newsize[2] = numz;
	printf("Dimensions of Output: %d m by %d ccs by %d z: %d total\n", newsize[0], newsize[1], newsize[2], newsize[0] * newsize[1] * newsize[2]);
	double *massaxis = NULL, *massvals = NULL, *ccsaxis = NULL, *ccsvals = NULL, *newgrid = NULL;
	if (newsize[0] * newsize[1] * newsize[2] * sizeof(double)>2E6) { printf("Warning: May exceed system memory capacity"); }

	massaxis = calloc(maaxle, sizeof(double));
	massvals = calloc(maaxle, sizeof(double));
	ccsaxis = calloc(ccaxle, sizeof(double));
	ccsvals = calloc(ccaxle, sizeof(double));
	newgrid = calloc(newsize[0] * newsize[1] * newsize[2], sizeof(double));
	memset(newgrid, 0, newsize[0] * newsize[1] * newsize[2]);
	makeaxis(massaxis, maaxle, ranges[0], config.massbins);
	makeaxis(ccsaxis, ccaxle, ranges[2], config.ccsbins);

	//Transform the Grid config.to Mass CCS space
	if (config.poolflag == 0) {
		//#pragma omp parallel for private (i,j,k,tempmass,tempccs,indm,indc), schedule(dynamic)
		for (i = 0; i<size[0]; i++)
		{
			for (j = 0; j<size[1]; j++)
			{
				for (k = 0; k<size[2]; k++)
				{
					tempmass = masstab[index2D(size[2], i, k)];
					tempccs = ccstab[index3D(size[1], size[2], i, j, k)];
					if (tempmass>massaxis[0] && tempmass<massaxis[newsize[0] - 1] && tempccs>ccsaxis[0] && tempccs<ccsaxis[newsize[1] - 1]) {
						indm = nearfast(massaxis, tempmass, newsize[0]);
						indc = nearfast(ccsaxis, tempccs, newsize[1]);
						//if(config.rawflag==0){
						//newgrid[index3D(newsize[1],newsize[2],indm,indc,k)]+=newblur[index3D(size[1],size[2],i,j,k)];}
						//else{newgrid[index3D(newsize[1],newsize[2],indm,indc,k)]+=blur[index3D(size[1],size[2],i,j,k)];}

						double val;
						if (config.rawflag == 0) {
							val = newblur[index3D(size[1], size[2], i, j, k)];
						}
						else { val = blur[index3D(size[1], size[2], i, j, k)]; }

						int indm2;
						int indc2;
						if (massaxis[indm] < tempmass) { indm2 = indm + 1; }
						else { indm2 = indm - 1; }
						if (ccsaxis[indc] < tempccs) { indc2 = indc + 1; }
						else { indc2 = indc - 1; }
						if (indm2 >= 0 && indm2 < maaxle&&indc2 >= 0 && indc2 < ccaxle) {
							double interposm = LinearInterpolatePosition(massaxis[indm], massaxis[indm2], tempmass);
							double interposc = LinearInterpolatePosition(ccsaxis[indc], ccsaxis[indc2], tempccs);
							double val1 = (1 - interposm)*(1 - interposc)*val;
							double val2 = (1 - interposm)*(interposc)*val;
							double val3 = (interposm)*(1 - interposc)*val;
							double val4 = (interposm)*(interposc)*val;
							newgrid[index3D(newsize[1], newsize[2], indm, indc, k)] += val1;
							newgrid[index3D(newsize[1], newsize[2], indm, indc2, k)] += val2;
							newgrid[index3D(newsize[1], newsize[2], indm2, indc, k)] += val3;
							newgrid[index3D(newsize[1], newsize[2], indm2, indc2, k)] += val4;
						}
					}
				}
			}
		}
		printf("Tranformed m/z and dt to mass and ccs by Integration\n");
	}
	else {
		double tempmz, tempdt, endval;
#pragma omp parallel for private (i,j,k,tempmass,tempmz,endval,tempccs,indm,indc), schedule(dynamic)
		for (i = 0; i<newsize[0]; i++)
		{
			for (j = 0; j<newsize[1]; j++)
			{
				for (k = 0; k<newsize[2]; k++)
				{
					tempmass = massaxis[i];
					tempmz = calcMz(tempmass, ztab[k], config.adductmass);
					tempccs = ccsaxis[j];

					if (config.twaveflag == 0) { tempdt = calcDt(tempmass, ztab[k], tempccs, ccsconst, config.hmass, config.to); }
					else if (config.twaveflag == 1) { tempdt = calcDtTwaveLog(tempmass, ztab[k], tempccs, config.tcal1, config.tcal2, config.hmass, config.edc); }
					else if (config.twaveflag == 2) { tempdt = calcDtTwaveLinear(tempmass, ztab[k], tempccs, config.tcal1, config.tcal2, config.hmass, config.edc); }
					else if (config.twaveflag == 3) { tempdt = calcDtTwavePower(tempmass, ztab[k], tempccs, config.tcal1, config.tcal2, config.hmass, config.edc); }
					else { TwaveError(config.twaveflag); }

					if (tempmz>mzext[1] && tempmz<mzext[size[0] - 2] && tempdt>dtext[1] && tempdt<dtext[size[1] - 2])
					{
						indm = nearfast(mzext, tempmz, size[0]);
						indc = nearfast(dtext, tempdt, size[1]);
						endval = bicubicinterpolation(size, indm, indc, k, mzext, dtext, tempmz, tempdt, config.rawflag, newblur, blur);//bicubic interpolation
																																 //endval=bilinearinterpolation(size,indm,indc,k,mzext,dtext,tempmz,tempdt,config.rawflag,newblur,blur);//bilinear interpolation
																																 //endval=newblur[index3D(size[1],size[2],indm,indc,k)]; //(Nearest Neighbor)
					}
					else { endval = 0; }
					newgrid[index3D(newsize[1], newsize[2], i, j, k)] = endval;
				}
			}
		}
		printf("Transformed m/z and dt to mass and ccs by Interpolation\n");
	}
	//Sum all config.to mass axis
	sum1D(newsize, massvals, newgrid, 0);
	sum1D(newsize, ccsvals, newgrid, 1);
	char *suffixmass = "mass";
	write2D(config.outfile, suffixmass, massaxis, massvals, maaxle);
	char *suffixccs = "ccs";
	write2D(config.outfile, suffixccs, ccsaxis, ccsvals, ccaxle);

	//Sum across charge and CCS
	double *massccsgrid = NULL, *masszgrid = NULL, *ccszgrid = NULL;//,*ccsaxis=NULL,*ccsvals=NULL,*newgrid=NULL;
	massccsgrid = calloc(maaxle*ccaxle, sizeof(double));
	masszgrid = calloc(maaxle*numz, sizeof(double));
	ccszgrid = calloc(ccaxle*numz, sizeof(double));
	sum2D(newsize, massccsgrid, newgrid, 2);
	sum2D(newsize, masszgrid, newgrid, 1);
	sum2D(newsize, ccszgrid, newgrid, 0);
	char *suffix3 = "massccs";
	write1D(config.outfile, suffix3, massccsgrid, maaxle*ccaxle);
	char *suffix4 = "massgrid";
	write1D(config.outfile, suffix4, masszgrid, maaxle*numz);
	char *suffix6 = "ccsz";
	write1D(config.outfile, suffix6, ccszgrid, ccaxle*numz);

	if (config.mflag == 1) {
		MFileInt(tmlen, massaxis, massvals, testmasses, testmassvals, maaxle);
		MFileCCS(tmlen, massaxis, ccsaxis, massccsgrid, testmasses, testmassCCSavg, testmassCCSstd, maaxle, ccaxle);
		MFileZ(tmlen, massaxis, ztab, masszgrid, testmasses, testmassZavg, testmassZstd, maaxle, numz);
		char *suffixmfileint = "mfileresults";
		writemfileres(config.outfile, suffixmfileint, testmasses, testmassvals, testmassCCSavg, testmassCCSstd, testmassZavg, testmassZstd, tmlen);
	}

	//Write out specific charge slices
	char *suffix5 = "zout";
	if (config.zout>0) {
		writezslice(newsize, config.outfile, suffix5, massaxis, ccsaxis, ztab, newgrid, nearint(ztab, config.zout, numz));
	}
	if (config.zout < 0) {
		for (int k = 0; k < numz; k++)
		{
			writezslice(newsize, config.outfile, suffix5, massaxis, ccsaxis, ztab, newgrid, k);
		}
	}
	char *suffixgrid = "mzgrid";
	writemzgrid(config.outfile, suffixgrid, mzext, dtext, ztab, newblur, size);

	char outstring[500];
	char *suffixerr = "error";
	FILE *out_ptrIM = NULL;
	sprintf(outstring, "%s_%s.txt", config.outfile, suffixerr);
	out_ptrIM = fopen(outstring, "w");
	fprintf(out_ptrIM, "error %f\n", error);
	fclose(out_ptrIM);
	printf("File written to: %s\n", outstring);


	printf("\nR Squared: %f\n", error);
	free(massaxis);
	free(massvals);
	free(ccsaxis);
	free(ccsvals);
	free(newgrid);
	free(barr);
	free(mzdat);
	free(dtdat);
	free(dataInt);
	free(ccstab);
	free(ztab);
	free(masstab);
	free(blur);
	free(newblur);
	free(peakshape);
	free(denom);
	free(deltas);
	free(mzext);
	free(dtext);
	free(closemind);
	free(closezind);
	free(closecind);
	free(closetab);
	free(testmasses);
	free(testmassvals);
	free(testmassCCSavg);
	free(testmassCCSstd);
	free(testmassZavg);
	free(testmassZstd);
	free(massccsgrid);
	free(masszgrid);
	free(ccszgrid);
	endtime = time(NULL);
	printf("\nDone in %ds!\n", (int)difftime(endtime, starttime));
	return 0;
}