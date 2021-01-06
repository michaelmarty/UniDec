/*
* UD_score.h
*
*  Created on : 20 August 2019
* Author : Michael.Marty
*/

//
// 
// Copyright 2019 University of Arizona
//
//

#include "UD_analysis.h"
//#include "UniDec.h"
//#include "UD_peak_width.h"

void WritePeaks(const Config config, const Decon* decon);


void get_fwhms(Config config, const int plen, const int mlen, const float* massaxis, const float* masssum, const float* peakx, float* fwhmlow, float* fwhmhigh, float* badfwhm)
{
	for (int i = 0; i < plen; i++)
	{
		float peak = peakx[i];
		int index = nearfast(massaxis, peak, mlen);
		float max = masssum[index];
		
		float halfmax = max / 2.0;
		int lindex = index;
		int hindex = index;
		
		bool lfound = 0;
		while (lfound == 0)
		{
			lindex -= 1;
			if (lindex < 0)
			{
				lfound = 1;
			}
			else{
				float lval = masssum[lindex];
				if (lval <= halfmax)
				{
					lfound = 1;
				}
			}
		}

		bool hfound = 0;
		while (hfound == 0)
		{
			hindex += 1;
			if (hindex >= mlen)
			{
				hfound = 1;
			}
			else{
				float hval = masssum[hindex];
				if (hval <= halfmax)
				{
					hfound = 1;
				}
			}
			
		}

		float mlow = massaxis[lindex];
		float mhigh = massaxis[hindex];

		//Catch peaks that are too assymetric. Fix and flag them.
		float highdiff = mhigh - peak;
		float lowdiff = peak - mlow;
		float fwhm = mhigh - mlow;

		float threshold = 0.75;
		float mult = threshold / (1 - threshold);
		if (fwhm != 0) {
			float rh = highdiff / fwhm;
			float rl = lowdiff / fwhm;
			if (rh > threshold)
			{
				mhigh = peak + lowdiff * mult;
				badfwhm[i] = 1;
			}
			else if (rl > threshold)
			{
				mlow = peak - highdiff * mult;
				badfwhm[i] = 1;
			}
		}
		else
		{
			printf("Warning: FWHM was 0.\n");
			//mlow = peak - config.massbins * 3;
			//mhigh = peak + config.massbins * 3;
		}

		fwhmlow[i] = mlow;
		fwhmhigh[i] = mhigh;
		//printf("%f %f %f\n", massaxis[index], mlow, mhigh);
	}
}

// To lower duplication, consider calling this in the loop of get_fwhms()?
float single_fwhm(Config config, const int mlen, const float* massaxis, const float* masssum, const float peak, int index, float max)
{
	float halfmax = max / 2.0;
	int lindex = index;
	int hindex = index;

	bool lfound = 0;
	while (lfound == 0)
	{
		lindex -= 1;
		if (lindex < 0)
		{
			lfound = 1;
		}
		else {
			float lval = masssum[lindex];
			if (lval <= halfmax)
			{
				lfound = 1;
			}
		}
	}

	bool hfound = 0;
	while (hfound == 0)
	{
		hindex += 1;
		if (hindex >= mlen)
		{
			hfound = 1;
		}
		else {
			float hval = masssum[hindex];
			if (hval <= halfmax)
			{
				hfound = 1;
			}
		}

	}

	float mlow = massaxis[lindex];
	float mhigh = massaxis[hindex];

	//Catch peaks that are too assymetric. Fix and flag them.
	float highdiff = mhigh - peak;
	float lowdiff = peak - mlow;
	float fwhm = mhigh - mlow;

	float threshold = 0.75;
	float mult = threshold / (1 - threshold);
	if (fwhm != 0) {
		float rh = highdiff / fwhm;
		float rl = lowdiff / fwhm;
		if (rh > threshold)
		{
			mhigh = peak + lowdiff * mult;
		}
		else if (rl > threshold)
		{
			mlow = peak - highdiff * mult;
		}
	}
	else
	{
		printf("Warning: FWHM was 0.\n");
	}
	fwhm = mhigh - mlow;
	// printf("FWHM is %f\n", fwhm);

	return fwhm;
}

float uscore(Config config, const float* dataMZ, const float* dataInt, const float* mzgrid, const int* nztab, 
	const float mlow, const float mhigh, const float peak)
{
	float power = 2;
	//float* errors = NULL;
	//float* weights = NULL;
	//errors = calloc(config.numz, sizeof(float));
	//weights = calloc(config.numz, sizeof(float));

	float numerator = 0;
	float denominator = 0;

	float datamin = dataMZ[0];
	float datamax = dataMZ[config.lengthmz - 1];
	
	for (int i = 0; i < config.numz; i++)
	{
		float z = (float)nztab[i];
		float lmz = (mlow + z * config.adductmass) / z;
		float hmz = (mhigh + z * config.adductmass) / z;

		if ( hmz > datamin && lmz < datamax)
		{
			int lindex = nearfast(dataMZ, lmz, config.lengthmz);
			int hindex = nearfast(dataMZ, hmz, config.lengthmz);

			if (dataMZ[lindex] < lmz) { lindex += 1; }
			if (dataMZ[hindex] > hmz) { hindex -= 1; }

			//printf("%f %f\n", dataMZ[lindex], dataMZ[hindex]);

			float sumdata = 0;
			float sumerrors = 0;
			float sumdecon = 0;

			for (int j = lindex; j <= hindex; j++)
			{
				float data = dataInt[j];
				float decon = mzgrid[index2D(config.numz, j, i)];
				if (config.orbimode == 1) { decon *= z; }
				sumerrors += fabs(data - decon);
				sumdata += data;
				sumdecon += decon;
			}

			float per = 0;
			if (sumdata > 0)
			{
				per = 1 - (sumerrors / sumdata);
			}
			//printf("%f %f\n", mzgrid[index2D(config.numz, lindex, i)], mzgrid[index2D(config.numz, hindex, i)]);
			//printf("%f %f %f %f %f\n", z, sumerrors, per, sumdecon, sumdata);
			sumdecon = pow(sumdecon, power);
			numerator += sumdecon * per;
			denominator += sumdecon;
		}
	}

	if (denominator != 0)
	{
		return numerator / denominator;
	}
	else { return 0; }

}


float mscore(Config config, const int mlen, const float* massaxis, const float* masssum, const float* massgrid, const float mlow, const float mhigh, const float peak)
{
	float power = 2;
	float mscore = 0;

	float numerator = 0;
	float denominator = 0;

	float datamin = massaxis[0];
	float datamax = massaxis[mlen - 1];

	if (mhigh > datamin && mlow < datamax)
	{
		int lindex = nearfast(massaxis, mlow, mlen);
		int hindex = nearfast(massaxis, mhigh, mlen);

		if (massaxis[lindex] < mlow) { lindex += 1; }
		if (massaxis[hindex] > mhigh) { hindex -= 1; }
		//printf("%f %f\n", massaxis[lindex], massaxis[hindex]);

		float mmax = 0;
		float msum = 0;
		for (int j = lindex; j <= hindex; j++)
		{
			float data = masssum[j];
			if (data > mmax) { mmax = data; }
			msum += data;
		}


		for (int i = 0; i < config.numz; i++)
		{
			float sumdata = 0;
			float sumerrors = 0;
			float sumdecon = 0;
			for (int j = lindex; j <= hindex; j++)
			{
				float decon = massgrid[index2D(config.numz, j, i)];
				sumdecon += decon;
			}

			//float sumdecon2 = 0;
			if (sumdecon > 0)
			{
				for (int j = lindex; j <= hindex; j++)
				{
					float data = masssum[j];
					float decon = massgrid[index2D(config.numz, j, i)] / sumdecon * msum;
					sumerrors += fabs(data - decon);
					sumdata += data;
					//sumdecon2 += decon;
				}

				float per = 0;
				if (sumdata > 0)
				{
					per = 1 - (sumerrors / sumdata);
				}

				//printf("%f %f\n", mzgrid[index2D(config.numz, lindex, i)], mzgrid[index2D(config.numz, hindex, i)]);
				//printf("%f %f %f %f %f\n", z, sumerrors, per, sumdecon, sumdata);
				sumdecon = pow(sumdecon, power);
				numerator += sumdecon * per;
				denominator += sumdecon;
			}
		}
	}
	if (denominator != 0)
	{
		return numerator / denominator;
	}
	else { return 0; }

}


float csscore(Config config, const int mlen, const float* massaxis, const float* masssum, const float* massgrid, const float mlow, const float mhigh, const float peak)
{
	float csscore = 0;
	
	float* zvals = NULL;
	zvals = calloc(config.numz, sizeof(float));

	float datamin = massaxis[0];
	float datamax = massaxis[mlen - 1];

	if (mhigh > datamin && mlow < datamax)
	{
		int lindex = nearfast(massaxis, mlow, mlen);
		int hindex = nearfast(massaxis, mhigh, mlen);

		if (massaxis[lindex] < mlow) { lindex += 1; }
		if (massaxis[hindex] > mhigh) { hindex -= 1; }

		float zmax = 0;
		float zsum = 0;
		int zmaxindex = -1;
		for (int i = 0; i < config.numz; i++)
		{
			float zval = 0;
			for (int j = lindex; j <= hindex; j++)
			{
				zval += massgrid[index2D(config.numz, j, i)];
			}
			zvals[i] = zval;
			zsum += zval;
			if (zval > zmax) { zmax = zval; zmaxindex = i; }
		}

		float badarea = 0;
		int index = zmaxindex;
		float lowval = zmax;
		while (index < config.numz - 1)
		{
			index += 1;
			float v = zvals[index];
			if (v < lowval)
			{
				lowval = v;
			}
			else {
				badarea += v - lowval;
			}
		}

		index = zmaxindex;
		lowval = zmax;
		while (index > 0)
		{
			index -= 1;
			float v = zvals[index];
			if (v < lowval)
			{
				lowval = v;
			}
			else {
				badarea += v - lowval;
			}
		}

		if (zsum != 0)
		{
			csscore = 1 - (badarea / zsum);
		}

	}
	return csscore;
}


float find_minimum(Config config, const int mlen, const float* massaxis, const float* masssum, const float lowpt, const float highpt)
{
	int lindex = nearfast(massaxis, lowpt, mlen);
	int hindex = nearfast(massaxis, highpt, mlen);

	float minval = masssum[hindex];
	for (int i = lindex; i < hindex; i++)
	{
		float val = masssum[i];
		if (val < minval) { minval = val; }
	}
	return minval;
}

float score_minimum(float height, float min)
{
	float hh = height / 2;
	if (min > hh) {
		return 1 - ((min - hh) / (height - hh));
	}
	else { return 1; }
}

float fscore(Config config, const int plen, const int mlen, const float* massaxis, const float* masssum, const float *peakx, const float height, 
	const float mlow, const float mhigh, const float peak, const int badfwhm)
{
	float fscore = 1;

	//Score down peaks that are highly assymetic.
	if (badfwhm == 1) {
		//printf("Badfwhm\n");
		float highdiff = mhigh - peak;
		float lowdiff = peak - mlow;
		float fwhm = mhigh - mlow;
		if (lowdiff > highdiff) {
			float min = find_minimum(config, mlen, massaxis, masssum, mlow - config.massbins, peak);
			float fsc = score_minimum(height, min);
			//printf("Fscore2 %f %f %f\n", fsc, min, height);
			fscore *= fsc;
		}
		else {
			float min = find_minimum(config, mlen, massaxis, masssum, peak, mhigh + config.massbins);
			float fsc = score_minimum(height, min);
			//printf("Fscore3 %f %f %f\n", fsc, min, height);
			fscore *= fsc;
		}
	}

	//Test if the peaks are actually separated with at least a half max height
	for (int i = 0; i < plen; i++)
	{
		float peak2 = peakx[i];
		if (peak != peak2)
		{
			if (peak2 < peak && peak2 > mlow)
			{
				float min = find_minimum(config, mlen, massaxis, masssum, peak2, peak);
				float fsc = score_minimum(height, min);
				//printf("Fscore4 %f %f %f %f\n", fsc, min, height, peak2);
				fscore *= fsc;
			}
			if (peak2 > peak && peak2 < mhigh)
			{
				float min = find_minimum(config, mlen, massaxis, masssum, peak, peak2);
				float fsc = score_minimum(height, min);
				//printf("Fscore5 %f %f %f %f\n", fsc, min, height, peak2);
				fscore *= fsc;
			}
		}
	}

	return fscore;
}


float score_from_peaks(const int plen, const float *peakx, const float *peaky, float *dscores, const Config config, Decon *decon, const Input inp, const float threshold) {

	float xfwhm = 2;
	float* fwhmlow = NULL;
	float* fwhmhigh = NULL;
	float* badfwhm = NULL;
	fwhmlow = calloc(plen, sizeof(float));
	fwhmhigh = calloc(plen, sizeof(float));
	badfwhm = calloc(plen, sizeof(float));

	get_fwhms(config, plen, decon->mlen, decon->massaxis, decon->massaxisval, peakx, fwhmlow, fwhmhigh, badfwhm);

	float numerator = 0;
	float denominator = 0;
	float uniscore = 0;

	for (int i = 0; i < plen; i++)
	{
		float m = peakx[i];
		float ival = peaky[i];
		float l = m - (m - fwhmlow[i]) * xfwhm;
		float h = m + (fwhmhigh[i] - m) * xfwhm;
		int index = nearfast(decon->massaxis, m, decon->mlen);
		float height = decon->massaxisval[index];

		float usc = uscore(config, inp.dataMZ, inp.dataInt, decon->newblur, inp.nztab, l, h, m);
		float msc = mscore(config, decon->mlen, decon->massaxis, decon->massaxisval, decon->massgrid, l, h, m);
		float cssc = csscore(config, decon->mlen, decon->massaxis, decon->massaxisval, decon->massgrid, l, h, m);
		float fsc = fscore(config, plen, decon->mlen, decon->massaxis, decon->massaxisval, peakx, height, fwhmlow[i], fwhmhigh[i], m, badfwhm[i]);

		float dsc = usc * msc * cssc * fsc;
		dscores[i] = dsc;
		if (dsc > threshold)
		{
			//printf("Peak: Mass:%f Int:%f M:%f U:%f CS:%f F:%f D: %f \n", peakx[i], peaky[i], msc, usc, cssc, fsc, dsc);
			numerator += ival * ival * dsc;
			denominator += ival * ival;
		}
	}

	//printf("R Squared: %f\n", decon->rsquared);

	if (denominator != 0) { uniscore = decon->rsquared * numerator / denominator; }
	return uniscore;
}

float score(Config config, Decon *decon, Input inp, const float threshold)
{
	//printf("Starting Score %f %f\n", config.peakwin, config.peakthresh);
	
	decon->peakx = calloc(decon->mlen, sizeof(float));
	decon->peaky = calloc(decon->mlen, sizeof(float));
	
	int plen = peak_detect(decon->massaxis, decon->massaxisval, decon->mlen, config.peakwin, config.peakthresh, decon->peakx, decon->peaky);
	decon->plen = plen;

	decon->peakx = realloc(decon->peakx, plen * sizeof(float));
	decon->peaky = realloc(decon->peaky, plen * sizeof(float));
	decon->dscores = calloc(plen, sizeof(float));

	peak_norm(decon->peaky, plen, config.peaknorm);

	float uniscore = score_from_peaks(plen, decon->peakx, decon->peaky, decon->dscores, config, decon, inp, threshold);
	printf("Average Peaks Score (UniScore): %f\n", uniscore);

	return uniscore;
}

int ReadDecon(Config* config, const Input inp, Decon* decon) 
{
	char outdat[1024];
	char strval[1024];

	//Import Rsquared
	decon->rsquared = float_attr(config->file_id, config->dataset, "rsquared", decon->rsquared);

	//Mass Data
	strjoin(config->dataset, "/mass_data", outdat);
	//printf("\tReading: %s\n", outdat);
	decon->mlen = mh5getfilelength(config->file_id, outdat);
	decon->massaxis = calloc(decon->mlen, sizeof(float));
	decon->massaxisval = calloc(decon->mlen, sizeof(float));
	mh5readfile2d(config->file_id, outdat, decon->mlen, decon->massaxis, decon->massaxisval);

	//MZ Grid
	strjoin(config->dataset, "/mz_grid", outdat);
	int status = check_group_noexit(config->file_id, outdat);
	if (status == 0) { return 0; }
	//printf("\tReading: %s\n", outdat);
	decon->newblur = calloc(config->lengthmz * config->numz, sizeof(float));
	mh5readfile1d(config->file_id, outdat, decon->newblur);

	//Mass Grid
	strjoin(config->dataset, "/mass_grid", outdat);
	status = check_group_noexit(config->file_id, outdat);
	if (status == 0) { return 0; }
	//printf("\tReading: %s\n", outdat);
	decon->massgrid = calloc(decon->mlen * config->numz, sizeof(float));
	mh5readfile1d(config->file_id, outdat, decon->massgrid);
	return 1;
}

void get_scan_scores(int argc, char* argv[], Config config)
{
	config.file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
	int num = 0;
	num = int_attr(config.file_id, "/ms_dataset", "num", num);

	for (int i = 0; i < num; i++) {
		config.metamode = i;
		
		Decon decon = SetupDecon();
		Input inp = SetupInputs();
		ReadInputs(argc, argv, &config, &inp);
		int status = ReadDecon(&config, inp, &decon);

		if(status ==1){
			score(config, &decon, inp, 0);
			WritePeaks(config, &decon);
		}
		else
		{
			printf("Missing deconvolution outputs. Turn off Fast Profile/Fast Centroid and try deconvolving again.");
		}

		FreeDecon(decon);
		FreeInputs(inp);
	}
	H5Fclose(config.file_id);
}


