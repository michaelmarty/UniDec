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



void get_fwhms(Config config, const int mlen, const double* massaxis, const double* masssum, const double* peakx, double* fwhmlow, double* fwhmhigh, double* badfwhm)
{
	for (int i = 0; i < config.plen; i++)
	{
		double peak = peakx[i];
		int index = nearfast(massaxis, peak, mlen);
		double max = masssum[index];
		
		double halfmax = max / 2.0;
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
				double lval = masssum[lindex];
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
				double hval = masssum[hindex];
				if (hval <= halfmax)
				{
					hfound = 1;
				}
			}
			
		}

		double mlow = massaxis[lindex];
		double mhigh = massaxis[hindex];

		//Catch peaks that are too assymetric. Fix and flag them.
		double highdiff = mhigh - peak;
		double lowdiff = peak - mlow;
		double fwhm = mhigh - mlow;

		double threshold = 0.75;
		double mult = threshold / (1 - threshold);
		if (fwhm != 0) {
			double rh = highdiff / fwhm;
			double rl = lowdiff / fwhm;
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
double single_fwhm(Config config, const int mlen, const double* massaxis, const double* masssum, const double peak, int index, double max)
{
	double halfmax = max / 2.0;
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
			double lval = masssum[lindex];
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
			double hval = masssum[hindex];
			if (hval <= halfmax)
			{
				hfound = 1;
			}
		}

	}

	double mlow = massaxis[lindex];
	double mhigh = massaxis[hindex];

	//Catch peaks that are too assymetric. Fix and flag them.
	double highdiff = mhigh - peak;
	double lowdiff = peak - mlow;
	double fwhm = mhigh - mlow;

	double threshold = 0.75;
	double mult = threshold / (1 - threshold);
	if (fwhm != 0) {
		double rh = highdiff / fwhm;
		double rl = lowdiff / fwhm;
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

double uscore(Config config, const double* dataMZ, const double* dataInt, const double* mzgrid, const int* nztab, 
	const double mlow, const double mhigh, const double peak)
{
	double power = 2;
	//double* errors = NULL;
	//double* weights = NULL;
	//errors = calloc(config.numz, sizeof(double));
	//weights = calloc(config.numz, sizeof(double));

	double numerator = 0;
	double denominator = 0;

	double datamin = dataMZ[0];
	double datamax = dataMZ[config.lengthmz - 1];
	
	for (int i = 0; i < config.numz; i++)
	{
		double z = (double)nztab[i];
		double lmz = (mlow + z * config.adductmass) / z;
		double hmz = (mhigh + z * config.adductmass) / z;

		if ( hmz > datamin && lmz < datamax)
		{
			int lindex = nearfast(dataMZ, lmz, config.lengthmz);
			int hindex = nearfast(dataMZ, hmz, config.lengthmz);

			if (dataMZ[lindex] < lmz) { lindex += 1; }
			if (dataMZ[hindex] > hmz) { hindex -= 1; }

			//printf("%f %f\n", dataMZ[lindex], dataMZ[hindex]);

			double sumdata = 0;
			double sumerrors = 0;
			double sumdecon = 0;

			for (int j = lindex; j <= hindex; j++)
			{
				double data = dataInt[j];
				double decon = mzgrid[index2D(config.numz, j, i)];
				sumerrors += fabs(data - decon);
				sumdata += data;
				sumdecon += decon;
			}

			double per = 0;
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


double mscore(Config config, const int mlen, const double* massaxis, const double* masssum, const double* massgrid, const double mlow, const double mhigh, const double peak)
{
	double power = 2;
	double mscore = 0;

	double numerator = 0;
	double denominator = 0;

	double datamin = massaxis[0];
	double datamax = massaxis[mlen - 1];

	if (mhigh > datamin && mlow < datamax)
	{
		int lindex = nearfast(massaxis, mlow, mlen);
		int hindex = nearfast(massaxis, mhigh, mlen);

		if (massaxis[lindex] < mlow) { lindex += 1; }
		if (massaxis[hindex] > mhigh) { hindex -= 1; }
		//printf("%f %f\n", massaxis[lindex], massaxis[hindex]);

		double mmax = 0;
		double msum = 0;
		for (int j = lindex; j <= hindex; j++)
		{
			double data = masssum[j];
			if (data > mmax) { mmax = data; }
			msum += data;
		}


		for (int i = 0; i < config.numz; i++)
		{
			double sumdata = 0;
			double sumerrors = 0;
			double sumdecon = 0;
			for (int j = lindex; j <= hindex; j++)
			{
				double decon = massgrid[index2D(config.numz, j, i)];
				sumdecon += decon;
			}

			//double sumdecon2 = 0;
			if (sumdecon > 0)
			{
				for (int j = lindex; j <= hindex; j++)
				{
					double data = masssum[j];
					double decon = massgrid[index2D(config.numz, j, i)] / sumdecon * msum;
					sumerrors += fabs(data - decon);
					sumdata += data;
					//sumdecon2 += decon;
				}

				double per = 0;
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


double csscore(Config config, const int mlen, const double* massaxis, const double* masssum, const double* massgrid, const double mlow, const double mhigh, const double peak)
{
	double csscore = 0;
	
	double* zvals = NULL;
	zvals = calloc(config.numz, sizeof(double));

	double datamin = massaxis[0];
	double datamax = massaxis[mlen - 1];

	if (mhigh > datamin && mlow < datamax)
	{
		int lindex = nearfast(massaxis, mlow, mlen);
		int hindex = nearfast(massaxis, mhigh, mlen);

		if (massaxis[lindex] < mlow) { lindex += 1; }
		if (massaxis[hindex] > mhigh) { hindex -= 1; }

		double zmax = 0;
		double zsum = 0;
		int zmaxindex = -1;
		for (int i = 0; i < config.numz; i++)
		{
			double zval = 0;
			for (int j = lindex; j <= hindex; j++)
			{
				zval += massgrid[index2D(config.numz, j, i)];
			}
			zvals[i] = zval;
			zsum += zval;
			if (zval > zmax) { zmax = zval; zmaxindex = i; }
		}

		double badarea = 0;
		int index = zmaxindex;
		double lowval = zmax;
		while (index < config.numz - 1)
		{
			index += 1;
			double v = zvals[index];
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
			double v = zvals[index];
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


double find_minimum(Config config, const int mlen, const double* massaxis, const double* masssum, const double lowpt, const double highpt)
{
	int lindex = nearfast(massaxis, lowpt, mlen);
	int hindex = nearfast(massaxis, highpt, mlen);

	double minval = masssum[hindex];
	for (int i = lindex; i < hindex; i++)
	{
		double val = masssum[i];
		if (val < minval) { minval = val; }
	}
	return minval;
}

double score_minimum(double height, double min)
{
	double hh = height / 2;
	if (min > hh) {
		return 1 - ((min - hh) / (height - hh));
	}
	else { return 1; }
}

double fscore(Config config, const int mlen, const double* massaxis, const double* masssum, const double *peakx, const double height, 
	const double mlow, const double mhigh, const double peak, const int badfwhm)
{
	double fscore = 1;

	//Score down peaks that are highly assymetic.
	if (badfwhm == 1) {
		//printf("Badfwhm\n");
		double highdiff = mhigh - peak;
		double lowdiff = peak - mlow;
		double fwhm = mhigh - mlow;
		if (lowdiff > highdiff) {
			double min = find_minimum(config, mlen, massaxis, masssum, mlow - config.massbins, peak);
			double fsc = score_minimum(height, min);
			//printf("Fscore2 %f %f %f\n", fsc, min, height);
			fscore *= fsc;
		}
		else {
			double min = find_minimum(config, mlen, massaxis, masssum, peak, mhigh + config.massbins);
			double fsc = score_minimum(height, min);
			//printf("Fscore3 %f %f %f\n", fsc, min, height);
			fscore *= fsc;
		}
	}

	//Test if the peaks are actually separated with at least a half max height
	for (int i = 0; i < config.plen; i++)
	{
		double peak2 = peakx[i];
		if (peak != peak2)
		{
			if (peak2 < peak && peak2 > mlow)
			{
				double min = find_minimum(config, mlen, massaxis, masssum, peak2, peak);
				double fsc = score_minimum(height, min);
				//printf("Fscore4 %f %f %f %f\n", fsc, min, height, peak2);
				fscore *= fsc;
			}
			if (peak2 > peak && peak2 < mhigh)
			{
				double min = find_minimum(config, mlen, massaxis, masssum, peak, peak2);
				double fsc = score_minimum(height, min);
				//printf("Fscore5 %f %f %f %f\n", fsc, min, height, peak2);
				fscore *= fsc;
			}
		}
	}

	return fscore;
}


double score(Config config, Decon decon, const double * dataMZ, const double * dataInt, const double *mzgrid, const double *massaxis, const double* masssum, const double *massgrid, const int *nztab, const double threshold)
{
	//printf("Starting Score %f %f\n", config.peakwin, config.peakthresh);
	double xfwhm = 2;

	double* peakx = NULL;
	double* peaky = NULL;
	peakx = calloc(decon.mlen, sizeof(double));
	peaky = calloc(decon.mlen, sizeof(double));

	int plen = peak_detect(massaxis, masssum, decon.mlen, config.peakwin, config.peakthresh, peakx, peaky);
	config.plen = plen;

	peakx = realloc(peakx, plen * sizeof(double));
	peaky = realloc(peaky, plen * sizeof(double));

	peak_norm(peaky, plen, config.peaknorm);

	double* fwhmlow = NULL;
	double* fwhmhigh = NULL;
	double* badfwhm = NULL;
	fwhmlow = calloc(plen, sizeof(double));
	fwhmhigh = calloc(plen, sizeof(double));
	badfwhm = calloc(plen, sizeof(double));

	get_fwhms(config, decon.mlen, massaxis, masssum, peakx, fwhmlow, fwhmhigh, badfwhm);

	double numerator = 0;
	double denominator = 0;
	double uniscore = 0;

	for (int i = 0; i < plen; i++)
	{
		double m = peakx[i];
		double ival = peaky[i];
		double l = m - (m-fwhmlow[i])*xfwhm;
		double h = m + (fwhmhigh[i]-m)*xfwhm;
		int index = nearfast(massaxis, m, decon.mlen);
		double height = masssum[index];
		
		double usc = uscore(config, dataMZ, dataInt, mzgrid, nztab, l, h, m);
		double msc = mscore(config, decon.mlen, massaxis, masssum, massgrid, l, h, m);
		double cssc = csscore(config, decon.mlen, massaxis, masssum, massgrid, l, h, m);
		double fsc = fscore(config, decon.mlen, massaxis, masssum, peakx, height, fwhmlow[i], fwhmhigh[i], m, badfwhm[i]);

		double dsc = usc * msc * cssc * fsc;
		if (dsc > threshold)
		{
			//printf("Peak: Mass:%f Int:%f M:%f U:%f CS:%f F:%f D: %f \n", peakx[i], peaky[i], msc, usc, cssc, fsc, dsc);
			numerator += ival * ival * dsc;
			denominator += ival * ival;
		}
	}

	printf("R Squared: %f\n", decon.rsquared);

	if (denominator != 0) { uniscore = decon.rsquared * numerator / denominator; }
	
	printf("Average Peaks Score (UniScore): %f\n", uniscore);
	return uniscore;
}