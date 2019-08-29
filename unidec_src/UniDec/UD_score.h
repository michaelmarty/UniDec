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



void get_fwhms(Config config, const double* massaxis, const double* masssum, const double* peakx, double* fwhmlow, double* fwhmhigh)
{
	for (int i = 0; i < config.plen; i++)
	{
		int index = nearfast(massaxis, peakx[i], config.mlen);
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
			if (hindex >= config.mlen)
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

		fwhmlow[i] = massaxis[lindex];
		fwhmhigh[i] = massaxis[hindex];
		//printf("%f %f %f\n", massaxis[index], massaxis[lindex], massaxis[hindex]);
	}
}

double uscore(Config config, const double* dataMZ, const double* dataInt, const double* mzgrid, const int* nztab, 
	const double mlow, const double mhigh, const double peak)
{
	double power = 2;
	double uscore = 0;
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


double mscore(Config config, const double* massaxis, const double* masssum, const double* massgrid, const double mlow, const double mhigh, const double peak)
{
	double power = 2;
	double mscore = 0;

	double numerator = 0;
	double denominator = 0;

	double datamin = massaxis[0];
	double datamax = massaxis[config.mlen - 1];

	if (mhigh > datamin && mlow < datamax)
	{
		int lindex = nearfast(massaxis, mlow, config.mlen);
		int hindex = nearfast(massaxis, mhigh, config.mlen);

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


double zscore(Config config, const double* massaxis, const double* masssum, const double* massgrid, const double mlow, const double mhigh, const double peak)
{
	double zscore = 0;
	
	double* zvals = NULL;
	zvals = calloc(config.numz, sizeof(double));

	double datamin = massaxis[0];
	double datamax = massaxis[config.mlen - 1];

	if (mhigh > datamin && mlow < datamax)
	{
		int lindex = nearfast(massaxis, mlow, config.mlen);
		int hindex = nearfast(massaxis, mhigh, config.mlen);

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
			zscore = 1 - (badarea / zsum);
		}

	}
	return zscore;
}


double find_minimum(Config config, const double* massaxis, const double* masssum, const double lowpt, const double highpt)
{
	int lindex = nearfast(massaxis, lowpt, config.mlen);
	int hindex = nearfast(massaxis, highpt, config.mlen);

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

double fscore(Config config, const double* massaxis, const double* masssum, const double *peakx, const double *peaky, const double mlow, const double mhigh, const double peak)
{
	double fscore = 1;
	//First, test if the FWHM interval is highly asymetric 
	double threshold = 0.8;

	double highdiff = mhigh - peak;
	double lowdiff = peak - mlow;
	double fwhm = mhigh - mlow;
	
	if (fwhm != 0) {
		double rh = highdiff / fwhm;
		double rl = lowdiff / fwhm;
		if (rh > threshold)
		{
			fscore = 1 - ((rh - threshold) / (1 - threshold));
		}
		else if (rl > threshold)
		{
			fscore = 1 - ((rl - threshold) / (1 - threshold));
		}
	}
	else
	{
		printf("Warning: FWHM was 0.");
		return 0;
	}
	
	//Now test if the peaks are actually separated with at least a half max height
	for (int i = 0; i < config.plen; i++)
	{
		double peak2 = peakx[i];
		if (peak != peak2)
		{
			if (peak2 < peak && peak2 > mlow)
			{
				int index = nearfast(massaxis, peak, config.mlen);
				double height = masssum[index];
				double min = find_minimum(config, massaxis, masssum, peak2, peak);
				double fsc = score_minimum(height, min);
				fscore *= fsc;
			}
			if (peak2 > peak && peak2 < mhigh)
			{
				int index = nearfast(massaxis, peak, config.mlen);
				double height = masssum[index];
				double min = find_minimum(config, massaxis, masssum, peak, peak2);
				double fsc = score_minimum(height, min);
				fscore *= fsc;
			}
		}
	}

	return fscore;
}


double score(Config config, const double * dataMZ, const double * dataInt, const double *mzgrid, const double *massaxis, const double* masssum, const double *massgrid, const int *nztab, const double threshold)
{
	printf("Starting Score %f %f\n", config.peakwin, config.peakthresh);
	double xfwhm = 3;

	double* peakx = NULL;
	double* peaky = NULL;
	int mlen = config.mlen;
	peakx = calloc(mlen, sizeof(double));
	peaky = calloc(mlen, sizeof(double));

	int plen = peak_detect(massaxis, masssum, mlen, config.peakwin, config.peakthresh, peakx, peaky);
	config.plen = plen;

	peakx = realloc(peakx, plen * sizeof(double));
	peaky = realloc(peaky, plen * sizeof(double));

	peak_norm(peaky, plen, config.peaknorm);

	double* fwhmlow = NULL;
	double* fwhmhigh = NULL;
	fwhmlow = calloc(plen, sizeof(double));
	fwhmhigh = calloc(plen, sizeof(double));

	get_fwhms(config, massaxis, masssum, peakx, fwhmlow, fwhmhigh);

	double numerator = 0;
	double denominator = 0;
	double avgscore = 0;

	for (int i = 0; i < plen; i++)
	{
		double m = peakx[i];
		double ival = peaky[i];
		double l = m - (m-fwhmlow[i])*xfwhm;
		double h = m + (fwhmhigh[i]-m)*xfwhm;
		
		double usc = uscore(config, dataMZ, dataInt, mzgrid, nztab, l, h, m);
		double msc = mscore(config, massaxis, masssum, massgrid, l, h, m);
		double zsc = zscore(config, massaxis, masssum, massgrid, l, h, m);
		double fsc = fscore(config, massaxis, masssum, peakx, peaky, l, h, m);

		double dsc = usc * msc * zsc * fsc;
		if (dsc > threshold)
		{
			printf("Peak: Mass:%f Int:%f M:%f U:%f Z:%f F:%f D: %f \n", peakx[i], peaky[i], msc, usc, zsc, fsc, dsc);
			numerator += ival * dsc;
			denominator += ival;
		}
	}

	if (denominator != 0) { avgscore = numerator / denominator; }

	printf("Average Peaks Score: %f\n", avgscore);
	return avgscore;
}