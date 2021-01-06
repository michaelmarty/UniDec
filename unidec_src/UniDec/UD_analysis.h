/*
* UD_analysis.h
*
*  Created on : 3 June 2017
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//

#include "UD_dataproc.h"

#ifndef ANALYSIS_HEADER
#define ANALYSIS_HEADER

float single_fwhm(Config config, const int mlen, const float* massaxis, const float* masssum, const float peak, int index, float max);
int ReadDecon(Config* config, const Input inp, Decon* decon);
float score_from_peaks(const int plen, const float* peakx, const float* peaky, float* dscores, const Config config, Decon* decon, const Input inp, const float threshold);
float score(Config config, Decon* decon, Input inp, const float threshold);

void interpolate_merge(const float *massaxis, float *outint, const float *tempaxis, const float *tempint, const int mlen, const int templen)
{
	float start = tempaxis[0];
	float end = tempaxis[templen - 1];

	for (int i = 0; i < mlen; i++)
	{
		float pos = massaxis[i];
		if (pos >= start && pos <= end)
		{
			int index = nearfast(tempaxis, pos, templen);
			int index2 = index;
			if (tempaxis[index] == pos)
			{
				outint[i] = tempint[index];
			}
			else
			{
				if (tempaxis[index]>pos&&index>1 && index<templen - 1)
				{
					index2 = index;
					index = index - 1;
				}
				else if (tempaxis[index]<pos&&index<templen - 2 && index>0)
				{
					index2 = index + 1;
				}
				if (index2>index && (tempaxis[index2] - tempaxis[index]) != 0)
				{
					float mu = (pos - tempaxis[index]) / (tempaxis[index2] - tempaxis[index]);
					float y0 = tempint[index - 1];
					float y1 = tempint[index];
					float y2 = tempint[index2];
					float y3 = tempint[index2 + 1];
					float newval = clip(CubicInterpolate(y0, y1, y2, y3, mu), 0);
					//newval=CRSplineInterpolate(y0,y1,y2,y3,mu);
					//newval=LinearInterpolate(y1,y2,mu);
					outint[i] = newval;
				}
			}
		}
	}

}



void make_grid(int argc, char *argv[], Config config, const char *dtype, const char *out1, const char *out2, const char *out3)
{
	clock_t starttime;
	starttime = clock();


	char dataset[1024];
	char outdat[1024];
	char strval[1024];
	
	config.file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
	int num = 0;
	num = int_attr(config.file_id, "/ms_dataset", "num", num);
	
	//Read In Data
	strcpy(dataset, "/ms_dataset");
	sprintf(strval, "/%d", 0);
	strcat(dataset, strval);
	strjoin(dataset, dtype, outdat);
	//printf("Processing HDF5 Data Set: %s\n", outdat);
	
	
	int mlen = mh5getfilelength(config.file_id, outdat);
	
	float *massgrid = NULL;
	float *massaxis = NULL;
	float *masssum = NULL;
	float *temp = NULL;
	
	massgrid = calloc(mlen*num, sizeof(float));
	massaxis = calloc(mlen, sizeof(float));
	masssum = calloc(mlen, sizeof(float));
	temp = calloc(mlen, sizeof(float));
	mh5readfile2d(config.file_id, outdat, mlen, massaxis, temp);
	for (int i = 0; i < mlen; i++)
	{
		massgrid[i] = temp[i];
		masssum[i] = temp[i];
	}
	free(temp);

	for (int i = 1; i < num; i++)
	{
		strcpy(dataset, "/ms_dataset");
		sprintf(strval, "/%d", i);
		strcat(dataset, strval);
		//printf("Processing HDF5 Data Set: %s\n", dataset);
		strjoin(dataset, dtype, outdat);

		int templen = mh5getfilelength(config.file_id, outdat);

		float *temp = NULL;
		float *tempaxis = NULL;
		temp = calloc(templen, sizeof(float));
		tempaxis = calloc(templen, sizeof(float));

		mh5readfile2d(config.file_id, outdat, templen, tempaxis, temp);

		float *outint = NULL;
		outint = calloc(mlen, sizeof(float));

		interpolate_merge(massaxis, outint, tempaxis, temp, mlen, templen);

		for (int j = 0; j < mlen; j++)
		{
			massgrid[i*mlen+j] = outint[j];
			masssum[j] += outint[j];
		}

		free(temp);
		free(tempaxis);
		free(outint);
	}

	//Write data to processed_data
	strcpy(dataset, "/ms_dataset");
	strjoin(dataset, out1, outdat);
	printf("\tWriting to: %s\n", outdat);
	mh5writefile1d(config.file_id, outdat, mlen*num, massgrid);

	strjoin(dataset, out2, outdat);
	printf("\tWriting to: %s\n", outdat);
	mh5writefile1d(config.file_id, outdat, mlen, massaxis);

	strjoin(dataset, out3, outdat);
	printf("\tWriting to: %s\n", outdat);
	if(config.datanorm ==1){
		norm1d(masssum, mlen);
	}
	else {
		float max1 = Max(massgrid, mlen*num);
		float max2 = Max(masssum, mlen);
		if(max2>0){
			for (int i = 0; i<mlen; i++)
			{
				masssum[i] = masssum[i] / max2 * max1;
			}
		}
	}
	mh5writefile1d(config.file_id, outdat, mlen, masssum);
	free(massgrid);
	free(massaxis);
	free(masssum);
	H5Fclose(config.file_id);
	clock_t end = clock();
	float totaltime = (float)(end - starttime) / CLOCKS_PER_SEC;
	printf("Done in %f seconds\n", totaltime);

}

int is_peak(const float *dataMZ, const float *dataInt, const int lengthmz, const float window, const float thresh, const int index)
{
	float xval = dataMZ[index];
	float yval = dataInt[index];
	//Check if below threshold
	if (yval < thresh) { return 0; }
	//Check if local max
	for (int i = 0; i < lengthmz; i++)
	{
		float temp = dataMZ[i];
		if (fabs(temp - xval) <= window)
		{
			float tempy = dataInt[i];
			//printf("%f %f %f %f\n", temp, tempy, xval, yval);
			if (tempy > yval)
			{
				return 0;
			}
			if (tempy == yval && i < index)
			{
				return 0;
			}
		}
	}
	
	return 1;
}

int peak_detect(const float *dataMZ, const float *dataInt, const int lengthmz, const float window, const float thresh, float * peakx, float *peaky)
{
	//printf("Detecting Peaks %d %f %f\n", lengthmz, window, thresh);
	int plen = 0;
	float max = Max(dataInt, lengthmz);
	for (int i = 0; i < lengthmz; i++)
	{
		if (is_peak(dataMZ, dataInt, lengthmz, window, thresh*max, i) == 1)
		{
			//printf("Peak %d: %f %f\n", plen, dataMZ[i], dataInt[i]);//
			peakx[plen] = dataMZ[i];
			peaky[plen] = dataInt[i];
			plen++;
		}
	}
	return plen;
}


void peak_norm(float *peaky, int plen, int peaknorm)
{
	float norm = 0;
	if (peaknorm==1){norm= Max(peaky, plen);}
	if (peaknorm==2){norm= Sum(peaky, plen);}
	if (norm != 0)
	{
		for (int i = 0; i < plen; i++)
		{
			peaky[i] /= norm;
		}
	}
}

float extract_height(Config config, const float peak, const float *xvals, const float *yvals, const int length)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[length - 1]) { return 0; }
	int pos = nearfast(xvals, peak, length);
	return yvals[pos];
}

float extract_localmax(Config config, const float peak, const float *xvals, const float *yvals, const int length)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[length - 1]) { return 0; }

	int pos1 = nearfast(xvals, peak - config.exwindow, length);
	int pos2 = nearfast(xvals, peak + config.exwindow, length);

	float localmax = 0;
	for (int i = pos1; i <= pos2; i++)
	{
		if (yvals[i] > localmax) { localmax = yvals[i]; }
	}

	return localmax;
}

float extract_localmax_position(Config config, const float peak, const float *xvals, const float *yvals, const int length)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[length - 1]) { return 0; }

	int pos1 = nearfast(xvals, peak - config.exwindow, length);
	int pos2 = nearfast(xvals, peak + config.exwindow, length);

	float localmax = 0;
	int localmaxpos = 0;
	for (int i = pos1; i <= pos2; i++)
	{
		if (yvals[i] > localmax) { localmax = yvals[i]; localmaxpos = i; }
	}

	return xvals[localmaxpos];
}

float extract_integral(Config config, const float peak, const float *xvals, const float *yvals, const int length, const float thresh)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[length - 1]) { return 0; }

	float thresh2 = 0;
	if (thresh > 0) {
		float max = Max(yvals, length);
		thresh2 = thresh * max;
	}
	//printf("thresh %f\n", thresh2);

	int pos1 = nearfast(xvals, peak - config.exwindow, length);
	int pos2 = nearfast(xvals, peak + config.exwindow, length);

	float integral = 0;
	for (int i = pos1+1; i <= pos2; i++)
	{
		float b = xvals[i];
		float a = xvals[i - 1];
		float fb = yvals[i];
		float fa = yvals[i - 1];
		if (fa > thresh2 && fb > thresh2) {
			integral += (b - a) * ((fa + fb) / 2.0);
		}
	}

	return integral;
}


float extract_center_of_mass(Config config, const float peak, const float *xvals, const float *yvals, const int length, const float thresh)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[length - 1]) { return 0; }
		
	float thresh2 = 0;
	if (thresh > 0) {
		float max = Max(yvals, length);
		thresh2 = thresh * max;
	}
	//printf("thresh %f\n", thresh2);

	int pos1 = nearfast(xvals, peak - config.exwindow, length);
	int pos2 = nearfast(xvals, peak + config.exwindow, length);
	float sum = 0;
	float sum_masses = 0;
	for (int i = pos1; i <= pos2; i++)
	{
		float x = xvals[i];
		float y = yvals[i];
		if (y > thresh2)
		{
			sum_masses += y;
			sum += x*y;
		}
	}
	if (sum_masses > 0) { sum /= sum_masses; }

	return sum;
}


float extract_estimated_area(Config config, const float peak, const float* xvals, const float* yvals, const int length)
{
	int pos1 = nearfast(xvals, peak - config.exwindow, length);
	// printf("Position of left bound is %d\n", pos1);
	int pos2 = nearfast(xvals, peak + config.exwindow, length);
	// printf("Position of right bound is %d\n", pos2);
	int mlen = pos2 - pos1 + 1;
	const float* xwin = xvals + pos1;
	// printf("xvals is %p, xwin is %p\n", xvals, xwin);
	const float* ywin = yvals + pos1;

	int index = nearfast(xwin, peak, mlen);
	float height = ywin[index];
	float fwhm = single_fwhm(config, mlen, xwin, ywin, peak, index, height);
	
	float pi = 3.14159265358979323846;
	float gauss_coeff = sqrt(pi / log(2.0)) / 2.0;
	float adjusted_coeff = ((0.5 * gauss_coeff) + (pi / 4.0));
	float area = 0;
	if (config.psfun == 0) { // Gaussian
		area = height * fwhm * gauss_coeff;
	}
	else if (config.psfun == 1) { // Lorentzian
		area = height * fwhm * pi / 2.0;
	}
	else if (config.psfun == 2) { // Split G/L
		area = height * fwhm * adjusted_coeff;
	}

	return area;
}


float extract_switch(Config config, const float peak, const float* xvals, const float* yvals, const int length)
{
	float output = 0;
	if (config.exwindow == 0) { config.exchoice = 0; }
	
	float thresh = config.exthresh / 100;

	switch (config.exchoice)
	{
	case 0:
		//printf("Extracting Height\n");
		output = extract_height(config, peak, xvals, yvals, length);
		break;
	case 1:
		//printf("Extracting Local Max. Window: %f\n", config.exwindow);
		output = extract_localmax(config, peak, xvals, yvals, length);
		break;
	case 2:
		//printf("Extracting Integral. Window: %f\n", config.exwindow);
		output = extract_integral(config, peak, xvals, yvals, length, thresh);
		break;
	case 3:
		//printf("Extracting Center of Mass 0. Window: %f\n", config.exwindow);
		output = extract_center_of_mass(config, peak, xvals, yvals, length, thresh);
		break;
	case 4:
		//printf("Extracting Local Max. Window: %f\n", config.exwindow);
		output = extract_localmax_position(config, peak, xvals, yvals, length);
		break;
	/*
	case 5:
		//printf("Extracting Center of Mass 50. Window: %f\n", config.exwindow);
		output = extract_center_of_mass(config, peak, xvals, yvals, length, 0.5*max);
		break;
	case 6:
		//printf("Extracting Center of Mass 10. Window: %f\n", config.exwindow);
		output = extract_center_of_mass(config, peak, xvals, yvals, length, 0.1*max);
		break;
	case 7:
		//printf("Extracting Integral. Window: %f\n", config.exwindow);
		output = extract_integral(config, peak, xvals, yvals, length, 0.5 * max);
		break;
	case 8:
		//printf("Extracting Integral. Window: %f\n", config.exwindow);
		output = extract_integral(config, peak, xvals, yvals, length, 0.1 * max);
		break;
	case 9:
		//printf("Extracting Integral. Window: %f\n", config.exwindow);
		output = extract_integral(config, peak, xvals, yvals, length, 0.05 * max);
	   	break;
	case 10:
		//printf("Extracting Integral. Window: %f\n", config.exwindow);
		output = extract_integral(config, peak, xvals, yvals, length, 0.025 * max);
		break;*/
	case 5: // 11
		// printf("Extracting Estimated Area. Window: %f\n", config.exwindow);
		output = extract_estimated_area(config, peak, xvals, yvals, length);
		break;
	default:
		printf("Invalid Extraction Choice: %d\n", config.exchoice);
		output = 0;
	}
	return output;
}


void peak_extracts(Config config, const float *peakx, hid_t file_id, const char *dtype, int plen, int ultra)
{
	char dataset[1024];
	char outdat[1024];
	char strval[1024];

	int num = 0;
	num = int_attr(file_id, "/ms_dataset", "num", num);

	float *extracts = NULL;
	extracts = calloc(plen*num, sizeof(float));

	for (int i = 0; i < num; i++)
	{
		strcpy(dataset, "/ms_dataset");
		sprintf(strval, "/%d", i);
		strcat(dataset, strval);
		//printf("Processing HDF5 Data Set: %s\n", dataset);
		strjoin(dataset, dtype, outdat);

		int templen = mh5getfilelength(file_id, outdat);

		float *temp = NULL;
		float *tempaxis = NULL;
		temp = calloc(templen, sizeof(float));
		tempaxis = calloc(templen, sizeof(float));

		mh5readfile2d(file_id, outdat, templen, tempaxis, temp);

		float sum = 0;
		float max = 0;
		for (int j = 0; j < plen; j++)
		{
			float peak = peakx[j];
			float val = extract_switch(config, peak, tempaxis, temp, templen);
			extracts[index2D(plen, i, j)] = val;
			sum += val;
			if (val > max) {max = val;}
			//printf("Extracts %d %d %f %f %f\n", i, j, max, val, sum);
		}

		//Normalize
		if (max > 0 && config.exnorm == 1)
		{
			for (int j = 0; j < plen; j++)
			{
				extracts[index2D(plen, i, j)] /= max;
				//printf("Extracts %d %d %f %f\n", i, j, extracts[index2D(plen, i, j)], max);
			}
		}
		if (sum > 0 && config.exnorm == 2)
		{
			for (int j = 0; j < plen; j++)
			{
				extracts[index2D(plen, i, j)] /= sum;
				//printf("Extracts %d %d %f\n", i, j, extracts[index2D(plen, i, j)]);
			}
		}

		free(temp);
		free(tempaxis);
	}

	if (config.exnorm == 3|| config.exnorm==4)
	{
		for (int j = 0; j < plen; j++)
		{
			float max = 0;
			float sum = 0;
			for (int i = 0; i < num; i++)
			{
				if (extracts[index2D(plen, i, j)] > max) { max = extracts[index2D(plen, i, j)]; }
				sum += extracts[index2D(plen, i, j)];
			}
			if (config.exnorm == 4) { max = sum; }
			if (max > 0)
			{
				for (int i = 0; i < num; i++)
				{
					extracts[index2D(plen, i, j)] /= max;
				}
			}
		}
	}

	strcpy(dataset, "/peaks");
	makegroup(file_id, dataset);
	if (ultra)
		strjoin(dataset, "/ultraextracts", outdat);
	else
		strjoin(dataset, "/extracts", outdat);
	printf("\tWriting Extracts to: %s\t%d %d %f %f\n", outdat, config.exchoice, config.exnorm, config.exwindow, config.exthresh);
	mh5writefile2d_grid(file_id, outdat, num, plen, extracts);

	free(extracts);
}


void get_all_peaks(int argc, char *argv[], Config config)
{
	hid_t file_id;
	char dataset[1024];
	char outdat[1024];
	char strval[1024];

	file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
	int num = 0;
	num = int_attr(file_id, "/ms_dataset", "num", num);

	for (int i = 0; i < num; i++) {
		//Read In Data
		strcpy(dataset, "/ms_dataset");
		sprintf(strval, "/%d", i);
		strcat(dataset, strval);

		strjoin(dataset, "/mass_data", outdat);
		printf("Processing HDF5 Data: %s\n", outdat);

		int mlen = mh5getfilelength(file_id, outdat);

		float *massaxis = NULL;
		float *masssum = NULL;

		massaxis = calloc(mlen, sizeof(float));
		masssum = calloc(mlen, sizeof(float));

		mh5readfile2d(file_id, outdat, mlen, massaxis, masssum);

		float *peakx = NULL;
		float *peaky = NULL;
		peakx = calloc(mlen, sizeof(float));
		peaky = calloc(mlen, sizeof(float));

		int plen = peak_detect(massaxis, masssum, mlen, config.peakwin, config.peakthresh, peakx, peaky);

		peakx = realloc(peakx, plen * sizeof(float));
		peaky = realloc(peaky, plen * sizeof(float));

		peak_norm(peaky, plen, config.peaknorm);

		strjoin(dataset, "/peakdata", outdat);
		printf("\tWriting %d Peaks to: %s\n", plen, outdat);
		mh5writefile2d(file_id, outdat, plen, peakx, peaky);

		free(peakx);
		free(peaky);
		free(massaxis);
		free(masssum);
	}
	H5Fclose(file_id);
}

void get_peaks(int argc, char *argv[], Config config, int ultra)
{
	char dataset[1024];
	char outdat[1024];
	char strval[1024];

	config.file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);

	if (!ultra){
		//Read In Data
		strcpy(dataset, "/ms_dataset");
		strjoin(dataset, "/mass_axis", outdat);
		printf("Processing HDF5 Data: %s\n", outdat);

		int mlen = mh5getfilelength(config.file_id, outdat);

		float *massaxis = NULL;
		float *masssum = NULL;

		massaxis = calloc(mlen, sizeof(float));
		masssum = calloc(mlen, sizeof(float));

		mh5readfile1d(config.file_id, outdat, massaxis);
		strcpy(dataset, "/ms_dataset");
		strjoin(dataset, "/mass_sum", outdat);
		mh5readfile1d(config.file_id, outdat, masssum);

		float *peakx=NULL;
		float *peaky=NULL;
		peakx = calloc(mlen, sizeof(float));
		peaky = calloc(mlen, sizeof(float));
		float* dscores = NULL;
		float* numerators = NULL;
		float* denominators = NULL;

		int plen = peak_detect(massaxis, masssum, mlen, config.peakwin, config.peakthresh, peakx, peaky);

		peakx = realloc(peakx, plen*sizeof(float));
		peaky = realloc(peaky, plen*sizeof(float));
		dscores = calloc(plen, sizeof(float));
		numerators = calloc(plen, sizeof(float));
		denominators = calloc(plen, sizeof(float));

		peak_norm(peaky, plen, config.peaknorm);

		int num = 0;
		num = int_attr(config.file_id, "/ms_dataset", "num", num);

		// Get the dscores by scoring the peak at each scan
		for (int i = 0; i < num; i++) {
			//Create a temp array
			float* tempdscores = NULL;
			tempdscores = calloc(plen, sizeof(float));
			
			//Import everything
			config.metamode = i;
			Decon decon = SetupDecon();
			Input inp = SetupInputs();
			ReadInputs(argc, argv, &config, &inp);
			int status = ReadDecon(&config, inp, &decon);

			//Check to make sure key deconvolution grids are there
			if (status == 1) {
				//Get the scores for each peak
				//score(config, &decon, inp, 0);
				score_from_peaks(plen, peakx, peaky, tempdscores, config, &decon, inp, 0);

				//Average In Dscores
				for (int j = 0; j < plen; j++)
				{
					//Get the peaky value locally
					int index = nearfast(decon.massaxis, peakx[j], decon.mlen);
					float yval = decon.massaxisval[index];
					numerators[j] += yval * tempdscores[j];
					denominators[j] += yval;
				}
			}
			else
			{
				printf("Missing deconvolution outputs. Turn off Fast Profile/Fast Centroid and try deconvolving again.");
			}

			//Free things
			FreeDecon(decon);
			FreeInputs(inp);
			free(tempdscores);
		}

		//Average In Dscores
		for (int j = 0; j < plen; j++)
		{
			if(denominators[j] != 0){dscores[j] = numerators[j]/denominators[j];}
		}

		//Write the outputs
		strcpy(dataset, "/peaks");
		makegroup(config.file_id, dataset);
		strjoin(dataset, "/peakdata", outdat);
		printf("\tWriting %d Peaks to: %s\n", plen, outdat);

		float* ptemp = NULL;
		ptemp = calloc(plen * 3, sizeof(float));

		for (int i = 0; i < plen; i++) {
			ptemp[i * 3] = peakx[i];
			ptemp[i * 3 + 1] = peaky[i];
			ptemp[i * 3 + 2] = dscores[i];
		}

		mh5writefile2d_grid(config.file_id, outdat, plen, 3, ptemp);
		free(ptemp);
		//mh5writefile2d(config.file_id, outdat, plen, peakx, peaky);
	
		peak_extracts(config, peakx, config.file_id, "/mass_data", plen, 0);

		free(peakx);
		free(peaky);
		free(massaxis);
		free(masssum);
		free(dscores);
		free(numerators);
		free(denominators);
	}
	else {
		strcpy(dataset, "/peaks");
		strjoin(dataset, "/ultrapeakdata", outdat);
		printf("Importing Peaks: %s\n", outdat);

		int plen = mh5getfilelength(config.file_id, outdat);
		float *peakx = NULL;
		peakx = calloc(plen, sizeof(float));
		mh5readfile1d(config.file_id, outdat, peakx);

		peak_extracts(config, peakx, config.file_id, "/mass_data", plen, 1);

		free(peakx);

	}
	H5Fclose(config.file_id);
}

// peak_extracts() extracts for each row in /peakdata and writes to /extracts

/*
int peak_combine(int argc, char* argv[], float* peakx, float* peaky, Config config)
{
	//Unfinished
	char dataset[1024];
	char outdat[1024];
	char strval[1024];

	config.file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);

	int plen = 1;
	float* peakx = NULL;
	float* peaky = NULL;
	float* dscores = NULL;

	peakx = calloc(plen, sizeof(float));
	peaky = calloc(plen, sizeof(float));
	dscores = calloc(plen, sizeof(float));
	
	float* numerators = NULL;
	float* denominators = NULL;

	//peakx = realloc(peakx, plen * sizeof(float));
	//peaky = realloc(peaky, plen * sizeof(float));
	
	//numerators = calloc(plen, sizeof(float));
	//denominators = calloc(plen, sizeof(float));

	//peak_norm(peaky, plen, config.peaknorm);

	int num = 0;
	num = int_attr(config.file_id, "/ms_dataset", "num", num);

	// Get the dscores by scoring the peak at each scan
	for (int i = 0; i < num; i++) {
		//Create a temp array
		float* tempdscores = NULL;
		tempdscores = calloc(plen, sizeof(float));

		//Import everything
		config.metamode = i;
		Decon decon = SetupDecon();
		Input inp = SetupInputs();
		ReadInputs(argc, argv, &config, &inp);
		int status = ReadDecon(&config, inp, &decon);

		//Check to make sure key deconvolution grids are there
		if (status == 1) {
			//Get the scores for each peak
			//score(config, &decon, inp, 0);
			score_from_peaks(plen, peakx, peaky, tempdscores, config, &decon, inp, 0);

			//Average In Dscores
			for (int j = 0; j < plen; j++)
			{
				//Get the peaky value locally
				int index = nearfast(decon.massaxis, peakx[j], decon.mlen);
				float yval = decon.massaxisval[index];
				numerators[j] += yval * tempdscores[j];
				denominators[j] += yval;
			}
		}
		else
		{
			printf("Missing deconvolution outputs. Turn off Fast Profile/Fast Centroid and try deconvolving again.");
		}

		//Free things
		FreeDecon(decon);
		FreeInputs(inp);
		free(tempdscores);
	}

	//Average In Dscores
	for (int j = 0; j < plen; j++)
	{
		if (denominators[j] != 0) { dscores[j] = numerators[j] / denominators[j]; }
	}

	//Write the outputs
	strcpy(dataset, "/peaks");
	makegroup(config.file_id, dataset);
	strjoin(dataset, "/peakdata", outdat);
	printf("\tWriting %d Peaks to: %s\n", plen, outdat);

	float* ptemp = NULL;
	ptemp = calloc(plen * 3, sizeof(float));

	for (int i = 0; i < plen; i++) {
		ptemp[i * 3] = peakx[i];
		ptemp[i * 3 + 1] = peaky[i];
		ptemp[i * 3 + 2] = dscores[i];
	}

	mh5writefile2d_grid(config.file_id, outdat, plen, 3, ptemp);
	free(ptemp);
	//mh5writefile2d(config.file_id, outdat, plen, peakx, peaky);

	peak_extracts(config, peakx, config.file_id, "/mass_data", plen, 0);

	free(peakx);
	free(peaky);
	free(massaxis);
	free(masssum);
	free(dscores);
	free(numerators);
	free(denominators);
	
	H5Fclose(config.file_id);
}*/



#endif