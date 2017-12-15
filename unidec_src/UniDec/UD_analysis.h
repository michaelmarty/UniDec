/*
* MetaUniDec_Main.h
*
*  Created on : 3 June 2017
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//

void interpolate_merge(const double *massaxis, double *outint, const double *tempaxis, const double *tempint, const int mlen, const int templen)
{
	double start = tempaxis[0];
	double end = tempaxis[templen - 1];

	for (int i = 0; i < mlen; i++)
	{
		double pos = massaxis[i];
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
					double mu = (pos - tempaxis[index]) / (tempaxis[index2] - tempaxis[index]);
					double y0 = tempint[index - 1];
					double y1 = tempint[index];
					double y2 = tempint[index2];
					double y3 = tempint[index2 + 1];
					double newval = clip(CubicInterpolate(y0, y1, y2, y3, mu), 0);
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

	hid_t file_id;
	char dataset[1024];
	char outdat[1024];
	char strval[1024];
	
	file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
	int num = 0;
	num = int_attr(file_id, "/ms_dataset", "num", num);
	
	//Read In Data
	strcpy(dataset, "/ms_dataset");
	sprintf(strval, "/%d", 0);
	strcat(dataset, strval);
	strjoin(dataset, dtype, outdat);
	//printf("Processing HDF5 Data Set: %s\n", outdat);
	
	
	int mlen = mh5getfilelength(file_id, outdat);
	
	double *massgrid = NULL;
	double *massaxis = NULL;
	double *masssum = NULL;
	double *temp = NULL;
	
	massgrid = calloc(mlen*num, sizeof(double));
	massaxis = calloc(mlen, sizeof(double));
	masssum = calloc(mlen, sizeof(double));
	temp = calloc(mlen, sizeof(double));
	mh5readfile2d(file_id, outdat, mlen, massaxis, temp);
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

		int templen = mh5getfilelength(file_id, outdat);

		double *temp = NULL;
		double *tempaxis = NULL;
		temp = calloc(templen, sizeof(double));
		tempaxis = calloc(templen, sizeof(double));

		mh5readfile2d(file_id, outdat, templen, tempaxis, temp);

		double *outint = NULL;
		outint = calloc(mlen, sizeof(double));

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
	mh5writefile1d(file_id, outdat, mlen*num, massgrid);

	strjoin(dataset, out2, outdat);
	printf("\tWriting to: %s\n", outdat);
	mh5writefile1d(file_id, outdat, mlen, massaxis);

	strjoin(dataset, out3, outdat);
	printf("\tWriting to: %s\n", outdat);
	if(config.datanorm ==1){
		norm1d(masssum, mlen);
	}
	else {
		double max1 = Max(massgrid, mlen*num);
		double max2 = Max(masssum, mlen);
		if(max2>0){
			for (int i = 0; i<mlen; i++)
			{
				masssum[i] = masssum[i] / max2 * max1;
			}
		}
	}
	mh5writefile1d(file_id, outdat, mlen, masssum);
	free(massgrid);
	free(massaxis);
	free(masssum);
	H5Fclose(file_id);
	clock_t end = clock();
	double totaltime = (double)(end - starttime) / CLOCKS_PER_SEC;
	printf("Done in %f seconds\n", totaltime);

}

int is_peak(const double *dataMZ, const double *dataInt, const int lengthmz, const double window, const double thresh, const int index)
{
	double xval = dataMZ[index];
	double yval = dataInt[index];
	//Check if below threshold
	if (yval < thresh) { return 0; }
	//Check if local max
	for (int i = 0; i < lengthmz; i++)
	{
		double temp = dataMZ[i];
		if (fabs(temp - xval) < window)
		{
			double tempy = dataInt[i];
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

int peak_detect(const double *dataMZ, const double *dataInt, const int lengthmz, const double window, const double thresh, double * peakx, double *peaky)
{
	//printf("Detecting Peaks %d %f %f\n", lengthmz, window, thresh);
	int plen = 0;
	double max = Max(dataInt, lengthmz);
	for (int i = 0; i < lengthmz; i++)
	{
		if (is_peak(dataMZ, dataInt, lengthmz, window, thresh*max, i) == 1)
		{
			printf("Peak %d: %f %f\n", plen, dataMZ[i], dataInt[i]);
			peakx[plen] = dataMZ[i];
			peaky[plen] = dataInt[i];
			plen++;
		}
	}
	return plen;
}


void peak_norm(double *peaky, int plen, int peaknorm)
{
	double norm = 0;
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

double extract_height(Config config, const double peak, const double *xvals, const double *yvals, const int length)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[length - 1]) { return 0; }
	int pos = nearfast(xvals, peak, length);
	return yvals[pos];
}

double extract_localmax(Config config, const double peak, const double *xvals, const double *yvals, const int length)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[length - 1]) { return 0; }

	int pos1 = nearfast(xvals, peak - config.exwindow, length);
	int pos2 = nearfast(xvals, peak + config.exwindow, length);

	double localmax = 0;
	for (int i = pos1; i <= pos2; i++)
	{
		if (yvals[i] > localmax) { localmax = yvals[i]; }
	}

	return localmax;
}

double extract_localmax_position(Config config, const double peak, const double *xvals, const double *yvals, const int length)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[length - 1]) { return 0; }

	int pos1 = nearfast(xvals, peak - config.exwindow, length);
	int pos2 = nearfast(xvals, peak + config.exwindow, length);

	double localmax = 0;
	int localmaxpos = 0;
	for (int i = pos1; i <= pos2; i++)
	{
		if (yvals[i] > localmax) { localmax = yvals[i]; localmaxpos = i; }
	}

	return xvals[localmaxpos];
}

double extract_integral(Config config, const double peak, const double *xvals, const double *yvals, const int length)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[length - 1]) { return 0; }

	int pos1 = nearfast(xvals, peak - config.exwindow, length);
	int pos2 = nearfast(xvals, peak + config.exwindow, length);

	double integral = 0;
	for (int i = pos1+1; i <= pos2; i++)
	{
		double b = xvals[i];
		double a = xvals[i - 1];
		double fb = yvals[i];
		double fa = yvals[i - 1];
		integral += (b - a)*((fa + fb) / 2.0);
	}

	return integral;
}

double extract_center_of_mass(Config config, const double peak, const double *xvals, const double *yvals, const int length, double thresh)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[length - 1]) { return 0; }

	int pos1 = nearfast(xvals, peak - config.exwindow, length);
	int pos2 = nearfast(xvals, peak + config.exwindow, length);
	double sum = 0;
	double sum_masses = 0;
	for (int i = pos1; i <= pos2; i++)
	{
		double x = xvals[i];
		double y = yvals[i];
		if (y > thresh)
		{
			sum_masses += y;
			sum += x*y;
		}
	}
	if (sum_masses > 0) { sum /= sum_masses; }

	return sum;
}


double extract_switch(Config config, const double peak, const double *xvals, const double *yvals, const int length)
{
	double output = 0;
	if (config.exwindow == 0) { config.exchoice = 0; }

	double max = 0;
	if (config.exchoice == 5 || config.exchoice == 6)
	{
		for (int i = 0; i < length; i++)
		{
			if (yvals[i] > max) { max = yvals[i]; }
		}
	}
		

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
		output = extract_integral(config, peak, xvals, yvals, length);
		break;
	case 3:
		//printf("Extracting Center of Mass 0. Window: %f\n", config.exwindow);
		output = extract_center_of_mass(config, peak, xvals, yvals, length, 0);
		break;
	case 4:
		//printf("Extracting Local Max. Window: %f\n", config.exwindow);
		output = extract_localmax_position(config, peak, xvals, yvals, length);
		break;
	case 5:
		//printf("Extracting Center of Mass 50. Window: %f\n", config.exwindow);
		output = extract_center_of_mass(config, peak, xvals, yvals, length, 0.5*max);
		break;
	case 6:
		//printf("Extracting Center of Mass 10. Window: %f\n", config.exwindow);
		output = extract_center_of_mass(config, peak, xvals, yvals, length, 0.1*max);
		break;
	default:
		//printf("Invalid Extraction Choice: %d\n", config.exchoice);
		output = 0;
	}
	return output;
}


void peak_extracts(Config config, const double *peakx, hid_t file_id, const char *dtype, int plen, int ultra)
{
	char dataset[1024];
	char outdat[1024];
	char strval[1024];

	int num = 0;
	num = int_attr(file_id, "/ms_dataset", "num", num);

	double *extracts = NULL;
	extracts = calloc(plen*num, sizeof(double));

	for (int i = 0; i < num; i++)
	{
		strcpy(dataset, "/ms_dataset");
		sprintf(strval, "/%d", i);
		strcat(dataset, strval);
		//printf("Processing HDF5 Data Set: %s\n", dataset);
		strjoin(dataset, dtype, outdat);

		int templen = mh5getfilelength(file_id, outdat);

		double *temp = NULL;
		double *tempaxis = NULL;
		temp = calloc(templen, sizeof(double));
		tempaxis = calloc(templen, sizeof(double));

		mh5readfile2d(file_id, outdat, templen, tempaxis, temp);

		double sum = 0;
		double max = 0;
		for (int j = 0; j < plen; j++)
		{
			double peak = peakx[j];
			double val = extract_switch(config, peak, tempaxis, temp, templen);
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
			double max = 0;
			double sum = 0;
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
	printf("\tWriting Extracts to: %s\t%d %d %f\n", outdat, config.exchoice, config.exnorm, config.exwindow);
	mh5writefile2d_grid(file_id, outdat, num, plen, extracts);

	free(extracts);
}

void get_peaks(int argc, char *argv[], Config config, int ultra)
{

	hid_t file_id;
	char dataset[1024];
	char outdat[1024];
	char strval[1024];

	file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);

	if (!ultra){
	//Read In Data
	strcpy(dataset, "/ms_dataset");
	strjoin(dataset, "/mass_axis", outdat);
	printf("Processing HDF5 Data: %s\n", outdat);

	int mlen = mh5getfilelength(file_id, outdat);

	double *massaxis = NULL;
	double *masssum = NULL;

	massaxis = calloc(mlen, sizeof(double));
	masssum = calloc(mlen, sizeof(double));

	mh5readfile1d(file_id, outdat, massaxis);
	strcpy(dataset, "/ms_dataset");
	strjoin(dataset, "/mass_sum", outdat);
	mh5readfile1d(file_id, outdat, masssum);

	double *peakx=NULL;
	double *peaky=NULL;
	peakx = calloc(mlen, sizeof(double));
	peaky = calloc(mlen, sizeof(double));

	int plen = peak_detect(massaxis, masssum, mlen, config.peakwin, config.peakthresh, peakx, peaky);

	peakx = realloc(peakx, plen*sizeof(double));
	peaky = realloc(peaky, plen*sizeof(double));

	peak_norm(peaky, plen, config.peaknorm);

	strcpy(dataset, "/peaks");
	makegroup(file_id, dataset);
	strjoin(dataset, "/peakdata", outdat);
	printf("\tWriting %d Peaks to: %s\n", plen, outdat);
	mh5writefile2d(file_id, outdat, plen, peakx, peaky);

	
	peak_extracts(config, peakx, file_id, "/mass_data", plen, 0);

	free(peakx);
	free(peaky);
	free(massaxis);
	free(masssum);
	}
	else {
		strcpy(dataset, "/peaks");
		strjoin(dataset, "/ultrapeakdata", outdat);
		printf("Importing Peaks: %s\n", outdat);

		int plen = mh5getfilelength(file_id, outdat);
		double *peakx = NULL;
		peakx = calloc(plen, sizeof(double));
		mh5readfile1d(file_id, outdat, peakx);

		peak_extracts(config, peakx, file_id, "/mass_data", plen, 1);

		free(peakx);

	}
	H5Fclose(file_id);
}

// peak_extracts() extracts for each row in /peakdata and writes to /extracts
//TODO: in peak_extracts() find way to only search for masses from the ultra mass list