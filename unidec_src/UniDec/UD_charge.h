/*
* UD_charge.h
*
*  Created on : 19 Dec 2017
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//

void z_slice(const float *ivals, const int mindex, float *zdat, const int zlen, const float thresh)
{
	for (int i = 0; i < zlen; i++)
	{
		float val= ivals[index2D(zlen, mindex, i)];
		if (val > thresh) { zdat[i] = val; }
		else {zdat[i] = 0;}
	}
}

void z_slice_range(const float *ivals, const int mindex1, const int mindex2, float *zdat, const int zlen, const float thresh)
{
	for (int i = mindex1; i < mindex2; i++)
	{
		for (int j = 0; j < zlen; j++)
		{
			float val = ivals[index2D(zlen, i, j)];
			if (val > thresh) { zdat[j] += val; }
		}
	}
}

float extract_zmax(Config config, const float peak, const float *xvals, const float *yvals, const float *zvals, const int mlen, const int zlen, const float thresh)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[mlen - 1]) { return 0; }

	float thresh2 = thresh * Max(zvals, mlen * zlen);
	//printf("thresh2 %f %f\n", thresh2, thresh);

	float *zdat = NULL;
	zdat = calloc(zlen, sizeof(float));
	if(config.exwindow<=0){
		int pos = nearfast(xvals, peak, mlen);
		z_slice(zvals, pos, zdat, zlen, thresh2);
	}
	else {
		int pos1 = nearfast(xvals, peak - config.exwindow, mlen);
		int pos2 = nearfast(xvals, peak + config.exwindow, mlen);
		z_slice_range(zvals, pos1, pos2, zdat, zlen, thresh2);
	}
	int maxpos = argmax(zdat, zlen);
	float maxval = yvals[maxpos];
	free(zdat);
	return maxval;
}

float extract_zcom(Config config, const float peak, const float *xvals, const float *yvals, const float *zvals, const int mlen, const int zlen, const float thresh)
{
	if (peak < xvals[0]) { return 0; }
	if (peak > xvals[mlen - 1]) { return 0; }

	float thresh2 = thresh*Max(zvals, mlen*zlen);
	//printf("thresh2 %f %f\n", thresh2, thresh);
	
	int pos = nearfast(xvals, peak, mlen);
	float *zdat = NULL;
	zdat = calloc(zlen, sizeof(float));
	
	if (config.exwindow <= 0) {
		int pos = nearfast(xvals, peak, mlen);
		z_slice(zvals, pos, zdat, zlen, thresh2);
	}
	else {
		int pos1 = nearfast(xvals, peak - config.exwindow, mlen);
		int pos2 = nearfast(xvals, peak + config.exwindow, mlen);
		z_slice_range(zvals, pos1, pos2, zdat, zlen, thresh2);
	}
	float sum = 0;
	float sum_z = 0;
	
	for (int i = 0; i < zlen; i++)
	{
		float x = yvals[i];
		float y = zdat[i];
		//printf("Test\n x %f\n y %f\n", x, y);
		sum_z += y;
		sum += x*y;
	}
	if (sum_z > 0) { sum /= sum_z; }
	free(zdat);
	return sum;
}

float charge_extract_switch(const Config config, const float peak, const float *xvals, const float *yvals, const float *zvals, const int mlen, const int zlen)
{
	float output = 0;
	//config.exwindow = 100000;
	//config.exchoice = 2;
	int swint = config.exchoicez;
	//if (config.exchoice == 4) swint = 0;
	//if (config.exchoice == 3) swint = 1;
	//if (config.exchoice == 5) swint = 3;
	//if (config.exchoice == 6) swint = 2;
	float thresh = config.exthresh/100;

	switch (swint)
	{
	case 0:
		//printf("Extracting Exact Z Max\n");
		output = extract_zmax(config, peak, xvals, yvals, zvals, mlen, zlen, thresh);
		break;
	case 1:
		//printf("Extracting Exact Z COM\n");
		output = extract_zcom(config, peak, xvals, yvals, zvals, mlen, zlen, thresh);
		break;
	/*
	case 2:
		//printf("Extracting Exact Z COM 10 \n");
		output = extract_zcom(config, peak, xvals, yvals, zvals, mlen, zlen, 0.1);
		break;
	case 3:
		//printf("Extracting Exact Z COM 50\n");
		output = extract_zcom(config, peak, xvals, yvals, zvals, mlen, zlen, 0.5);
		break;
	case 4:
		//printf("Extracting Exact Z COM 5\n");
		output = extract_zcom(config, peak, xvals, yvals, zvals, mlen, zlen, 0.05);
		break;
	case 5:
		//printf("Extracting Exact Z COM 2.5\n");
		output = extract_zcom(config, peak, xvals, yvals, zvals, mlen, zlen, 0.025);
		break;
	case 6:
		//printf("Extracting Exact Z COM 1\n");
		output = extract_zcom(config, peak, xvals, yvals, zvals, mlen, zlen, 0.01);
		break;*/
	default:
		printf("Invalid Extraction Choice: %d\n", config.exchoicez);
		output = 0;
	}
	//printf("test %f\n", output);
	return output;
}


void charge_peak_extracts(int argc, char *argv[], Config config, const int ultra)
{
	hid_t file_id;
	file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);

	char dataset[1024];
	char outdat[1024];
	char strval[1024];
	char outdat2[1024];
	char outdat3[1024];
	char pdataset[1024];
	char poutdat[1024];

	strcpy(pdataset, "/peaks");
	if(ultra){ strjoin(pdataset, "/ultrapeakdata", poutdat); }
	else{ strjoin(pdataset, "/peakdata", poutdat); }
	printf("Importing Charge Peaks: %s\n", poutdat);

	int plen = mh5getfilelength(file_id, poutdat);
	float *peakx = NULL;
	peakx = calloc(plen, sizeof(float));
	mh5readfile1d(file_id, poutdat, peakx);
	
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
		strjoin(dataset, "/mass_data", outdat);
		strjoin(dataset, "/charge_data", outdat2);
		strjoin(dataset, "/mass_grid", outdat3);

		int masslen = mh5getfilelength(file_id, outdat);
		int zlen = mh5getfilelength(file_id, outdat2);
		int glen = mh5getfilelength(file_id, outdat3);

		float *mass_ints = NULL;
		float *massaxis = NULL;
		mass_ints = calloc(masslen, sizeof(float));
		massaxis = calloc(masslen, sizeof(float));
		mh5readfile2d(file_id, outdat, masslen, massaxis, mass_ints);
		free(mass_ints);

		float *zints = NULL;
		float *zaxis = NULL;
		zints = calloc(zlen, sizeof(float));
		zaxis = calloc(zlen, sizeof(float));
		mh5readfile2d(file_id, outdat2, zlen, zaxis, zints);
		free(zints);

		float *grid = NULL;
		grid = calloc(glen, sizeof(float));
		mh5readfile1d(file_id, outdat3, grid);
		float sum = 0;
		float max = 0;
		for (int j = 0; j < plen; j++)
		{
			float peak = peakx[j];
			float val = charge_extract_switch(config, peak, massaxis, zaxis, grid, masslen, zlen);
			extracts[index2D(plen, i, j)] = val;
			sum += val;
			if (val > max) { max = val; }
			//printf("Extracts %d %d %f %f %f\n", i, j, max, val, sum);
		}
		//Normalize
		if (max > 0 && config.exnormz == 1)
		{
			for (int j = 0; j < plen; j++)
			{
				extracts[index2D(plen, i, j)] /= max;
				//printf("Extracts %d %d %f %f\n", i, j, extracts[index2D(plen, i, j)], max);
			}
		}
		if (sum > 0 && config.exnormz == 2)
		{
			for (int j = 0; j < plen; j++)
			{
				extracts[index2D(plen, i, j)] /= sum;
				//printf("Extracts %d %d %f\n", i, j, extracts[index2D(plen, i, j)]);
			}
		}
		free(massaxis);
		free(zaxis);
		free(grid);
		
	}
	
	if (config.exnormz == 3 || config.exnormz == 4)
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
			if (config.exnormz == 4) { max = sum; }
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
		strjoin(dataset, "/ultrazextracts", outdat);
	else
		strjoin(dataset, "/zextracts", outdat);
	printf("\tWriting Charge Extracts to: %s\t%d %d %f %f\n", outdat, config.exchoicez, config.exnormz, config.exwindow, config.exthresh);
	mh5writefile2d_grid(file_id, outdat, num, plen, extracts);
	free(peakx);
	free(extracts);
	H5Fclose(file_id);
	//printf("Success!\n");
}