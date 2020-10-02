#pragma once
/*
* UD_dataproc.h
*
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//
#include <stdlib.h>

#ifndef DATAPROC_HEADER
#define DATAPROC_HEADER

//Pool data
int pool1d(double *oldmz, double *oldint, const int oldlen, const int mzbins)
{
	double *newmz, *newint;
	int newlen = oldlen/mzbins;

	if (oldlen > newlen*mzbins)
	{
		newlen++;
	}
	
	//declare memory for new arrays
	newmz = calloc(newlen, sizeof(double));
	newint = calloc(newlen, sizeof(double));
	
	// printf("%d %d %d\n", oldlen, newlen, mzbins);
	for (int i = 0; i < newlen; i++)
	{
		double mz = 0;
		double val = 0;
		double bins = 0;
		for(int j=0; j<mzbins; j++)
		{
			int index = i*mzbins + j;
			if (index < oldlen)
			{
				mz += oldmz[index];
				val += oldint[index];
				bins += 1;
			}
		}
		newmz[i] = mz / bins;
		newint[i] = val/ bins;
	}

	//Bookkeeping 
	realloc(oldmz, newlen * sizeof(double));
	realloc(oldint, newlen * sizeof(double));

	memcpy(oldmz, newmz, sizeof(double)*newlen);
	memcpy(oldint, newint, sizeof(double)*newlen);

	free(newmz);
	free(newint);

	return newlen;
}


//Chop data
int chop1d(double *oldmz, double *oldint, int oldlen, const double min, const double max)
{
	double *newmz, *newint;
	int newlen = 0;

	double maxval = max1d(oldmz, oldlen);
	double minval = min1d(oldmz, oldlen);
	if (maxval<max && minval>min)
	{
		return oldlen;
	}

	//get new length
	for (int i = 0; i < oldlen; i++)
	{
		if (oldmz[i] >= min && oldmz[i] <= max)
		{
			newlen++;
		}
	}

	if (newlen >= oldlen)
	{
		return oldlen;
	}
	
	//declare memory for new arrays
	newmz = calloc(newlen, sizeof(double));
	newint = calloc(newlen, sizeof(double));
	
	//Fill new arrays
	newlen = 0;
	for (int i = 0; i < oldlen; i++)
	{
		if (oldmz[i] >= min &&oldmz[i] <= max)
		{
			newmz[newlen] = oldmz[i];
			newint[newlen] = oldint[i];
			newlen++;
		}
	}
	
	//Bookkeeping 
	realloc(oldmz, newlen * sizeof(double));
	realloc(oldint, newlen * sizeof(double));
	
	memcpy(oldmz, newmz, sizeof(double)*newlen);
	memcpy(oldint, newint, sizeof(double)*newlen);
	
	free(newmz);
	free(newint);
	
	return newlen;	
}

//Remove duplicates
int remove_duplicates(double *oldmz, double *oldint, int oldlen)
{
	double *newmz, *newint;
	int newlen = 1;
	double val = 0;


	//get new length
	for (int i = 0; i < oldlen-1; i++)
	{
		if (oldmz[i] != oldmz[i+1])
		{
			newlen++;
		}
	}

	if (newlen >= oldlen)
	{
		return oldlen;
	}

	//declare memory for new arrays
	newmz = calloc(newlen, sizeof(double));
	newint = calloc(newlen, sizeof(double));

	//Fill new arrays
	newlen = 0;
	for (int i = 0; i < oldlen - 1; i++)
	{
		if (oldmz[i] != oldmz[i + 1])
		{
			newmz[newlen] = oldmz[i];
			newint[newlen] = oldint[i]+val;
			newlen++;
			val = 0;
		}
		else
		{
			val += oldint[i];
		}
	}
	newmz[newlen] = oldmz[oldlen-1];
	newint[newlen] = oldint[oldlen-1] + val;
	newlen++;

	//Bookkeeping 
	realloc(oldmz, newlen * sizeof(double));
	realloc(oldint, newlen * sizeof(double));

	memcpy(oldmz, newmz, sizeof(double)*newlen);
	memcpy(oldint, newint, sizeof(double)*newlen);

	//free(newmz);
	//free(newint);
	//printf("Removed %d Duplicates\n", oldlen - newlen);
	return newlen;
}

//Remove duplicates
int remove_middle_zeros(double *oldmz, double *oldint, int oldlen)
{
	double *newmz, *newint;
	int newlen = 2;

	//get new length
	for (int i = 1; i < oldlen - 1; i++)
	{
		if (oldint[i] != 0 ||  oldint[i - 1] != 0 || oldint[i + 1] != 0)
		{
			newlen++;
		}
	}

	if (newlen >= oldlen)
	{
		return oldlen;
	}

	//declare memory for new arrays
	newmz = calloc(newlen, sizeof(double));
	newint = calloc(newlen, sizeof(double));

	//Fill new arrays
	newlen = 0;
	newmz[newlen] = oldmz[0];
	newint[newlen] = oldint[0];
	newlen++;
	for (int i = 1; i < oldlen - 1; i++)
	{
		if (oldint[i] != 0 || oldint[i - 1] != 0 || oldint[i + 1] != 0)
		{
			newmz[newlen] = oldmz[i];
			newint[newlen] = oldint[i];
			newlen++;

		}
	}
	newmz[newlen] = oldmz[oldlen-1];
	newint[newlen] = oldint[oldlen-1];
	newlen++;

	//Bookkeeping 
	realloc(oldmz, newlen * sizeof(double));
	realloc(oldint, newlen * sizeof(double));

	memcpy(oldmz, newmz, sizeof(double)*newlen);
	memcpy(oldint, newint, sizeof(double)*newlen);

	//free(newmz);
	//free(newint);
	//printf("Removed %d Middle Zeros\n", oldlen - newlen);
	return newlen;
}

//Normalize
void norm1d(double *dataInt, const int lengthmz)
{
	double maxval = max1d(dataInt, lengthmz);
	if(maxval>0){
		for (int i = 0; i < lengthmz; i++)
		{
			dataInt[i] /= maxval;
		}
	}
}

void background_subtract(double *dataInt, const double bsub, const int lengthmz)
{
	double *background;
	background = calloc(lengthmz, sizeof(double));
	// Get Local Mins
	for (int i=0; i < lengthmz; i++)
	{
		int start = i - abs(bsub);
		int end = i + abs(bsub);
		if (start < 0) { start = 0; }
		if (end > lengthmz) { end = lengthmz; }
		double min = dataInt[start];
		for (int j = start+1; j < end; j++)
		{
			if (dataInt[j] < min) { min = dataInt[j]; }
		}
		background[i] = min;
	}
	// Smooth Local Mins and Subtract
	int threshold = abs(bsub) * 6;
	double *gaussian;
	gaussian = calloc(threshold, sizeof(double));
	for (int i = 0; i < threshold; i++)
	{
		gaussian[i] = ndis(i, abs(bsub) * 3, abs(bsub));
	}
	int i;
	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i < lengthmz; i++)
	{
		int start = i - abs(bsub)*3;
		int end = i + abs(bsub)*3;
		if (start < 0) { start = 0; }
		if (end > lengthmz) { end = lengthmz; }
		double val = 0;
		for (int j = start; j < end; j++)
		{
			val += background[j] * gaussian[j - start];
		}
		dataInt[i] -= val;
	}

	free(background);
	free(gaussian);
}

// A comparator function used by qsort 
inline int compare(const void* a, const void* b)
{
	return (*(double*)a - *(double*)b);
}

//Cut out the lowest x percent of the data
int data_reduction(double* oldmz, double* oldint, int oldlen, const double redper)
{
	double* newmz, * newint, *sortint, cutoff;
	int newlen = 0;
	
	//Create the sort array and load it with intensities
	sortint = calloc(oldlen, sizeof(double));
	memcpy(sortint, oldint, sizeof(double) * oldlen);
	//Sort it
	qsort(sortint, oldlen, sizeof(double), compare);
	
	//Find the index redper percent of the way through the list
	int index = (int)redper * oldlen / 100.;
	if (index < 0) { index = 0; }
	if (index >= oldlen) { index = oldlen - 2; }
	
	// Get the cutoff at that index
	cutoff = sortint[index];
	//printf("Cutoff: %f\n", cutoff);
	
	//get new length
	for (int i = 0; i < oldlen; i++)
	{
		if (oldint[i] >= cutoff)
		{
			newlen++;
		}
	}

	//Some checks
	if (newlen >= oldlen){return oldlen;}
	if (newlen < 3) { printf("Warning: Aggressive Data Reduction!!! Lower the data reduction percentage\n\n"); }

	//declare memory for new arrays
	newmz = calloc(newlen, sizeof(double));
	newint = calloc(newlen, sizeof(double));

	//Fill new arrays
	newlen = 0;
	for (int i = 0; i < oldlen; i++)
	{
		if (oldint[i] >= cutoff)
		{
			newmz[newlen] = oldmz[i];
			newint[newlen] = oldint[i];
			newlen++;
		}
	}

	//Bookkeeping 
	realloc(oldmz, newlen * sizeof(double));
	realloc(oldint, newlen * sizeof(double));

	memcpy(oldmz, newmz, sizeof(double) * newlen);
	memcpy(oldint, newint, sizeof(double) * newlen);

	free(newmz);
	free(newint);
	free(sortint);
	return newlen;
}


void process_data(int argc, char *argv[], Config config)
{	
	hid_t file_id;
	double *dataMZ = NULL;
	double *dataInt = NULL;
	char dataset[1024];
	char outdat[1024];
	char strval[1024];

	//Read In Data
	strcpy(dataset, "/ms_dataset");
	sprintf(strval, "/%d", config.metamode);
	strcat(dataset, strval);
	printf("Processing HDF5 Data Set: %s\n", dataset);
	file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
	strjoin(dataset, "/raw_data", outdat);
	int lengthmz = mh5getfilelength(file_id, outdat);
	dataMZ = calloc(lengthmz, sizeof(double));
	dataInt = calloc(lengthmz, sizeof(double));
	mh5readfile2d(file_id, outdat, lengthmz, dataMZ, dataInt);
	//printf("Length of Data: %d \n", lengthmz);


	//Bin
	if (config.mzbins > 1)
	{
		//printf("Binning data every %d data points\n", config.mzbins);
		lengthmz = pool1d(dataMZ, dataInt, lengthmz, config.mzbins);
	}

	//Chop
	if (config.minmz >= 0 && config.maxmz > 0)
	{
		//printf("Chopping data to range: %f %f\n", config.minmz, config.maxmz);
		lengthmz=chop1d(dataMZ, dataInt, lengthmz, config.minmz, config.maxmz);
	}

	//Remove Zeros
	//printf("Removing zeros %d\n", lengthmz);
	lengthmz = remove_middle_zeros(dataMZ, dataInt, lengthmz);

	//Remove Duplicates
	lengthmz = remove_duplicates(dataMZ, dataInt, lengthmz);
	//printf("New Length %d\n", lengthmz);

	//Background subtraction
	if (config.bsub != 0)
	{
		background_subtract(dataInt, config.bsub, lengthmz);
	}

	if (config.datareduction != 0)
	{
		//printf("Data Reduction: %f\n", config.datareduction);
		lengthmz = data_reduction(dataMZ, dataInt, lengthmz, config.datareduction);
	}

	//Normalize
	if (config.datanorm == 1)
	{
		//printf("Normalizing\n");
		norm1d(dataInt, lengthmz);
	}

	//Remove any weird negatives
	ApplyCutoff1D(dataInt, 0, lengthmz);

	config.speedyflag = 0;
	//Write data to processed_data
	strjoin(dataset, "/processed_data", outdat);
	//printf("\tWriting to: %s...", outdat);
	mh5writefile2d(file_id, outdat, lengthmz, dataMZ, dataInt);
	free(dataMZ);
	free(dataInt);
	H5Fclose(file_id);
	//printf("Done\n");
}

#endif