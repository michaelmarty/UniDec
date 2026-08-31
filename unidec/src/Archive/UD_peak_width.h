/*
* MetaUniDec_Main.h
*
*  Created on : 12 June 2017
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//


int autocorr_and_stop(const float * dataInt, const int lengthmz)
{
	float val = 0;
	int index = 0;
	//get pak of autocorrelation
	for (int j = 0; j < lengthmz; j++)
	{
		val += dataInt[j] * dataInt[j];
	}
	float initval = val;

	//find local minimum
	for (int i = 1; i < lengthmz; i++)
	{
		float sum = 0;
		for (int j = 0; j < lengthmz; j++)
		{
			sum += dataInt[j] * dataInt[(i + j) % lengthmz];//CHECK THAT THIS DIRECTION IS CORRECT MIGHT NEED INVERSE
		}
		if (sum < val) { val = sum; }
		else {
			val = sum;
			index = i;
			break;
		}
	}
	printf("Min: %d %f %f", index, initval, val);
	return index;
}

void local_peak_widths(const float * dataMZ, const float * dataInt, const int lengthmz)
{
	float *dmin = NULL;
	float *dmax = NULL;
	dmin = calloc(lengthmz, sizeof(float));
	dmax = calloc(lengthmz, sizeof(float));

	for (int i = 0; i < lengthmz; i++)
	{
		float xval = dataMZ[i];
		float yval = dataInt[i];
		float minval = yval*yval;
		float maxval = 0;
		int minindex = 0;
		int maxindex = 0;
		float minx = 0;
		float maxx = 0;
		for (int j = 1; j < lengthmz; j++)
		{
			int mini = i - j;
			int maxi = i + j;
			if (mini < 0) { mini = 0; }// mini = -mini;}
			if (maxi >= lengthmz) { maxi = lengthmz - 1; }//{ maxi = 2 * lengthmz - maxi-1; }
			float val = dataInt[mini] * dataInt[maxi];
			if (val < minval) { minval = val; }
			else if (maxval == 0)
			{
				minindex = j;
				minx = dataMZ[j] - dataMZ[i];
				maxval = minval;
			}
			else if (maxval != 0 && val > maxval) { maxval = val; }
			else if (maxval != 0 && val <= maxval)
			{
				maxindex = j;
				maxx = dataMZ[j] - dataMZ[i];
				break;
			}
		}

		dmin[i] = minval;
		dmax[i] = maxval;
	}



	free(dmin);
	free(dmax);

}


void get_peak_widths(int argc, char *argv[], Config config)
{
	char dataset[1024];
	char outdat[1024];
	char strval[1024];

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

	int index = autocorr_and_stop(masssum, mlen);

	free(massaxis);
	free(masssum);
}

/*
The goal here is to determine the local peak width and feed that into the deconvolution.
From the local peak width, you might be able to get whether it was noise (peak width is on the order of data density)
or whether it was baseline (peak width much larger than the standard autocorrelation).
It may be possible from the auto peak width to determine whether it is a peak or not,
which would readily allow baseline calculation and dimensional reduction.

The local peak width can be in some way determined by the minimum of the local autocorrelation. 
For the secondary maximum of the local autocorrelation, it would reveal the potential isotopic spacing and thus the charge. 
If the isotopic spacing could be determined, it could lead to a much smarter isotope mode that ignored noise and any non-isotopic peaks.

*/