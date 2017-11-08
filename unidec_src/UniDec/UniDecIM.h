/*
* UniDecIM.h
*
*  Created on : 23 December 2016
* Author : Michael.Marty
*/

//
// Copyright 2015 University of Oxford
// Copyright 2016 University of Arizona
//
//

#if HAVE_MALLOC_H
#include <malloc.h>
#else
#include <stdlib.h>
#endif


#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>
//#include <omp.h>

void readfilemanual(char *infile, int lengthmz, double *array1, double *array2, double *array3, double *array4, double *array5)
{
	FILE *file_ptr;
	int i;
	char a[500];
	char b[500];
	char c[500];
	char d[500];
	char e[500];
	file_ptr = fopen(infile, "r");
	if (file_ptr == 0)
	{
		printf("Could not open %s file\n", infile);
		exit(10);
	}
	else {

		for (i = 0; i<lengthmz; i++)
		{
			fscanf(file_ptr, "%s %s %s %s %s", a, b, c, d, e);
			array1[i] = atof(a);
			array2[i] = atof(b);
			array3[i] = atof(c);
			array4[i] = atof(d);
			array5[i] = atof(e);
		}
	}
	fclose(file_ptr);
}

int GetSize0(double *array, int lines)
{
	int cnt = 1;
	double test = array[0];
	int i;
	for (i = 1; i<lines; i++) //look for the number of times dH decrements
	{
		if (array[i] != test)
		{
			test = array[i];
			cnt++;
		}
	}
	return cnt;
}

int GetSize1(double *array, int lines)
{
	int cnt = 1;
	int flag = 0;
	double test = array[0];
	int i;
	for (i = 1; i<lines; i++)
	{
		if (array[i] != test && flag == 0)
		{
			cnt++;
		}
		else
		{
			flag = 1;
		}
	}
	return cnt;
}

void Extract(double *mzext, double *dtext, double *mzdat, double *dtdat, int *size)
{
	int i;
	//Extract mz data
	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i<size[0]; i++)
	{
		mzext[i] = mzdat[i*size[1]];
	}
	//Extract dt data
	#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i<size[1]; i++)
	{
		dtext[i] = dtdat[i];
	}
	return;
}

double Gaus(double x, double x0, double sig)
{
	return exp(-pow(x - x0, 2.) / (2.*sig*sig));
}


//function to return a gaussian at a specific point
double Lorentz(double x, double x0, double sig)
{
	return pow(sig / 2, 2) / (pow(x - x0, 2) + pow(sig / 2, 2));
}

double SplitGL(double x, double y, double sig)
{
	if (y<x)
	{
		return exp(-(pow(x - y, 2)) / (2 * sig*sig*0.180337));
	}
	else
	{
		return (sig / 2)*(sig / 2) / (pow((x - y), 2) + pow((sig / 2), 2));
	}
}

//function to return a 2D peak value at a given point
double Peak(double x1, double x2, double y1, double y2, double sig1, double sig2, int psfun)
{
	if (psfun == 0) {
		return Gaus(x1, x2, sig1)*Gaus(y1, y2, sig2);
	}
	if (psfun == 1) {
		return Lorentz(x1, x2, sig1)*Gaus(y1, y2, sig2);
	}
	if (psfun == 2) {
		return SplitGL(x1, x2, sig1)*Gaus(y1, y2, sig2);
	}
	return 0;
}


static int index4D(int *size, int r, int c, int d, int f)
{
	return r*size[1] * size[2] * size[3] + c*size[2] * size[3] + d*size[3] + f;
}

void GetPeaks(double *peak, int *Size, double *dH, double *dC, double sig1, double sig2, int psfun)
{
	int i, j;
	double hmax = 2 * dH[Size[0] - 1] - dH[Size[0] - 2];
	double cmax = 2 * dC[Size[1] - 1] - dC[Size[1] - 2];
	double hmin = dH[0];
	double cmin = dC[0];
	double Hval, Cval, v1, v2, v3, v4, val;
	#pragma omp parallel for private (i,j,Hval,Cval,v1,v2,v3,v4,val), schedule(dynamic)
	for (i = 0; i<Size[0]; i++)
	{
		for (j = 0; j<Size[1]; j++)
		{
			Hval = dH[i];
			Cval = dC[j];

			//put the peakj in the four corners
			v1 = Peak(hmin, Hval, cmin, Cval, sig1, sig2, psfun);
			v2 = Peak(hmax, Hval, cmax, Cval, sig1, sig2, psfun);
			v3 = Peak(hmin, Hval, cmax, Cval, sig1, sig2, psfun);
			v4 = Peak(hmax, Hval, cmin, Cval, sig1, sig2, psfun);

			val = v1 + v2 + v3 + v4;
			peak[index2D(Size[1], i, j)] = val; //store peak in complete peak list
		}
	}
	//printf("mzmax: %f\n",hmax);
	//printf("dtmax: %f\n",cmax);
	return;
}

void fftconvolve2D(double *corr, double *list1, double *list2, int *size)
{
	fftw_complex *in1, *in2, *out1, *out2;
	fftw_plan p1, p2, p3;
	double a, b, c, d;
	int i;
	int length = size[0] * size[1];
	in1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
	in2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
	out1 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
	out2 = (fftw_complex*)fftw_malloc(sizeof(fftw_complex) * length);
#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i<length; i++)
	{
		in1[i][0] = list1[i];
		in1[i][1] = 0;
		in2[i][0] = list2[i];
		in2[i][1] = 0;
	}
	p1 = fftw_plan_dft_2d(size[0], size[1], in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftw_plan_dft_2d(size[0], size[1], in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftw_execute(p1);
	fftw_execute(p2);
#pragma omp parallel for private (i,a,b,c,d), schedule(dynamic)
	for (i = 0; i<length; i++)
	{
		a = out1[i][0];
		b = out1[i][1];
		c = out2[i][0];
		d = out2[i][1];
		in1[i][0] = a*c - b*d;
		in1[i][1] = b*c + a*d;
	}

	p3 = fftw_plan_dft_2d(size[0], size[1], in1, out1, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftw_execute(p3);
#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i<length; i++)
	{
		corr[i] = out1[i][0] / length;
	}
	//
	fftw_destroy_plan(p1);
	fftw_destroy_plan(p2);
	fftw_destroy_plan(p3);
	fftw_free(in1); fftw_free(in2); fftw_free(out1); fftw_free(out2);
}

void convolve3D(double *out, double *in, double *peakshape, int *size)
{
	//Initialize memory
	double *temp = NULL, *tempout = NULL;
	temp = calloc(size[0] * size[1], sizeof(double));
	tempout = calloc(size[0] * size[1], sizeof(double));
	int i, j, k;
	//Convolve for each charge slice
	for (k = 0; k<size[2]; k++)
	{
		//Cut out a slice
#pragma omp parallel for private (i,j), schedule(dynamic)
		for (i = 0; i<size[0]; i++)
		{
			for (j = 0; j<size[1]; j++)
			{
				temp[index2D(size[1], i, j)] = in[index3D(size[1], size[2], i, j, k)];
			}
		}
		//Convolve it
		fftconvolve2D(tempout, temp, peakshape, size);
		//Repack slice into output array
#pragma omp parallel for private (i,j), schedule(dynamic)
		for (i = 0; i<size[0]; i++)
		{
			for (j = 0; j<size[1]; j++)
			{
				out[index3D(size[1], size[2], i, j, k)] = tempout[index2D(size[1], i, j)];
			}
		}
	}
	//Free Memory
	free(temp);
	free(tempout);
}

double calcCCS(double mass, int z, double dt, double ccsconst, double hmass, double to)
{
	double ac = 1.6605389E-27;
	double rmass = ac*(mass*hmass) / (mass + hmass);
	double td = (dt - to)*0.001;
	double ccs = (double)z*(td)*sqrt(1 / rmass)*ccsconst;
	return ccs;
}

double calcCCSTwaveLog(double mass, int z, double dt, double tcal1, double tcal2, double hmass, double edc)
{
	double rdt = dt - (edc*sqrt(mass / z) / 1.000E3);
	double rccs;
	if (rdt <= 0) { return 0; }
	else { rccs = exp(tcal1*log(rdt) + tcal2); }
	double rmass = (mass*hmass) / (mass + hmass);
	return rccs*z*sqrt(1 / rmass);
}

double calcCCSTwaveLinear(double mass, int z, double dt, double tcal1, double tcal2, double hmass, double edc)
{
	double rdt = dt - (edc*sqrt(mass / z) / 1.000E3);
	double rccs;
	if (rdt <= 0) { return 0; }
	else { rccs = tcal1*rdt + tcal2; }
	double rmass = (mass*hmass) / (mass + hmass);
	return rccs*z*sqrt(1 / rmass);
}

double calcCCSTwavePower(double mass, int z, double dt, double tcal1, double tcal2, double hmass, double edc)
{
	double rdt = dt - (edc*sqrt(mass / z) / 1.000E3);
	double rccs;
	if (rdt <= 0) { return 0; }
	else { rccs = tcal1*pow(rdt, tcal2); }
	double rmass = (mass*hmass) / (mass + hmass);
	return rccs*z*sqrt(1 / rmass);
}

double calcDtTwaveLog(double mass, int z, double ccs, double tcal1, double tcal2, double hmass, double edc)
{
	//double ac=1.6605389E-27;
	double rmass = (mass*hmass) / (mass + hmass);
	double rccs = ccs / (z*sqrt(1 / rmass));
	double rdt;
	if (rccs <= 0) { return 0; }
	else { rdt = exp((log(rccs) - tcal2) / tcal1); }
	return rdt + (edc*sqrt(mass / z) / 1000);
}

double calcDtTwaveLinear(double mass, int z, double ccs, double tcal1, double tcal2, double hmass, double edc)
{
	double rmass = (mass*hmass) / (mass + hmass);
	double rccs = ccs / (z*sqrt(1 / rmass));
	double rdt;
	if (rccs <= 0) { return 0; }
	else { rdt = (rccs - tcal2) / tcal1; }
	return rdt + (edc*sqrt(mass / z) / 1000);
}

double calcDtTwavePower(double mass, int z, double ccs, double tcal1, double tcal2, double hmass, double edc)
{
	double rmass = (mass*hmass) / (mass + hmass);
	double rccs = ccs / (z*sqrt(1 / rmass));
	double rdt;
	if (rccs <= 0) { return 0; }
	else { rdt = pow(rccs / tcal1, 1.0 / tcal2); }
	return rdt + (edc*sqrt(mass / z) / 1000);
}

double calcDt(double mass, int z, double ccs, double ccsconst, double hmass, double to)
{

	double ac = 1.6605389E-27;
	double rmass = ac*(mass*hmass) / (mass + hmass);
	double td = ccs / (ccsconst*(double)z*sqrt(1 / rmass));
	double dt = (td * 1000) + to;
	return dt;
}

double calcMass(double mz, int z, double adductmass)
{
	return mz*(double)z - adductmass*(double)z;
}

double calcMz(double mass, int z, double adductmass)
{
	return (mass + z*adductmass) / (double)z;
}

int nearint(int *array, int point, int numdat)
{
	int i;
	int temp = -1;
	for (i = 0; i<numdat; i++)
	{
		if (array[i] == point) { temp = i; }
	}
	if (temp == -1) { printf("Problem finding integer in list\n"); }
	return temp;
}

int nearfastpoint(double *dataMZ, double point, int numdat)
{
	int start = 0;
	int length = numdat - 1;
	int end = 0;

	int diff = length - start;
	while (diff>1)
	{
		if (point<dataMZ[start + (length - start) / 2])
		{
			length = start + (length - start) / 2;
		}
		else if (point == dataMZ[start + (length - start) / 2])
		{
			end = start + (length - start) / 2;
			length = start + (length - start) / 2;
			start = start + (length - start) / 2;
			return end;
		}
		else if (point>dataMZ[start + (length - start) / 2])
		{
			start = start + (length - start) / 2;
		}
		diff = length - start;
	}
	if (fabs(point - dataMZ[start]) >= fabs(point - dataMZ[length]))
	{
		end = length;
	}
	else
	{
		end = start;
	}
	return dataMZ[end];
}

double euclid(double a, double b, double c, double d)
{
	return sqrt((a - b)*(a - b) + (c - d)*(c - d));
}

//Slow nearest unsorted.
int nearunsorted_IM(double *testmasses, double *testtimes, double mz, double dt, int lengthtest)
{
	double minval = euclid(mz, testmasses[0], dt, testtimes[0]);
	double difftest;
	int pos = 0;
	for (int i = 1; i<lengthtest; i++)
	{
		difftest = euclid(mz, testmasses[i], dt, testtimes[i]);
		if (difftest<minval)
		{
			minval = difftest;
			pos = i;
		}
	}
	return pos;
}

void sum2D(int *size, double *sum, double *array, int axis)
{
	int i, j, k;
	double temp;
	if (axis == 2)
	{
		for (i = 0; i<size[0]; i++)
		{
			for (j = 0; j<size[1]; j++)
			{
				temp = 0;
				for (k = 0; k<size[2]; k++)
				{
					temp += array[index3D(size[1], size[2], i, j, k)];
				}
				sum[index2D(size[1], i, j)] = temp;
			}
		}
	}
	else if (axis == 1)
	{
		for (i = 0; i<size[0]; i++)
		{
			for (k = 0; k<size[2]; k++)
			{
				temp = 0;
				for (j = 0; j<size[1]; j++)
				{
					temp += array[index3D(size[1], size[2], i, j, k)];
				}
				sum[index2D(size[2], i, k)] = temp;
			}
		}
	}
	else if (axis == 0)
	{
		for (k = 0; k<size[2]; k++)
		{
			for (j = 0; j<size[1]; j++)
			{
				temp = 0;
				for (i = 0; i<size[0]; i++)
				{
					temp += array[index3D(size[1], size[2], i, j, k)];
				}
				sum[index2D(size[1], k, j)] = temp;
			}
		}
	}
	return;
}

void blur_it_IM(int *size, double *blur, double *newblur, int *closetab, int *barr, double csig)
{

	if (csig < 0) {
		double *mzzsumgrid = NULL;
		mzzsumgrid = calloc(size[0] * size[2], sizeof(double));
		sum2D(size, mzzsumgrid, blur, 1);
#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < size[0]; i++)
		{
			for (int j = 0; j < size[1]; j++)
			{
				for (int k = 0; k < size[2]; k++)
				{
					if (barr[index3D(size[1], size[2], i, j, k)] == 1)
					{
						blur[index3D(size[1], size[2], i, j, k)] = mzzsumgrid[index2D(size[2], i, k)];
					}
					else {
						blur[index3D(size[1], size[2], i, j, k)] = 0;
					}
				}
			}
		}
		free(mzzsumgrid);


#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i<size[0]; i++)
		{
			for (int j = 0; j<size[1]; j++)
			{
				for (int k = 0; k<size[2]; k++)
				{
					if (barr[index3D(size[1], size[2], i, j, k)] == 1)
					{
						double junk2 = 0;
						for (int l = 0; l<size[3]; l++)
						{
							int index = closetab[index4D(size, i, j, k, l)];
							if (index == -1) {
								junk2 += -12;
							}
							else
							{
								double junk = blur[index];
								if (junk>0) { junk2 += log(junk); }
								else { junk2 += -12; }
							}
						}
						newblur[index3D(size[1], size[2], i, j, k)] = exp(junk2 / size[3]);
					}
					else
					{
						newblur[index3D(size[1], size[2], i, j, k)] = 0;
					}
				}
			}
		}
	}

	else {
		double junk, junk2;
		int i, j, k, l, index;
		int logswitch = 1;
#pragma omp parallel for private (i,j,k,l,index,junk,junk2), schedule(dynamic)
		for (i = 0; i < size[0]; i++)
		{
			for (j = 0; j < size[1]; j++)
			{
				for (k = 0; k < size[2]; k++)
				{
					if (barr[index3D(size[1], size[2], i, j, k)] == 1)
					{
						junk2 = 0;
						for (l = 0; l < size[3]; l++)
						{
							index = closetab[index4D(size, i, j, k, l)];
							if (index == -1) {
								if (logswitch == 1) { junk2 += -12; }
								else { junk2 += 0; }
							}
							else
							{
								junk = blur[index];
								if (logswitch == 0) { junk2 += junk; }
								else {
									if (junk > 0) { junk2 += log(junk); }
									else { junk2 += -12; }
								}
							}
						}
						if (logswitch == 0) { newblur[index3D(size[1], size[2], i, j, k)] = junk2 / size[3]; }
						else { newblur[index3D(size[1], size[2], i, j, k)] = exp(junk2 / size[3]); }
					}
					else
					{
						newblur[index3D(size[1], size[2], i, j, k)] = 0;
					}
				}
			}
		}
	}
}

void sumdeltas(int *size, double *deltas, double *newblur)
{
	int i, j, k;
	double temp;
#pragma omp parallel for private (i,j,temp), schedule(dynamic)
	for (i = 0; i<size[0]; i++)
	{
		for (j = 0; j<size[1]; j++)
		{
			temp = 0;
			for (k = 0; k<size[2]; k++)
			{
				temp += newblur[index3D(size[1], size[2], i, j, k)];
			}
			deltas[index2D(size[1], i, j)] = temp;
		}
	}

}

void sum1D(int *size, double *sum, double *array, int axis)
{
	int i, j, k;
	double temp;
	if (axis == 0) {
		for (i = 0; i<size[0]; i++)
		{
			temp = 0;
			for (j = 0; j<size[1]; j++)
			{

				for (k = 0; k<size[2]; k++)
				{
					temp += array[index3D(size[1], size[2], i, j, k)];
				}

			}
			sum[i] = temp;
		}
	}
	else if (axis == 1)
	{
		for (j = 0; j<size[1]; j++)
		{
			temp = 0;
			for (i = 0; i<size[0]; i++)
			{

				for (k = 0; k<size[2]; k++)
				{
					temp += array[index3D(size[1], size[2], i, j, k)];
				}

			}
			sum[j] = temp;
		}
	}
	else if (axis == 2)
	{
		for (k = 0; k<size[2]; k++)
		{
			temp = 0;
			for (j = 0; j<size[1]; j++)
			{

				for (i = 0; i<size[0]; i++)
				{
					temp += array[index3D(size[1], size[2], i, j, k)];
				}

			}
			sum[k] = temp;
		}
	}
	return;
}



void bayes(int *size, double*denom, double *blur, double *newblur, double *dataInt, int *barr)
{
	int i, j, k;
#pragma omp parallel for private (i,j,k), schedule(dynamic)
	for (i = 0; i<size[0]; i++)
	{
		for (j = 0; j<size[1]; j++)
		{
			for (k = 0; k<size[2]; k++)
			{
				if (denom[index2D(size[1], i, j)] != 0 && barr[index3D(size[1], size[2], i, j, k)] == 1)
				{
					blur[index3D(size[1], size[2], i, j, k)] = newblur[index3D(size[1], size[2], i, j, k)] * dataInt[index2D(size[1], i, j)] / denom[index2D(size[1], i, j)];
				}
				else
				{
					blur[index3D(size[1], size[2], i, j, k)] = 0;
				}
			}
		}
	}
	return;
}

void getranges(int *size, double *blur, double *masstab, double *ccstab, double *range, int *barr)
{
	int i, j, k;
	int mass, ccs;
	range[0] = 100000000;
	range[1] = 0;
	range[2] = 100000000;
	range[3] = 0;
	for (i = 0; i<size[0]; i++)
	{
		for (j = 0; j<size[1]; j++)
		{
			for (k = 0; k<size[2]; k++)
			{
				if (blur[index3D(size[1], size[2], i, j, k)]>1E-12&&barr[index3D(size[1], size[2], i, j, k)] == 1)
				{
					mass = masstab[index2D(size[2], i, k)];
					ccs = ccstab[index3D(size[1], size[2], i, j, k)];
					if (mass<range[0]) { range[0] = mass; }
					if (mass>range[1]) { range[1] = mass; }
					if (ccs<range[2]) { range[2] = ccs; }
					if (ccs>range[3]) { range[3] = ccs; }
				}
			}
		}
	}
	return;
}

void makeaxis(double *axis, int length, double min, double binsize)
{
	int i;
#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i<length; i++)
	{
		axis[i] = min + i*binsize;
	}
	return;
}

void write1D(char *outfile, char *suffix, double *array, int length)
{
	char outstring[500];
	int i;
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s.bin", outfile, suffix);
	out_ptr = fopen(outstring, "wb");
	fwrite(array, sizeof(double), length, out_ptr);

	//sprintf(outstring,"%s_%s.txt",outfile,suffix);
	//out_ptr=fopen(outstring,"w");
	//for(i=0;i<length;i++)
	//{
	//	  fprintf(out_ptr,"%f\n",array[i]);
	//}
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}

void write2D(char *outfile, char *suffix, double *array1, double *array2, int length)
{
	char outstring[500];
	int i;
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s.txt", outfile, suffix);
	out_ptr = fopen(outstring, "w");
	for (i = 0; i<length; i++)
	{
		fprintf(out_ptr, "%f %f\n", array1[i], array2[i]);
	}
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}

void write3D(char *outfile, char *suffix, double *array1, double *array2, double *array3, int length)
{
	char outstring[500];
	int i;
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s.txt", outfile, suffix);
	out_ptr = fopen(outstring, "w");
	for (i = 0; i<length; i++)
	{
		fprintf(out_ptr, "%f %f %f\n", array1[i], array2[i], array3[i]);
	}
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}

void writemfileres(char *outfile, char *suffix, double *array1, double *array2, double *array3, double *array4, double *array5, double *array6, int length)
{
	char outstring[500];
	int i;
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s.txt", outfile, suffix);
	out_ptr = fopen(outstring, "w");
	//fprintf(out_ptr,"Mass\tIntensity\tCCS Avg.\tCCS Std. Dev.\tZ avg.\tZ Std. Dev.\n");
	for (i = 0; i<length; i++)
	{
		fprintf(out_ptr, "%f\t%f\t%f\t%f\t%f\t%f\n", array1[i], array2[i], array3[i], array4[i], array5[i], array6[i]);
	}
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}

void writemzgrid(char *outfile, char *suffix, double *mzext, double *dtext, int *ztab, double *blur, int *size)
{
	char outstring[500];
	//int i,j,k;
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s.bin", outfile, suffix);
	//printf("%s\n",outstring);
	out_ptr = fopen(outstring, "wb");
	/*
	for(i=0;i<size[0];i++)
	{
	for(j=0;j<size[1];j++)
	{
	for(k=0;k<size[2];k++)
	{
	fprintf(out_ptr,"%f %f %d %f\n",mzext[i],dtext[j],ztab[k],blur[index3D(size[1],size[2],i,j,k)]);
	}
	}
	}
	*/
	fwrite(blur, sizeof(double), size[0] * size[1] * size[2], out_ptr);
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}

double errfun(int length, double *dataInt, double *fitdat)
{
	int i;
	double tot = 0;
	double mean;

	for (i = 0; i<length; i++)
	{
		tot += fitdat[i];
	}
	mean = tot / ((double)length);

	double error = 0;
	double error2 = 0;
	for (i = 0; i<length; i++)
	{
		error += pow((fitdat[i] - dataInt[i]), 2);
		error2 += pow((mean - dataInt[i]), 2);
	}
	if (error2 == 0) { printf("Error function divide by 0!"); error2 = 1; }
	return 1.0 - error / error2;
}

double getmax(int length, double *data)
{
	double max = 0;
	for (int i = 0; i<length; i++)
	{
		if (data[i]>max) { max = data[i]; }
	}
	return max;
}

void normalize(int length, double *data, double norm)
{
	if (norm == 0) { printf("Tried to Normalize by Zero! Aborting!\n"); exit(50); }
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<length; i++)
	{
		data[i] = data[i] / norm;
	}
	return;
}

void writezslice(int *size, char *outfile, char *suffix, double *massaxis, double *ccsaxis, int *ztab, double *array, int k)
{
	int i, j;
	char outstring[500];
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s_%d.bin", outfile, suffix, ztab[k]);
	out_ptr = fopen(outstring, "wb");
	double *temp = NULL;
	temp = calloc(size[0] * size[1], sizeof(double));
	for (i = 0; i<size[0]; i++)
	{
		for (j = 0; j<size[1]; j++)
		{
			//fprintf(out_ptr,"%f %f %f\n",massaxis[i],ccsaxis[j],array[index3D(size[1],size[2],i,j,k)]);
			temp[index2D(size[1], i, j)] = array[index3D(size[1], size[2], i, j, k)];
		}
	}
	fwrite(temp, sizeof(double), size[0] * size[1], out_ptr);
	fclose(out_ptr);
	free(temp);
	printf("File written to: %s\n", outstring);
}


double nativeCCS(double mass, double fudge, double gasmass)
{
	double a, b;
	if (gasmass<10)
	{
		a = 4.06739;
		b = 0.629424;
	}
	else {
		a = 5.33311;
		b = 0.613072;
	}
	//double c=-1481.7;
	return a*pow(mass, b) + fudge;
}

double bilinearinterpolation(int *size, int indm, int indc, int k, double *mzext, double *dtext, double tempmz, double tempdt, int rawflag, double *newblur, double *blur)
{
	double y11, y12, y21, y22, mu1, mu2;
	int im1, im2, ic1, ic2;
	if (mzext[indm]>tempmz) { im2 = indm; im1 = indm - 1; }
	else { im1 = indm; im2 = indm + 1; }
	if (dtext[indc]>tempdt) { ic2 = indc; ic1 = indc - 1; }
	else { ic1 = indc; ic2 = indc + 1; }
	mu1 = (tempmz - mzext[im1]) / (mzext[im2] - mzext[im1]);
	mu2 = (tempdt - dtext[ic1]) / (dtext[ic2] - dtext[ic1]);
	if (rawflag == 0) {
		y11 = newblur[index3D(size[1], size[2], im1, ic1, k)];
		y12 = newblur[index3D(size[1], size[2], im1, ic2, k)];
		y21 = newblur[index3D(size[1], size[2], im2, ic1, k)];
		y22 = newblur[index3D(size[1], size[2], im2, ic2, k)];
	}
	else {
		y11 = blur[index3D(size[1], size[2], im1, ic1, k)];
		y12 = blur[index3D(size[1], size[2], im1, ic2, k)];
		y21 = blur[index3D(size[1], size[2], im2, ic1, k)];
		y22 = blur[index3D(size[1], size[2], im2, ic2, k)];
	}
	return y11*(1 - mu1)*(1 - mu2) + y12*mu1*(1 - mu2) + y21*mu2*(1 - mu1) + y22*mu1*mu2;
}

//cubicInterpolate_IM and bicubicInterpolate_IM were taken from http://www.paulinternet.nl/?page=bicubic
double cubicInterpolate_IM(double p[4], double x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}
double bicubicInterpolate_IM(double p[4][4], double x, double y) {
	double arr[4];
	arr[0] = cubicInterpolate_IM(p[0], y);
	arr[1] = cubicInterpolate_IM(p[1], y);
	arr[2] = cubicInterpolate_IM(p[2], y);
	arr[3] = cubicInterpolate_IM(p[3], y);
	return cubicInterpolate_IM(arr, x);
}

double bicubicinterpolation(int *size, int indm, int indc, int k, double *mzext, double *dtext, double tempmz, double tempdt, int rawflag, double *newblur, double *blur)
{
	double mu1, mu2;
	double p[4][4];
	int im[4], ic[4];
	int i, j;
	if (mzext[indm]>tempmz) { im[2] = indm; im[1] = indm - 1; im[0] = indm - 2; im[3] = indm + 1; }
	else { im[1] = indm; im[2] = indm + 1; im[0] = indm - 1; im[3] = indm + 2; }
	if (dtext[indc]>tempdt) { ic[2] = indc; ic[1] = indc - 1; ic[0] = indc - 2; ic[3] = indc + 1; }
	else { ic[1] = indc; ic[2] = indc + 1; ic[0] = indc - 1; ic[3] = indc + 2; }
	mu1 = (tempmz - mzext[im[1]]) / (mzext[im[2]] - mzext[im[1]]);
	mu2 = (tempdt - dtext[ic[1]]) / (dtext[ic[2]] - dtext[ic[1]]);
	for (i = 0; i<4; i++)
	{
		for (j = 0; j<4; j++)
		{
			if (rawflag == 0) { p[i][j] = newblur[index3D(size[1], size[2], im[i], ic[j], k)]; }
			else { p[i][j] = blur[index3D(size[1], size[2], im[i], ic[j], k)]; }
		}
	}
	return clip(bicubicInterpolate_IM(p, mu1, mu2), 0);
}

double WeightedAverage(int length, double *xarray, double *warray)
{
	double temp1 = 0;
	double temp2 = 0;
	int i;
	for (i = 0; i<length; i++)
	{
		temp1 += xarray[i] * warray[i];
		temp2 += warray[i];
	}
	return temp1 / temp2;
}

double WeightedAverageInt(int length, int *xarray, double *warray)
{
	double temp1 = 0;
	double temp2 = 0;
	int i;
	for (i = 0; i<length; i++)
	{
		temp1 += xarray[i] * warray[i];
		temp2 += warray[i];
	}
	return temp1 / temp2;
}

double StdDev2(int length, double *xarray, double *warray, double wmean)
{
	double temp1 = 0;
	double temp2 = 0;
	int i;
	for (i = 0; i<length; i++)
	{
		temp1 += pow(xarray[i] - wmean, 2)*warray[i];
		temp2 += warray[i];
	}
	return sqrt(temp1 / temp2);
}

double StdDev2Int(int length, int *xarray, double *warray, double wmean)
{
	double temp1 = 0;
	double temp2 = 0;
	int i;
	for (i = 0; i<length; i++)
	{
		temp1 += pow(xarray[i] - wmean, 2)*warray[i];
		temp2 += warray[i];
	}
	return sqrt(temp1 / temp2);
}


void MFileInt(int mfilelen, double *massaxis, double *massaxisval, double *testmasses, double *testmassint, double maaxle)
{
	int i, index;
	for (i = 0; i<mfilelen; i++)
	{
		index = nearfast(massaxis, testmasses[i], maaxle);
		testmassint[i] = massaxisval[index];
	}
}

//Works for both CCS and Z to get weighted average and Standard Deviation.
void MFileCCS(int mfilelen, double *massaxis, double *ccsaxis, double *massccsgrid, double *testmasses, double *testmassCCSavg, double *testmassCCSstddev, int maaxle, int ccaxle)
{
	int i, j, index;
	double *temp;
	double wmean;
	temp = calloc(ccaxle, sizeof(double));
	for (i = 0; i<mfilelen; i++)
	{
		index = nearfast(massaxis, testmasses[i], maaxle);
		for (j = 0; j<ccaxle; j++)
		{
			temp[j] = massccsgrid[index2D(ccaxle, index, j)];
		}
		wmean = WeightedAverage(ccaxle, ccsaxis, temp);
		testmassCCSavg[i] = wmean;
		testmassCCSstddev[i] = StdDev2(ccaxle, ccsaxis, temp, wmean);
	}
	free(temp);
}

void MFileZ(int mfilelen, double *massaxis, int *ztab, double *masszgrid, double *testmasses, double *testmassZavg, double *testmassZstddev, int maaxle, int numz)
{
	int i, j, index;
	double *temp;
	double wmean;
	temp = calloc(numz, sizeof(double));
	for (i = 0; i<mfilelen; i++)
	{
		index = nearfast(massaxis, testmasses[i], maaxle);
		for (j = 0; j<numz; j++)
		{
			temp[j] = masszgrid[index2D(numz, index, j)];
		}
		wmean = WeightedAverageInt(numz, ztab, temp);
		testmassZavg[i] = wmean;
		testmassZstddev[i] = StdDev2Int(numz, ztab, temp, wmean);
	}
	free(temp);
}

void KillB_IM(double *I, int *B, int * size, double intthresh)
{
	for (int i = 0; i<size[0]; i++)
	{
		for (int j = 0; j<size[1]; j++)
		{
			for (int k = 0; k<size[2]; k++)
			{
				if (I[index2D(size[1], i, j)] < intthresh) {
					B[index3D(size[1], size[2], i, j, k)] = 0;
				}
			}
		}
	}
}



void ManualAssign_IM(char *manualfile, int * size, double *mzdat, double * dtdat, int *ztab, int *barr)
{

	//Manual Assignments
	int manlength;
	double *manmz = NULL,
		*manmzwin = NULL,
		*mandt = NULL,
		*mandtwin = NULL,
		*mancharge = NULL;
	manlength = getfilelength(manualfile);
	manmz = calloc(manlength, sizeof(double));
	manmzwin = calloc(manlength, sizeof(double));
	mandt = calloc(manlength, sizeof(double));
	mandtwin = calloc(manlength, sizeof(double));
	mancharge = calloc(manlength, sizeof(double));
	readfilemanual(manualfile, manlength, manmz, manmzwin, mandt, mandtwin, mancharge);
	printf("Read manual file of length: %d\n", manlength);

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<size[0]; i++)
	{
		for (int j = 0; j<size[1]; j++)
		{
			double mz = mzdat[index2D(size[1], i, j)];
			double dt = dtdat[index2D(size[1], i, j)];
			int close = nearunsorted_IM(manmz, mandt, mz, dt, manlength);
			int charge = mancharge[close];
			if (fabs(manmz[close] - mz) < manmzwin[close] && fabs(mandt[close] - dt) < mandtwin[close])
			{
				for (int k = 0; k < size[2]; k++)
				{
					if (ztab[k] != charge)
					{
						barr[index3D(size[1], size[2], i, j, k)] = 0;
					}
				}
			}
		}
	}

	free(manmz);
	free(manmzwin);
	free(mandt);
	free(mandtwin);
	free(mancharge);
	printf("Using Manual Assignments for Some Peaks\n");
}

void TwaveError(int twaveflag)
{
	printf("Error: Undefined twaveflag. Value was %d. Aborting.", twaveflag);
	exit(20);
}
