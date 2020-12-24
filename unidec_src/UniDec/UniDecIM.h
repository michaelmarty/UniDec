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

void readfilemanual(char *infile, int lengthmz, float *array1, float *array2, float *array3, float *array4, float *array5)
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

int GetSize0(float *array, int lines)
{
	int cnt = 1;
	float test = array[0];
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

int GetSize1(float *array, int lines)
{
	int cnt = 1;
	int flag = 0;
	float test = array[0];
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

void Extract(float *mzext, float *dtext, float *mzdat, float *dtdat, int *size)
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

float Gaus(float x, float x0, float sig)
{
	return exp(-pow(x - x0, 2.) / (2.*sig*sig));
}


//function to return a gaussian at a specific point
float Lorentz(float x, float x0, float sig)
{
	return pow(sig / 2, 2) / (pow(x - x0, 2) + pow(sig / 2, 2));
}

float SplitGL(float x, float y, float sig)
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
float Peak(float x1, float x2, float y1, float y2, float sig1, float sig2, int psfun)
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

void GetPeaks(float *peak, int *Size, float *dH, float *dC, float sig1, float sig2, int psfun)
{
	int i, j;
	float hmax = 2 * dH[Size[0] - 1] - dH[Size[0] - 2];
	float cmax = 2 * dC[Size[1] - 1] - dC[Size[1] - 2];
	float hmin = dH[0];
	float cmin = dC[0];
	float Hval, Cval, v1, v2, v3, v4, val;
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

void fftconvolve2D(float *corr, float *list1, float *list2, int *size)
{
	fftw_complex *in1, *in2, *out1, *out2;
	fftw_plan p1, p2, p3;
	float a, b, c, d;
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

void convolve3D(float *out, float *in, float *peakshape, int *size)
{
	//Initialize memory
	float *temp = NULL, *tempout = NULL;
	temp = calloc(size[0] * size[1], sizeof(float));
	tempout = calloc(size[0] * size[1], sizeof(float));
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

float calcCCS(float mass, int z, float dt, float ccsconst, float hmass, float to)
{
	float ac = 1.6605389E-27;
	float rmass = ac*(mass*hmass) / (mass + hmass);
	float td = (dt - to)*0.001;
	float ccs = (float)z*(td)*sqrt(1 / rmass)*ccsconst;
	return ccs;
}

float calcCCSTwaveLog(float mass, int z, float dt, float tcal1, float tcal2, float hmass, float edc)
{
	float rdt = dt - (edc*sqrt(mass / z) / 1.000E3);
	float rccs;
	if (rdt <= 0) { return 0; }
	else { rccs = exp(tcal1*log(rdt) + tcal2); }
	float rmass = (mass*hmass) / (mass + hmass);
	return rccs*z*sqrt(1 / rmass);
}

float calcCCSTwaveLinear(float mass, int z, float dt, float tcal1, float tcal2, float hmass, float edc)
{
	float rdt = dt - (edc*sqrt(mass / z) / 1.000E3);
	float rccs;
	if (rdt <= 0) { return 0; }
	else { rccs = tcal1*rdt + tcal2; }
	float rmass = (mass*hmass) / (mass + hmass);
	return rccs*z*sqrt(1 / rmass);
}

float calcCCSTwavePower(float mass, int z, float dt, float tcal1, float tcal2, float hmass, float edc)
{
	float rdt = dt - (edc*sqrt(mass / z) / 1.000E3);
	float rccs;
	if (rdt <= 0) { return 0; }
	else { rccs = tcal1*pow(rdt, tcal2); }
	float rmass = (mass*hmass) / (mass + hmass);
	return rccs*z*sqrt(1 / rmass);
}

float calcDtTwaveLog(float mass, int z, float ccs, float tcal1, float tcal2, float hmass, float edc)
{
	//float ac=1.6605389E-27;
	float rmass = (mass*hmass) / (mass + hmass);
	float rccs = ccs / (z*sqrt(1 / rmass));
	float rdt;
	if (rccs <= 0) { return 0; }
	else { rdt = exp((log(rccs) - tcal2) / tcal1); }
	return rdt + (edc*sqrt(mass / z) / 1000);
}

float calcDtTwaveLinear(float mass, int z, float ccs, float tcal1, float tcal2, float hmass, float edc)
{
	float rmass = (mass*hmass) / (mass + hmass);
	float rccs = ccs / (z*sqrt(1 / rmass));
	float rdt;
	if (rccs <= 0) { return 0; }
	else { rdt = (rccs - tcal2) / tcal1; }
	return rdt + (edc*sqrt(mass / z) / 1000);
}

float calcDtTwavePower(float mass, int z, float ccs, float tcal1, float tcal2, float hmass, float edc)
{
	float rmass = (mass*hmass) / (mass + hmass);
	float rccs = ccs / (z*sqrt(1 / rmass));
	float rdt;
	if (rccs <= 0) { return 0; }
	else { rdt = pow(rccs / tcal1, 1.0 / tcal2); }
	return rdt + (edc*sqrt(mass / z) / 1000);
}

float calcDt(float mass, int z, float ccs, float ccsconst, float hmass, float to)
{

	float ac = 1.6605389E-27;
	float rmass = ac*(mass*hmass) / (mass + hmass);
	float td = ccs / (ccsconst*(float)z*sqrt(1 / rmass));
	float dt = (td * 1000) + to;
	return dt;
}

float calcMass(float mz, int z, float adductmass)
{
	return mz*(float)z - adductmass*(float)z;
}

float calcMz(float mass, int z, float adductmass)
{
	return (mass + z*adductmass) / (float)z;
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

int nearfastpoint(float *dataMZ, float point, int numdat)
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

float euclid(float a, float b, float c, float d)
{
	return sqrt((a - b)*(a - b) + (c - d)*(c - d));
}

//Slow nearest unsorted.
int nearunsorted_IM(float *testmasses, float *testtimes, float mz, float dt, int lengthtest)
{
	float minval = euclid(mz, testmasses[0], dt, testtimes[0]);
	float difftest;
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

void sum2D(int *size, float *sum, float *array, int axis)
{
	int i, j, k;
	float temp;
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

void blur_it_IM(int *size, float *blur, float *newblur, int *closetab, int *barr, float csig)
{

	if (csig < 0) {
		float *mzzsumgrid = NULL;
		mzzsumgrid = calloc(size[0] * size[2], sizeof(float));
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
						float junk2 = 0;
						for (int l = 0; l<size[3]; l++)
						{
							int index = closetab[index4D(size, i, j, k, l)];
							if (index == -1) {
								junk2 += -12;
							}
							else
							{
								float junk = blur[index];
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
		float junk, junk2;
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

void sumdeltas(int *size, float *deltas, float *newblur)
{
	int i, j, k;
	float temp;
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

void sum1D(int *size, float *sum, float *array, int axis)
{
	int i, j, k;
	float temp;
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



void bayes(int *size, float*denom, float *blur, float *newblur, float *dataInt, int *barr)
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

void getranges(int *size, float *blur, float *masstab, float *ccstab, float *range, int *barr)
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

void makeaxis(float *axis, int length, float min, float binsize)
{
	int i;
#pragma omp parallel for private (i), schedule(dynamic)
	for (i = 0; i<length; i++)
	{
		axis[i] = min + i*binsize;
	}
	return;
}

void write1D(char *outfile, char *suffix, float *array, int length)
{
	char outstring[500];
	int i;
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s.bin", outfile, suffix);
	out_ptr = fopen(outstring, "wb");
	fwrite(array, sizeof(float), length, out_ptr);

	//sprintf(outstring,"%s_%s.txt",outfile,suffix);
	//out_ptr=fopen(outstring,"w");
	//for(i=0;i<length;i++)
	//{
	//	  fprintf(out_ptr,"%f\n",array[i]);
	//}
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}

void write2D(char *outfile, char *suffix, float *array1, float *array2, int length)
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

void write3D(char *outfile, char *suffix, float *array1, float *array2, float *array3, int length)
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

void writemfileres(char *outfile, char *suffix, float *array1, float *array2, float *array3, float *array4, float *array5, float *array6, int length)
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

void writemzgrid(char *outfile, char *suffix, float *mzext, float *dtext, int *ztab, float *blur, int *size)
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
	fwrite(blur, sizeof(float), size[0] * size[1] * size[2], out_ptr);
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}

float errfun(int length, float *dataInt, float *fitdat)
{
	int i;
	float tot = 0;
	float mean;

	for (i = 0; i<length; i++)
	{
		tot += fitdat[i];
	}
	mean = tot / ((float)length);

	float error = 0;
	float error2 = 0;
	for (i = 0; i<length; i++)
	{
		error += pow((fitdat[i] - dataInt[i]), 2);
		error2 += pow((mean - dataInt[i]), 2);
	}
	if (error2 == 0) { printf("Error function divide by 0!"); error2 = 1; }
	return 1.0 - error / error2;
}

float getmax(int length, float *data)
{
	float max = 0;
	for (int i = 0; i<length; i++)
	{
		if (data[i]>max) { max = data[i]; }
	}
	return max;
}

void normalize(int length, float *data, float norm)
{
	if (norm == 0) { printf("Tried to Normalize by Zero! Aborting!\n"); exit(50); }
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<length; i++)
	{
		data[i] = data[i] / norm;
	}
	return;
}

void writezslice(int *size, char *outfile, char *suffix, float *massaxis, float *ccsaxis, int *ztab, float *array, int k)
{
	int i, j;
	char outstring[500];
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s_%d.bin", outfile, suffix, ztab[k]);
	out_ptr = fopen(outstring, "wb");
	float *temp = NULL;
	temp = calloc(size[0] * size[1], sizeof(float));
	for (i = 0; i<size[0]; i++)
	{
		for (j = 0; j<size[1]; j++)
		{
			//fprintf(out_ptr,"%f %f %f\n",massaxis[i],ccsaxis[j],array[index3D(size[1],size[2],i,j,k)]);
			temp[index2D(size[1], i, j)] = array[index3D(size[1], size[2], i, j, k)];
		}
	}
	fwrite(temp, sizeof(float), size[0] * size[1], out_ptr);
	fclose(out_ptr);
	free(temp);
	printf("File written to: %s\n", outstring);
}


float nativeCCS(float mass, float fudge, float gasmass)
{
	float a, b;
	if (gasmass<10)
	{
		a = 4.06739;
		b = 0.629424;
	}
	else {
		a = 5.33311;
		b = 0.613072;
	}
	//float c=-1481.7;
	return a*pow(mass, b) + fudge;
}

float bilinearinterpolation(int *size, int indm, int indc, int k, float *mzext, float *dtext, float tempmz, float tempdt, int rawflag, float *newblur, float *blur)
{
	float y11, y12, y21, y22, mu1, mu2;
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
float cubicInterpolate_IM(float p[4], float x) {
	return p[1] + 0.5 * x*(p[2] - p[0] + x*(2.0*p[0] - 5.0*p[1] + 4.0*p[2] - p[3] + x*(3.0*(p[1] - p[2]) + p[3] - p[0])));
}
float bicubicInterpolate_IM(float p[4][4], float x, float y) {
	float arr[4];
	arr[0] = cubicInterpolate_IM(p[0], y);
	arr[1] = cubicInterpolate_IM(p[1], y);
	arr[2] = cubicInterpolate_IM(p[2], y);
	arr[3] = cubicInterpolate_IM(p[3], y);
	return cubicInterpolate_IM(arr, x);
}

float bicubicinterpolation(int *size, int indm, int indc, int k, float *mzext, float *dtext, float tempmz, float tempdt, int rawflag, float *newblur, float *blur)
{
	float mu1, mu2;
	float p[4][4];
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

float WeightedAverage(int length, float *xarray, float *warray)
{
	float temp1 = 0;
	float temp2 = 0;
	int i;
	for (i = 0; i<length; i++)
	{
		temp1 += xarray[i] * warray[i];
		temp2 += warray[i];
	}
	return temp1 / temp2;
}

float WeightedAverageInt(int length, int *xarray, float *warray)
{
	float temp1 = 0;
	float temp2 = 0;
	int i;
	for (i = 0; i<length; i++)
	{
		temp1 += xarray[i] * warray[i];
		temp2 += warray[i];
	}
	return temp1 / temp2;
}

float StdDev2(int length, float *xarray, float *warray, float wmean)
{
	float temp1 = 0;
	float temp2 = 0;
	int i;
	for (i = 0; i<length; i++)
	{
		temp1 += pow(xarray[i] - wmean, 2)*warray[i];
		temp2 += warray[i];
	}
	return sqrt(temp1 / temp2);
}

float StdDev2Int(int length, int *xarray, float *warray, float wmean)
{
	float temp1 = 0;
	float temp2 = 0;
	int i;
	for (i = 0; i<length; i++)
	{
		temp1 += pow(xarray[i] - wmean, 2)*warray[i];
		temp2 += warray[i];
	}
	return sqrt(temp1 / temp2);
}


void MFileInt(int mfilelen, float *massaxis, float *massaxisval, float *testmasses, float *testmassint, float maaxle)
{
	int i, index;
	for (i = 0; i<mfilelen; i++)
	{
		index = nearfast(massaxis, testmasses[i], maaxle);
		testmassint[i] = massaxisval[index];
	}
}

//Works for both CCS and Z to get weighted average and Standard Deviation.
void MFileCCS(int mfilelen, float *massaxis, float *ccsaxis, float *massccsgrid, float *testmasses, float *testmassCCSavg, float *testmassCCSstddev, int maaxle, int ccaxle)
{
	int i, j, index;
	float *temp;
	float wmean;
	temp = calloc(ccaxle, sizeof(float));
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

void MFileZ(int mfilelen, float *massaxis, int *ztab, float *masszgrid, float *testmasses, float *testmassZavg, float *testmassZstddev, int maaxle, int numz)
{
	int i, j, index;
	float *temp;
	float wmean;
	temp = calloc(numz, sizeof(float));
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

void KillB_IM(float *I, int *B, int * size, float intthresh)
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



void ManualAssign_IM(char *manualfile, int * size, float *mzdat, float * dtdat, int *ztab, int *barr)
{

	//Manual Assignments
	int manlength;
	float *manmz = NULL,
		*manmzwin = NULL,
		*mandt = NULL,
		*mandtwin = NULL,
		*mancharge = NULL;
	manlength = getfilelength(manualfile);
	manmz = calloc(manlength, sizeof(float));
	manmzwin = calloc(manlength, sizeof(float));
	mandt = calloc(manlength, sizeof(float));
	mandtwin = calloc(manlength, sizeof(float));
	mancharge = calloc(manlength, sizeof(float));
	readfilemanual(manualfile, manlength, manmz, manmzwin, mandt, mandtwin, mancharge);
	printf("Read manual file of length: %d\n", manlength);

#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<size[0]; i++)
	{
		for (int j = 0; j<size[1]; j++)
		{
			float mz = mzdat[index2D(size[1], i, j)];
			float dt = dtdat[index2D(size[1], i, j)];
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
