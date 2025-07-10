//
// Created by mm96978 on 7/10/2025.
//

#include "UniDecIM.h"


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
	int l = size[0] * size[1];
	temp = calloc(l, sizeof(float));
	tempout = calloc(l, sizeof(float));
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

float calcCCS(const float mass, const int z, const float dt, const float ccsconst, const float hmass, const float to, const int type)
{
	float ac = 1.6605389E-27;
	float factor = 0.001;
	float rmass = ac * (mass * hmass) / (mass + hmass);
	if (type == 1) {
		ac = 1;
		factor = 1;
		rmass = mass / (mass + hmass);
	}

	float td = (dt - to)*factor;
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

float calcDt(const float mass, const int z, const float ccs, const float ccsconst, const float hmass, const float to, const int type)
{

	float ac = 1.6605389E-27;
	float factor = 1000;
	float rmass = ac * (mass * hmass) / (mass + hmass);
	if (type == 1)
	{
		ac = 1;
		factor = 1;
		rmass = mass / (mass + hmass);
	}

	float td = ccs / (ccsconst*(float)z*sqrt(1 / rmass));
	float dt = (td * factor) + to;
	return dt;
}

float calcCCSSLIMpoly3(float mass, int z, float dt, float tcal1, float tcal2, float tcal3, float tcal4, float hmass, float edc)
{
	float rdt = dt -(edc * sqrtf(mass / z) / 1000.);
	float rccs;
	if (rdt <= 0) { return 0; }
	rccs = tcal4 * powf(rdt,3) + tcal3 * powf(rdt, 2) + tcal2 * rdt+tcal1;
	float gamma = z;// (1 / (float)z)* sqrtf(mass / (mass + hmass));
	if (gamma > 0) { rccs *= gamma; }
	else { printf("CCS Gamma is negative %d %f %f %f\n", z, mass, hmass, gamma); exit(93);}
	return rccs;
}


float calcDtSLIMpoly3(float mass, int z, float ccs, float tcal1, float tcal2, float tcal3, float tcal4, float hmass, float edc)
{
	float gamma = z;// (1 / (float)z)* sqrtf(mass / (mass + hmass));
	float rccs = ccs;
	if (gamma > 0) { rccs /= gamma; }
	else { printf("Dt Gamma is negative %d %f %f %f\n", z, mass, hmass, gamma); exit(94);}
	if (rccs <= 0) { return 0; }

	float a = tcal4;
	float b = tcal3;
	float c = tcal2;
	float d = tcal1 - rccs;

	//printf("%.5e x^3 + %fx^2 + %fx + %f = 0 \n", a, b, c, d);
	float d0 = powf(b, 2) - 3 * a * c;
	float d1 = 2 * powf(b, 3) - 9 * a * b * c + 27 * powf(a, 2) * d;
	float d2 = powf(d1, 2) - 4 * powf(d0, 3);
	if (d2 < 0) { return 0; }

	float cst = cbrtf((d1+sqrtf(d2))/2.);
	if (cst == 0) { return 0; }
	float rdt = -1 / (3 * a) * (b + cst + d0 / cst);
	//printf("x = %f \n", rdt);
	/*
	//An early attempt. Doesn't seem to work very well...
	float p = -b / (3*a);
	float q = powf(p, 3) + ((b * c) - (3 * a * d))/(6 * powf(a, 2));
	float r = c / (3 * a);

	float p2 = powf(p, 2);
	float q2 = powf(q, 2);
	float rmp23 = powf(r - p2, 3);

	float _Complex c1 = csqrtf(q2 + rmp23);
	float _Complex c2 = cpowf(q + c1, 1. / 3.);
	float _Complex c3 = cpowf(q - c1, 1. / 3.);

	float _Complex crdt = c2 + c3 + p;
	float rdt = cabsf(crdt);
	//rdt = pow(q + pow(pow(q, 2) + pow(r - pow(p, 2), 3), 0.5), 1 / 3) + pow(q - pow(pow(q, 2) + pow(r - pow(p, 2), 3), 0.5), 1 / 3) + p;
	*/

	return rdt +(edc * sqrtf(mass / z) / 1000);
}

float calcCCSSLIMpoly2(float mass, int z, float dt, float tcal1, float tcal2, float tcal3,  float hmass, float edc)
{
	float rdt = dt -(edc * sqrtf(mass / z) / 1000.);
	float rccs;
	if (rdt <= 0) { return 0; }
	rccs = tcal3 * powf(rdt,2) + tcal2 *rdt+tcal1;
	float gamma = z;// (1 / (float)z)* sqrtf(mass / (mass + hmass));
	if (gamma > 0) { rccs *= gamma; }
	else { printf("CCS Gamma is negative %d %f %f %f\n", z, mass, hmass, gamma); exit(91); }
	return rccs;
}

float calcDtSLIMpoly2(float mass, int z, float ccs, float tcal1, float tcal2, float tcal3,  float hmass, float edc)
{
	float gamma = z;// (1 / (float)z)* sqrtf(mass / (mass + hmass));
	float rccs = ccs;
	if (gamma > 0) { rccs /= gamma; }
	else { printf("DT Gamma is negative %d %f %f %f\n", z, mass, hmass, gamma); exit(92); }
	if (rccs <= 0) { return 0; }

	float a = tcal3;
	float b = tcal2;
	float c = tcal1-rccs;

	float d0 = powf(b, 2) - 4 * a * c;
	if (d0 < 0) { return 0; }
	float d1 = sqrtf(d0);

	//printf(" %fx^2 + %fx + %f = 0 \n", a, b, c);

	float rdt = (-b + d1)/(2*a);
	//printf("x = %f \n", rdt);
	if (rdt > 0) { return rdt; }// + (edc * sqrt(mass / z) / 1000);}

	rdt = (-b - d1) / (2 * a);
	//printf("x = %f \n", rdt);
	if (rdt > 0) { return rdt; }// +(edc * sqrt(mass / z) / 1000); }

	return 0;
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

void blur_it_IM(int *size, float *blur, float *newblur, int *closetab, int *barr, float csig)
{

	if (csig < 0) {
		float *mzzsumgrid = NULL;
		int l = size[0] * size[2];
		mzzsumgrid = calloc(l, sizeof(float));
		Sum2D(size, mzzsumgrid, blur, 1);
		int i, j, k, h;
#pragma omp parallel for private (i,j,k)
		for (i = 0; i < size[0]; i++)
		{
			for (j = 0; j < size[1]; j++)
			{
				for (k = 0; k < size[2]; k++)
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

#pragma omp parallel for private (i,j,k,h)
		for (i = 0; i<size[0]; i++)
		{
			for (j = 0; j<size[1]; j++)
			{
				for (k = 0; k<size[2]; k++)
				{
					if (barr[index3D(size[1], size[2], i, j, k)] == 1)
					{
						float junk2 = 0;
						for (h = 0; h<size[3]; h++)
						{
							int index = closetab[index4D(size, i, j, k, h)];
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
#pragma omp parallel for private (i,j,k,l,index,junk,junk2)
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
		exit(1);
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


void writemfileres(char *outfile, char *suffix, float *array1, float *array2, float *array3, float *array4, float *array5, float *array6, int length)
{
	char outstring[500];
	int i;
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s.txt", outfile, suffix);
	out_ptr = fopen(outstring, "w");
	if (out_ptr == 0) { printf("Error Opening %s\n", outstring); exit(1); }
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
	if (out_ptr == 0) { printf("Error Opening %s\n", outstring); exit(1); }
	int l = size[0] * size[1] * size[2];
	fwrite(blur, sizeof(float), l, out_ptr);
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

void writezslice(int *size, char *outfile, char *suffix, float *massaxis, float *ccsaxis, int *ztab, float *array, int k)
{
	int i, j;
	char outstring[500];
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s_%d.bin", outfile, suffix, ztab[k]);
	out_ptr = fopen(outstring, "wb");
	if (out_ptr == 0) { printf("Error Opening %s\n", outstring); exit(1); }
	float *temp = NULL;
	int newlen = size[0] * size[1];
	temp = calloc(newlen, sizeof(float));
	for (i = 0; i<size[0]; i++)
	{
		for (j = 0; j<size[1]; j++)
		{
			//fprintf(out_ptr,"%f %f %f\n",massaxis[i],ccsaxis[j],array[index3D(size[1],size[2],i,j,k)]);
			temp[index2D(size[1], i, j)] = array[index3D(size[1], size[2], i, j, k)];
		}
	}
	int l = size[0] * size[1];
	fwrite(temp, sizeof(float), l, out_ptr);
	fclose(out_ptr);
	free(temp);
	printf("File written to: %s\n", outstring);
}


float nativeCCS(const float mass, const float fudge, const float gasmass)
{
	float a, b;
	if (gasmass<10)
	{
		a = 4.06739f;
		b = 0.629424f;
	}
	else {
		a = 5.33311f;
		b = 0.613072f;
	}
	//float c=-1481.7;
	return a*powf(mass, b) + fudge;
}

float bilinearinterpolation(const int *size, const int indm, const int indc, const int k, const float *mzext,
	const float *dtext, const float tempmz, const float tempdt, const int rawflag, const float *newblur, const float *blur)
{
	float y11, y12, y21, y22;
	int im1, im2, ic1, ic2;
	if (mzext[indm]>tempmz) { im2 = indm; im1 = indm - 1; }
	else { im1 = indm; im2 = indm + 1; }
	if (dtext[indc]>tempdt) { ic2 = indc; ic1 = indc - 1; }
	else { ic1 = indc; ic2 = indc + 1; }
	const float mu1 = (tempmz - mzext[im1]) / (mzext[im2] - mzext[im1]);
	const float mu2 = (tempdt - dtext[ic1]) / (dtext[ic2] - dtext[ic1]);
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
	return clip(y11*(1 - mu1)*(1 - mu2) + y12*mu1*(1 - mu2) + y21*mu2*(1 - mu1) + y22*mu1*mu2, 0);
}

//cubicInterpolate_IM and bicubicInterpolate_IM were taken from http://www.paulinternet.nl/?page=bicubic
float cubicInterpolate_IM(float p[4], const float x) {
	return p[1] + 0.5f * x*(p[2] - p[0] + x*(2.0f*p[0] - 5.0f*p[1] + 4.0f*p[2] - p[3] + x*(3.0f*(p[1] - p[2]) + p[3] - p[0])));
}

float bicubicInterpolate_IM(float p[4][4], const float x, const float y) {
	float arr[4] = { 0, 0 ,0, 0 };
	arr[0] = cubicInterpolate_IM(p[0], y);
	arr[1] = cubicInterpolate_IM(p[1], y);
	arr[2] = cubicInterpolate_IM(p[2], y);
	arr[3] = cubicInterpolate_IM(p[3], y);
	return cubicInterpolate_IM(arr, x);
}

float bicubicinterpolation(const int *size, const int indm, const int indc, const int k, const float *mzext,
	const float *dtext, const float tempmz, const float tempdt, const int rawflag, const float *newblur, const float *blur)
{
	float p[4][4] = { { 0,0,0,0 } , { 0,0,0,0 } , { 0,0,0,0 } , { 0,0,0,0 } };
	int im[4] = { 0,0,0,0 }, ic[4] = { 0,0,0,0 };
	if (mzext[indm]>tempmz) { im[2] = indm; im[1] = indm - 1; im[0] = indm - 2; im[3] = indm + 1; }
	else { im[1] = indm; im[2] = indm + 1; im[0] = indm - 1; im[3] = indm + 2; }
	if (dtext[indc]>tempdt) { ic[2] = indc; ic[1] = indc - 1; ic[0] = indc - 2; ic[3] = indc + 1; }
	else { ic[1] = indc; ic[2] = indc + 1; ic[0] = indc - 1; ic[3] = indc + 2; }
	const float mu1 = (tempmz - mzext[im[1]]) / (mzext[im[2]] - mzext[im[1]]);
	const float mu2 = (tempdt - dtext[ic[1]]) / (dtext[ic[2]] - dtext[ic[1]]);
	for (int i = 0; i<4; i++)
	{
		for (int j = 0; j<4; j++)
		{
			if (rawflag == 0) { p[i][j] = newblur[index3D(size[1], size[2], im[i], ic[j], k)]; }
			else { p[i][j] = blur[index3D(size[1], size[2], im[i], ic[j], k)]; }
		}
	}
	return clip(bicubicInterpolate_IM(p, mu1, mu2), 0);
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
		testmassCCSstddev[i] = WeightedStdDev(ccaxle, ccsaxis, temp, wmean);
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
		testmassZstddev[i] = WeightedStdDevInt(numz, ztab, temp, wmean);
	}
	free(temp);
}

void KillB_IM(float *IntArray, int *B, int * size, float intthresh)
{
	for (int i = 0; i<size[0]; i++)
	{
		for (int j = 0; j<size[1]; j++)
		{
			for (int k = 0; k<size[2]; k++)
			{
				if (IntArray[index2D(size[1], i, j)] < intthresh) {
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
	int i,j;
#pragma omp parallel for private(i,j)
	for (i = 0; i<size[0]; i++)
	{
		for (j = 0; j<size[1]; j++)
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

void TwaveError(const int twaveflag)
{
	printf("Error: Undefined twaveflag. Value was %d. Aborting.", twaveflag);
	exit(20);
}