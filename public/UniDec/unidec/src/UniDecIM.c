//
// Created by mm96978 on 7/10/2025.
//

#include "UniDecIM.h"


void fftconvolve2D(float *corr, const float *list1, const float *list2, const int *size)
{
	const int length = size[0] * size[1];
	fftwf_complex *in1 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * length);
	fftwf_complex *in2 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * length);
	fftwf_complex *out1 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * length);
	fftwf_complex *out2 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * length);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<length; i++)
	{
		in1[i][0] = list1[i];
		in1[i][1] = 0;
		in2[i][0] = list2[i];
		in2[i][1] = 0;
	}
	fftwf_plan p1 = fftwf_plan_dft_2d(size[0], size[1], in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_plan p2 = fftwf_plan_dft_2d(size[0], size[1], in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(p1);
	fftwf_execute(p2);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<length; i++)
	{
		const float a = out1[i][0];
		const float b = out1[i][1];
		const float c = out2[i][0];
		const float d = out2[i][1];
		in1[i][0] = a*c - b*d;
		in1[i][1] = b*c + a*d;
	}

	fftwf_plan p3 = fftwf_plan_dft_2d(size[0], size[1], in1, out1, FFTW_BACKWARD, FFTW_ESTIMATE);
	fftwf_execute(p3);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<length; i++)
	{
		corr[i] = out1[i][0] / (float) length;
	}
	//
	fftwf_destroy_plan(p1);
	fftwf_destroy_plan(p2);
	fftwf_destroy_plan(p3);
	fftwf_free(in1); fftwf_free(in2); fftwf_free(out1); fftwf_free(out2);
}

void convolve3D(float *out, const float *in, const float *peakshape, const int *size)
{
	//Initialize memory
	const int l = size[0] * size[1];
	float *temp = calloc(l, sizeof(float));
	float *tempout = calloc(l, sizeof(float));
	if (temp == NULL || tempout == NULL) {
		fprintf(stderr, "Memory allocation failed in convolve3D\n");
		exit(EXIT_FAILURE);
	}
	//Convolve for each charge slice
	for (int k = 0; k<size[2]; k++)
	{
		//Cut out a slice
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i<size[0]; i++)
		{
			for (int j = 0; j<size[1]; j++)
			{
				temp[index2D(size[1], i, j)] = in[index3D(size[1], size[2], i, j, k)];
			}
		}
		//Convolve it
		fftconvolve2D(tempout, temp, peakshape, size);
		//Repack slice into output array
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i<size[0]; i++)
		{
			for (int j = 0; j<size[1]; j++)
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
	const float ac = 1.6605389E-27f;
	float factor = 0.001f;
	float rmass = ac * (mass * hmass) / (mass + hmass);
	if (type == 1) {
		factor = 1;
		rmass = mass / (mass + hmass);
	}

	const float td = (dt - to)*factor;
	const float ccs = (float)z*(td)*sqrtf(1 / rmass)*ccsconst;
	return ccs;
}

float calcCCSTwaveLog(const float mass, const int z, const float dt, const float tcal1, const float tcal2, const float hmass, const float edc)
{
	const float rdt = dt - (edc*sqrtf(mass / (float) z) / 1.000E3f);
	if (rdt <= 0) { return 0; }
	const float rccs = expf(tcal1*logf(rdt) + tcal2);
	const float rmass = (mass*hmass) / (mass + hmass);
	return rccs* (float) z *sqrtf(1 / rmass);
}

float calcCCSTwaveLinear(const float mass, const int z, const float dt, const float tcal1, const float tcal2, const float hmass, const float edc)
{
	const float rdt = dt - (edc*sqrtf(mass / (float) z) / 1.000E3f);
	if (rdt <= 0) { return 0; }
	const float rccs = tcal1*rdt + tcal2;
	const float rmass = (mass*hmass) / (mass + hmass);
	return rccs* (float) z*sqrtf(1 / rmass);
}

float calcCCSTwavePower(const float mass, const int z, const float dt, const float tcal1, const float tcal2, const float hmass, const float edc)
{
	const float rdt = dt - (edc*sqrtf(mass / (float) z) / 1.000E3f);
	if (rdt <= 0) { return 0; }
	const float rccs = tcal1*powf(rdt, tcal2);
	const float rmass = (mass*hmass) / (mass + hmass);
	return rccs*(float) z*sqrtf(1 / rmass);
}

float calcDtTwaveLog(const float mass, const int z, const float ccs, const float tcal1, const float tcal2, const float hmass, const float edc)
{
	//float ac=1.6605389E-27;
	const float rmass = (mass*hmass) / (mass + hmass);
	const float rccs = ccs / ((float) z*sqrtf(1 / rmass));
	if (rccs <= 0) { return 0; }
	const float rdt = expf((logf(rccs) - tcal2) / tcal1);
	return rdt + (edc*sqrtf(mass / (float) z) / 1000);
}

float calcDtTwaveLinear(const float mass, const int z, const float ccs, const float tcal1, const float tcal2, const float hmass, const float edc)
{
	const float rmass = (mass*hmass) / (mass + hmass);
	const float rccs = ccs / ((float) z*sqrtf(1 / rmass));

	if (rccs <= 0) { return 0; }
	const float rdt = (rccs - tcal2) / tcal1;
	return rdt + (edc*sqrtf(mass / (float) z) / 1000);
}

float calcDtTwavePower(const float mass, const int z, const float ccs, const float tcal1, const float tcal2, const float hmass, const float edc)
{
	const float rmass = (mass*hmass) / (mass + hmass);
	const float rccs = ccs / ((float) z*sqrtf(1 / rmass));

	if (rccs <= 0) { return 0; }
	const float rdt = powf(rccs / tcal1, 1.0f / tcal2);
	return rdt + (edc*sqrtf(mass / (float) z) / 1000);
}

float calcDt(const float mass, const int z, const float ccs, const float ccsconst, const float hmass, const float to, const int type)
{
	//printf("Test: %f %d %f %f %f %f %d\n", mass, z, ccs, ccsconst, hmass, to, type);
	const float ac = 1.6605389E-27f;
	float factor = 1000;
	const float denom1 = mass + hmass;
	if (denom1 ==0){return -1;}
	float rmass = ac * (mass * hmass) / denom1;
	if (type == 1)
	{
		factor = 1;
		rmass = mass / denom1;
	}
	if (rmass==0){return -1;}
	const float denom2 = ccsconst*(float)z*sqrtf(1 / rmass);
	if (denom2 == 0) { return -1; }
	const float td = ccs / denom2;
	const float dt = (td * factor) + to;
	return dt;
}

float calcCCSSLIMpoly3(const float mass, const int z, const float dt, const float tcal1, const float tcal2, const float tcal3, const float tcal4, const float hmass, const float edc)
{
	const float rdt = dt -(edc * sqrtf(mass / (float) z) / 1000.f);
	if (rdt <= 0) { return 0; }
	float rccs = tcal4 * powf(rdt,3) + tcal3 * powf(rdt, 2) + tcal2 * rdt+tcal1;
	const float gamma = (float) z;// (1 / (float)z)* sqrtf(mass / (mass + hmass));
	if (gamma > 0) { rccs *= gamma; }
	else { printf("CCS Gamma is negative %d %f %f %f\n", z, mass, hmass, gamma); exit(93);}
	return rccs;
}


float calcDtSLIMpoly3(const float mass, const int z, const float ccs, const float tcal1, const float tcal2, const float tcal3, const float tcal4, const float hmass, const float edc)
{
	const float gamma = (float) z;// (1 / (float)z)* sqrtf(mass / (mass + hmass));
	float rccs = ccs;
	if (gamma > 0) { rccs /= gamma; }
	else { printf("Dt Gamma is negative %d %f %f %f\n", z, mass, hmass, gamma); exit(94);}
	if (rccs <= 0) { return 0; }

	const float a = tcal4;
	const float b = tcal3;
	const float c = tcal2;
	const float d = tcal1 - rccs;

	//printf("%.5e x^3 + %fx^2 + %fx + %f = 0 \n", a, b, c, d);
	const float d0 = powf(b, 2) - 3 * a * c;
	const float d1 = 2 * powf(b, 3) - 9 * a * b * c + 27 * powf(a, 2) * d;
	const float d2 = powf(d1, 2) - 4 * powf(d0, 3);
	if (d2 < 0) { return 0; }

	const float cst = cbrtf((d1+sqrtf(d2))/2.f);
	if (cst == 0) { return 0; }
	const float rdt = -1 / (3 * a) * (b + cst + d0 / cst);
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

	return rdt +(edc * sqrtf(mass / (float) z) / 1000);
}

float calcCCSSLIMpoly2(const float mass, const int z, const float dt, const float tcal1, const float tcal2, const float tcal3, const float hmass, const float edc)
{
	const float rdt = dt -(edc * sqrtf(mass / (float) z) / 1000.f);
	if (rdt <= 0) { return 0; }
	float rccs = tcal3 * powf(rdt,2) + tcal2 *rdt+tcal1;
	float gamma = (float) z;// (1 / (float)z)* sqrtf(mass / (mass + hmass));
	if (gamma > 0) { rccs *= gamma; }
	else { printf("CCS Gamma is negative %d %f %f %f\n", z, mass, hmass, gamma); exit(91); }
	return rccs;
}

float calcDtSLIMpoly2(const float mass, const int z, const float ccs, const float tcal1, const float tcal2, const float tcal3, const float hmass, float edc)
{
	const float gamma = (float) z;// (1 / (float)z)* sqrtf(mass / (mass + hmass));
	float rccs = ccs;
	if (gamma > 0) { rccs /= gamma; }
	else { printf("DT Gamma is negative %d %f %f %f\n", z, mass, hmass, gamma); exit(92); }
	if (rccs <= 0) { return 0; }

	const float a = tcal3;
	const float b = tcal2;
	const float c = tcal1-rccs;

	const float d0 = powf(b, 2) - 4 * a * c;
	if (d0 < 0) { return 0; }
	const float d1 = sqrtf(d0);

	//printf(" %fx^2 + %fx + %f = 0 \n", a, b, c);

	float rdt = (-b + d1)/(2*a);
	//printf("x = %f \n", rdt);
	if (rdt > 0) { return rdt; }// + (edc * sqrt(mass / z) / 1000);}

	rdt = (-b - d1) / (2 * a);
	//printf("x = %f \n", rdt);
	if (rdt > 0) { return rdt; }// +(edc * sqrt(mass / z) / 1000); }

	return 0;
}


void blur_it_IM(const int *size, float *blur, float *newblur, const int *closetab, const int *barr, const float csig)
{
	if (csig < 0) {
		const int l = size[0] * size[2];
		float *mzzsumgrid = calloc(l, sizeof(float));
		Sum2D(size, mzzsumgrid, blur, 1);
		#pragma omp parallel for
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

		#pragma omp parallel for
		for (int i = 0; i<size[0]; i++)
		{
			for (int j = 0; j<size[1]; j++)
			{
				for (int k = 0; k<size[2]; k++)
				{
					if (barr[index3D(size[1], size[2], i, j, k)] == 1)
					{
						float junk2 = 0;
						for (int h = 0; h<size[3]; h++)
						{
							const int index = closetab[index4D(size, i, j, k, h)];
							if (index == -1) {
								junk2 += -12;
							}
							else
							{
								const float junk = blur[index];
								if (junk>0) { junk2 += logf(junk); }
								else { junk2 += -12; }
							}
						}
						newblur[index3D(size[1], size[2], i, j, k)] = expf(junk2 / (float) size[3]);
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
		#pragma omp parallel for
		for (int i = 0; i < size[0]; i++)
		{
			for (int j = 0; j < size[1]; j++)
			{
				for (int k = 0; k < size[2]; k++)
				{
					if (barr[index3D(size[1], size[2], i, j, k)] == 1)
					{
						float junk2 = 0;
						for (int l = 0; l < size[3]; l++)
						{
							const int index = closetab[index4D(size, i, j, k, l)];
							if (index == -1) {
								junk2 += -12;
								// else { junk2 += 0; }
							}
							else
							{
								const float junk = blur[index];
								// if (logswitch == 0) { junk2 += junk; }
								// else {
								if (junk > 0) { junk2 += logf(junk); }
								else { junk2 += -12; }

							}
						}
						// if (logswitch == 0) { newblur[index3D(size[1], size[2], i, j, k)] = junk2 / size[3]; }
						newblur[index3D(size[1], size[2], i, j, k)] = expf(junk2 / (float) size[3]);
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

void sumdeltas(const int *size, float *deltas, const float *newblur)
{
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<size[0]; i++)
	{
		for (int j = 0; j<size[1]; j++)
		{
			float temp = 0;
			for (int k = 0; k<size[2]; k++)
			{
				temp += newblur[index3D(size[1], size[2], i, j, k)];
			}
			deltas[index2D(size[1], i, j)] = temp;
		}
	}
}


void bayes(const int *size, const float*denom, float *blur, const float *newblur, const float *dataInt, const int *barr)
{
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<size[0]; i++)
	{
		for (int j = 0; j<size[1]; j++)
		{
			for (int k = 0; k<size[2]; k++)
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
}

void getranges(const int *size, const float *blur, const float *masstab, const float *ccstab, float *range, const int *barr)
{
	range[0] = 100000000.f;
	range[1] = 0;
	range[2] = 100000000.f;
	range[3] = 0;
	for (int i = 0; i<size[0]; i++)
	{
		for (int j = 0; j<size[1]; j++)
		{
			for (int k = 0; k<size[2]; k++)
			{
				if (blur[index3D(size[1], size[2], i, j, k)]>1E-12&&barr[index3D(size[1], size[2], i, j, k)] == 1)
				{
					const float mass = masstab[index2D(size[2], i, k)];
					const float ccs = ccstab[index3D(size[1], size[2], i, j, k)];
					if (mass<range[0]) { range[0] = mass; }
					if (mass>range[1]) { range[1] = mass; }
					if (ccs<range[2]) { range[2] = ccs; }
					if (ccs>range[3]) { range[3] = ccs; }
				}
			}
		}
	}
}

void makeaxis(float *axis, const int length, const float min, const float binsize)
{
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i<length; i++)
	{
		axis[i] = min + (float) i*binsize;
	}
}


void readfilemanual(char *infile, const int lengthmz, float *array1, float *array2, float *array3, float *array4, float *array5)
{
	char a[500];
	char b[500];
	char c[500];
	char d[500];
	char e[500];
	FILE *file_ptr = fopen(infile, "r");
	if (file_ptr == 0)
	{
		printf("Could not open %s file\n", infile);
		exit(1);
	}

	for (int i = 0; i<lengthmz; i++)
	{
		// replacement for the fscanf line in `UniDecIM.c`
		const int ret = fscanf(file_ptr, "%499s %499s %499s %499s %499s", a, b, c, d, e);
		if (ret != 5) {
			if (ret == EOF) {
				fprintf(stderr, "Unexpected end of file while reading `infile` at entry %d\n", i);
				break; // or handle EOF appropriately
			}
			else {
				fprintf(stderr, "Parse error in `infile` at entry %d (fields read: %d)\n", i, ret);
				a[0]=b[0]=c[0]=d[0]=e[0]='\0'; // set defaults so strtof yields 0.0f
			}
		}
		char * endptr;
		array1[i] = strtof(a, &endptr);
		array2[i] = strtof(b, &endptr);
		array3[i] = strtof(c, &endptr);
		array4[i] = strtof(d, &endptr);
		array5[i] = strtof(e, &endptr);
	}

	fclose(file_ptr);
}


void writemfileres(char *outfile, char *suffix, const float *array1, const float *array2, const float *array3, const float *array4, const float *array5, const float *array6, const int length)
{
	char outstring[500];
	sprintf(outstring, "%s_%s.txt", outfile, suffix);
	FILE *out_ptr = fopen(outstring, "w");
	if (out_ptr == 0) { printf("Error Opening %s\n", outstring); exit(1); }
	//fprintf(out_ptr,"Mass\tIntensity\tCCS Avg.\tCCS Std. Dev.\tZ avg.\tZ Std. Dev.\n");
	for (int i = 0; i<length; i++)
	{
		fprintf(out_ptr, "%f\t%f\t%f\t%f\t%f\t%f\n", array1[i], array2[i], array3[i], array4[i], array5[i], array6[i]);
	}
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}

void writemzgrid(char *outfile, char *suffix, const float *blur, const int *size)
{
	char outstring[500];
	FILE *out_ptr = NULL;
	sprintf(outstring, "%s_%s.bin", outfile, suffix);
	//printf("%s\n",outstring);
	out_ptr = fopen(outstring, "wb");
	if (out_ptr == 0) { printf("Error Opening %s\n", outstring); exit(1); }
	const int l = size[0] * size[1] * size[2];
	fwrite(blur, sizeof(float), l, out_ptr);
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}

float errfun(const int length, const float *dataInt, const float *fitdat)
{
	float tot = 0;
	for (int i = 0; i<length; i++)
	{
		tot += fitdat[i];
	}
	const float mean = tot / ((float)length);

	float error = 0;
	float error2 = 0;
	for (int i = 0; i<length; i++)
	{
		error += powf((fitdat[i] - dataInt[i]), 2);
		error2 += powf((mean - dataInt[i]), 2);
	}
	if (error2 == 0) { printf("Error function divide by 0!"); error2 = 1; }
	return 1.0f - error / error2;
}

void writezslice(const int *size, char *outfile, char *suffix, const int *ztab, const float *array, const int k)
{
	char outstring[500];
	sprintf(outstring, "%s_%s_%d.bin", outfile, suffix, ztab[k]);
	FILE *out_ptr = fopen(outstring, "wb");
	if (out_ptr == 0) { printf("Error Opening %s\n", outstring); exit(1); }
	const int newlen = size[0] * size[1];
	float *temp = calloc(newlen, sizeof(float));
	for (int i = 0; i<size[0]; i++)
	{
		for (int j = 0; j<size[1]; j++)
		{
			//fprintf(out_ptr,"%f %f %f\n",massaxis[i],ccsaxis[j],array[index3D(size[1],size[2],i,j,k)]);
			temp[index2D(size[1], i, j)] = array[index3D(size[1], size[2], i, j, k)];
		}
	}
	const int l = size[0] * size[1];
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


void MFileInt(const int mfilelen, const float *massaxis, const float *massaxisval, const float *testmasses, float *testmassint, const int maaxle)
{
	for (int i = 0; i<mfilelen; i++)
	{
		const int index = nearfast(massaxis, testmasses[i], maaxle);
		testmassint[i] = massaxisval[index];
	}
}

//Works for both CCS and Z to get weighted average and Standard Deviation.
void MFileCCS(const int mfilelen, const float *massaxis, const float *ccsaxis, const float *massccsgrid,
	const float *testmasses, float *testmassCCSavg, float *testmassCCSstddev, const int maaxle, const int ccaxle)
{
	float *temp = calloc(ccaxle, sizeof(float));
	for (int i = 0; i<mfilelen; i++)
	{
		const int index = nearfast(massaxis, testmasses[i], maaxle);
		for (int j = 0; j<ccaxle; j++)
		{
			temp[j] = massccsgrid[index2D(ccaxle, index, j)];
		}
		const float wmean = WeightedAverage(ccaxle, ccsaxis, temp);
		testmassCCSavg[i] = wmean;
		testmassCCSstddev[i] = WeightedStdDev(ccaxle, ccsaxis, temp, wmean);
	}
	free(temp);
}

void MFileZ(const int mfilelen, const float *massaxis, const int *ztab, const float *masszgrid, const float *testmasses,
	float *testmassZavg, float *testmassZstddev, const int maaxle, const int numz)
{
	float *temp = calloc(numz, sizeof(float));
	for (int i = 0; i<mfilelen; i++)
	{
		const int index = nearfast(massaxis, testmasses[i], maaxle);
		for (int j = 0; j<numz; j++)
		{
			temp[j] = masszgrid[index2D(numz, index, j)];
		}
		const float wmean = WeightedAverageInt(numz, ztab, temp);
		testmassZavg[i] = wmean;
		testmassZstddev[i] = WeightedStdDevInt(numz, ztab, temp, wmean);
	}
	free(temp);
}

void KillB_IM(const float *IntArray, int *B, const int * size, const float intthresh)
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



void ManualAssign_IM(char *manualfile, const int * size, const float *mzdat, const float * dtdat, const int *ztab, int *barr)
{
	//Manual Assignments
	const int manlength = getfilelength(manualfile);
	float *manmz = calloc(manlength, sizeof(float));
	float *manmzwin = calloc(manlength, sizeof(float));
	float *mandt = calloc(manlength, sizeof(float));
	float *mandtwin = calloc(manlength, sizeof(float));
	float *mancharge = calloc(manlength, sizeof(float));
	if (manmz == NULL || manmzwin == NULL || mandt == NULL || mandtwin == NULL || mancharge == NULL)
	{
		printf("Error allocating memory for manual assignments.\n");
		exit(1);
	}
	readfilemanual(manualfile, manlength, manmz, manmzwin, mandt, mandtwin, mancharge);
	printf("Read manual file of length: %d\n", manlength);
	#pragma omp parallel for schedule(auto)
	for (int i = 0; i<size[0]; i++)
	{
		for (int j = 0; j<size[1]; j++)
		{
			const float mz = mzdat[index2D(size[1], i, j)];
			const float dt = dtdat[index2D(size[1], i, j)];
			const int close = nearunsorted2D(manmz, mandt, mz, dt, manlength);
			const int charge = (int) mancharge[close];
			if (fabsf(manmz[close] - mz) < manmzwin[close] && fabsf(mandt[close] - dt) < mandtwin[close])
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