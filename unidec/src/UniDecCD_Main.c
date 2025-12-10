//
// Created by mm96978 on 7/9/2025.
//

#include "UniDecCD_Main.h"


void blur_it_CD(float * output, const float * input, const int* upinds, const int *loinds, const int length, const float floor)
{
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < length; i++) {
		const int lowindex = loinds[i];
		const int uppindex = upinds[i];
		float i1 = input[i];
		float i2 = input[lowindex];
		float i3 = input[uppindex];
		if (floor > 0) {
			i1 = logf(i1 + floor);
			i2 = logf(i2 + floor);
			i3 = logf(i3 + floor);
			float newval = (i1 + i2 + i3) / 3;
			newval = expf(newval) - floor;
			if (newval < 0) {
				newval = 0;
			}
			output[i] = newval;
		}
		else {
			const float ratio = fabsf(floor);
			output[i] = (i1 + i2 * ratio + i3 * ratio) / 3;
		}
	}
}


void setup_blur_z(int * zupind, int * zloind, const float *mzdat, const float *zdat, const int lines, const float adductmass, const float mzranges[4], const int size[3])
{
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < lines; i++) {
		const float mz = mzdat[i];
		const float z = zdat[i];
		const int zind = (int) (z - mzranges[2]);
		const int zint = (int) z;
		const float mass = calcmass(mz, zint, adductmass);

		float uppermz = calcmz(mass, zint + 1, adductmass);

		float lowermz = 0;
		if (z - 1 != 0) {
			lowermz = calcmz(mass, zint - 1, adductmass);
		}

		const int mzupind = (int) roundf((uppermz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (float) (size[0] - 1));
		if (z == mzranges[3] || mzupind < 0 || mzupind > size[0]) {
			zupind[i] = i;
		}
		else
		{
			int newind = index2D(size[1], mzupind, zind + 1);
			if (newind < 0) { newind = 0; }
			if (newind >= lines) { newind = lines - 1; }
			zupind[i] = newind;
		}

		const int mzloind = (int) roundf((lowermz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (float) (size[0] - 1));
		if (zind == 0 || mzloind < 0 || mzloind > size[0]) {
			zloind[i] = i;
		}
		else
		{
			int newind = index2D(size[1], mzloind, zind - 1);
			if (newind < 0) { newind = 0; }
			if (newind >= lines) { newind = lines - 1; }
			zloind[i] = newind;
		}
	}
}


void setup_blur_m(int* mupind, int* mloind, const float* mzdat, const float* zdat, const int lines, const float adductmass, const float mzranges[4], const int size[3], const float molig)
{
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < lines; i++) {
		const float mz = mzdat[i];
		const float z = zdat[i];
		const int zint = (int) z;
		const int zind = (int) (z - mzranges[2]);
		const float mass = calcmass(mz, zint, adductmass);

		float uppermz = calcmz(mass+molig, zint, adductmass);
		float lowermz = calcmz(mass-molig, zint, adductmass);

		const int mzupind = (int)roundf((uppermz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (float) (size[0] - 1));
		if (mzupind < 0 || mzupind > size[0]) {
			mupind[i] = i;
		}
		else
		{
			int newind = index2D(size[1], mzupind, zind);
			if (newind < 0) { newind = 0; }
			if (newind >= lines) { newind = lines - 1; }
			mupind[i] = newind;
		}

		const int mzloind = (int)roundf((lowermz - mzranges[0]) / (mzranges[1] - mzranges[0]) * (float) (size[0] - 1));
		if (mzloind < 0 || mzloind > size[0]) {
			mloind[i] = i;
		}
		else
		{
			int newind = index2D(size[1], mzloind, zind);
			if (newind < 0) { newind = 0; }
			if (newind >= lines) { newind = lines - 1; }
			mloind[i] = newind;
		}
	}
}

void precompute_fft2D(const float* list2, const int* size, fftwf_complex * out2) {
	const int length = size[0] * size[1];
	fftwf_complex *in2 = (fftwf_complex *) fftwf_malloc(sizeof(fftwf_complex) * length);
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < length; i++)
	{
		in2[i][0] = list2[i];
		in2[i][1] = 0;
	}
	fftwf_plan p2 = fftwf_plan_dft_2d(size[0], size[1], in2, out2, FFTW_FORWARD, FFTW_ESTIMATE);
	fftwf_execute(p2);
	fftwf_destroy_plan(p2);
	fftwf_free(in2);
}

void fftconvolve2D_precomputed(float* corr, const float* list1, const fftwf_complex * out2, const int* size, fftwf_plan p1, fftwf_plan p3, fftwf_complex * in1, const fftwf_complex * out1)
{
	const int length = size[0] * size[1];

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < length; i++)
	{
		in1[i][0] = list1[i];
		in1[i][1] = 0;
	}

	fftwf_execute(p1);

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < length; i++)
	{
		const float a = out1[i][0];
		const float b = out1[i][1];
		const float c = out2[i][0];
		const float d = out2[i][1];
		in1[i][0] = a * c - b * d;
		in1[i][1] = b * c + a * d;
	}

	fftwf_execute(p3);

	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < length; i++)
	{
		corr[i] = out1[i][0] / (float) length;
	}
}

void complex_conjugate(const fftwf_complex* in, fftwf_complex* out, const int length)
{
	#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < length; i++)
	{
		out[i][0] = in[i][0];
		out[i][1] = -in[i][1];
	}
}


int run_unidec_CD(int argc, char* argv[], Config config) {

	time_t starttime, endtime;
	starttime = time(NULL);

	printf("Opening File: %s\n", config.infile);
	int lines=0;
	lines = getfilelength(config.infile);

	printf("Length of data: %d\n", lines);

	int size[3] = { 0,0,0 };

	int* zupind = NULL, * zloind = NULL;
	int* mupind = NULL, * mloind = NULL;
	char* barr = NULL;

	float* mzdat = NULL, * mzext = NULL,
		* zdat = NULL, * zext = NULL,
		* dataInt = NULL,
		* peakshape = NULL,
		*mkernel = NULL,
		*blur = NULL, *newblur=NULL, *newblur2=NULL, *oldblur=NULL
		;

	//Reading In Data
	mzdat = calloc(lines, sizeof(float));
	zdat = calloc(lines, sizeof(float));
	dataInt = calloc(lines, sizeof(float));
	if (mzdat == NULL || zdat == NULL || dataInt == NULL) {
		printf("Error allocating memory for data arrays\n");
		exit(1);
	}
	readfile3(config.infile, lines, mzdat, zdat, dataInt);

	//Determine the maximum intensity in the data
	float dmax = Max(dataInt, lines);
	//float betafactor = 1;
	//if (dmax > 1) { betafactor = dmax; }


	//Getting Dimension Sizes of Data
	size[0] = GetSize0(mzdat, lines);
	size[1] = GetSize1(zdat, lines);
	int totlen = size[0] * size[1] ;
	size[2] = totlen;
	printf("Dimensions of data: %d mz by %d z: %d total\n", size[0], size[1], totlen);
	if (totlen > 155E6) { printf("Warning: May exceed system memory capacity\n"); }

	// Set Up FFT stuff
	fftwf_complex* in1, * out1;
	fftwf_plan p1, p3;
	in1 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * lines);
	out1 = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * lines);
	p1 = fftwf_plan_dft_2d(size[0], size[1], in1, out1, FFTW_FORWARD, FFTW_ESTIMATE);
	p3 = fftwf_plan_dft_2d(size[0], size[1], in1, out1, FFTW_BACKWARD, FFTW_ESTIMATE);

	//Extracting mz and z ranges
	mzext = calloc(size[0], sizeof(float));
	zext = calloc(size[1], sizeof(float));
	if (mzext == NULL || zext == NULL) {
		printf("Error allocating memory for mz and z extents\n");
		exit(1);
	}
	PullXY(mzext, zext, mzdat, zdat, size);
	float mzranges[4] = { 0,0,0,0 };
	mzranges[0] = mzext[0];
	mzranges[1] = mzext[size[0] - 1];
	mzranges[2] = zext[0];
	mzranges[3] = zext[size[1] - 1];
	printf("MZ Range: %f to %f\n", mzranges[0], mzranges[1]);
	printf("Z Range: %f to %f\n", mzranges[2], mzranges[3]);

	// Set Up Peak Shapes
	peakshape = calloc(lines, sizeof(float));
	mkernel = calloc(lines, sizeof(float));
	if (peakshape == NULL || mkernel == NULL) {
		printf("Error allocating memory for peak shape and kernel\n");
		exit(1);
	}
	fftwf_complex* peakshape_FFT, *inverse_peakshape_FFT, *mkernel_FFT;
	peakshape_FFT = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * lines);
	inverse_peakshape_FFT = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * lines);
	mkernel_FFT = (fftwf_complex*)fftwf_malloc(sizeof(fftwf_complex) * lines);

	//Make Peak Shape Kernel
	MakeKernel2D(peakshape, size, mzext, zext, config.mzsig, config.csig, config.psfun, config.zpsfun);
	// Precompute FFTs
	precompute_fft2D(peakshape, size, peakshape_FFT);
	complex_conjugate(peakshape_FFT, inverse_peakshape_FFT, lines);

	printf("Peak Shape Set\n");


	//Set up charge blur
	zupind = calloc(lines, sizeof(int));
	zloind = calloc(lines, sizeof(int));
	if (zupind == NULL || zloind == NULL) {
		printf("Error allocating memory for charge blur indices\n");
		exit(1);
	}
	if (config.zsig != 0) { setup_blur_z(zupind, zloind, mzdat, zdat, lines, config.adductmass, mzranges, size);
	printf("Charge Blur Set\n");
	}

	//Set up mass blur
	mupind = calloc(lines, sizeof(int));
	mloind = calloc(lines, sizeof(int));
	if (mupind == NULL || mloind == NULL) {
		printf("Error allocating memory for mass blur indices\n");
		exit(1);
	}
	if (config.msig != 0) { setup_blur_m(mupind, mloind, mzdat, zdat, lines, config.adductmass, mzranges, size, config.molig);
	printf("Mass Blur Set: %f %f\n", config.molig, config.msig);
	}


	//Set up iteration
	blur = calloc(lines, sizeof(float));
	newblur = calloc(lines, sizeof(float));
	newblur2 = calloc(lines, sizeof(float));
	oldblur = calloc(lines, sizeof(float));
	barr = calloc(lines, sizeof(char));
	if (blur == NULL || newblur == NULL || newblur2 == NULL || oldblur == NULL || barr == NULL) {
		printf("Error allocating memory for iteration arrays\n");
		exit(1);
	}
	for (int i = 0; i < lines; i++) { barr[i] = 1; }
	size_t matsize = lines * sizeof(float);
	memcpy(blur, dataInt, matsize);
	memcpy(oldblur, blur, matsize);

	printf("Iterating: \n");
	//Iterating
	float conv=0;
	int off = 0;
	for (int m = 0; m < config.numit; m++) {
		// Apply softmax
		if (config.beta > 0) {
			softargmax(blur, size[0], size[1], config.beta);
		}
		// Apply point smoothing
		if (config.psig > 0) {
			point_smoothing(blur, barr, size[0], size[1], abs((int)config.psig));
		}
		// Apply charge smoothing
		if (config.zsig!=0) {
			blur_it_CD(newblur, blur, zupind, zloind, lines, config.zsig);
			memcpy(blur, newblur, matsize);
		}
		// Apply mass smoothing
		if (config.msig != 0) {
			blur_it_CD(newblur, blur, mupind, mloind, lines, config.msig);
			memcpy(blur, newblur, matsize);
		}

		// Richardson Lucy
		fftconvolve2D_precomputed(newblur, blur, peakshape_FFT, size, p1, p3, in1, out1);
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < lines; i++) {
			if(newblur[i]!=0){
				newblur2[i] = dataInt[i] / newblur[i];
			}
			else {
				newblur2[i] = 0;
			}
		}
		fftconvolve2D_precomputed(newblur, newblur2, inverse_peakshape_FFT, size, p1, p3, in1, out1);
		#pragma omp parallel for schedule(dynamic)
		for (int i = 0; i < lines; i++) {
			blur[i] = blur[i] * newblur[i];
		}

		//Determine the metrics for conversion. Only do this every 10% to speed up. Stop loop if converged
		if ((config.numit < 10 || m % 10 == 0 || m % 10 == 1 || m>0.9 * config.numit)) {
			float diff = 0;
			float tot = 0;
			for (int i = 0; i < lines; i++)
			{
				diff += powf((blur[i] - oldblur[i]), 2);
				tot += blur[i];
			}
			if (tot != 0) { conv = (diff / tot); }
			else {
				if (conv == 12345678) { printf("m/z vs. charge grid is zero. Iteration: %d\n", m); break; }
				else { conv = 12345678; }
			}

			if (conv < 0.000001) {
				if (off == 1 && config.numit > 0) {
					printf("Converged in %d iterations.\n", m);
					break;
				}
				off = 1;
			}
			memcpy(oldblur, blur, matsize);
		}

	}
	printf("Completed Iterations\n");
	//Writing outputs

	//Outputting Fit Reconvolved Data as newblur2
	memcpy(newblur, blur, matsize);
	fftconvolve2D_precomputed(newblur2, newblur, peakshape_FFT, size, p1, p3, in1, out1);
	//Normalize if necessary
	if (config.datanorm == 1) {
		float blurmax = Max(newblur2, lines);
		if (dmax != 0) {
			Normalize(lines, newblur2, blurmax / dmax);
		}
	}
	ApplyCutoff(newblur2, 0, lines);
	//Write the fitdat file
	char* suffixfit = "fitdat";
	write1D(config.outfile, suffixfit, newblur2, lines);


	//Writing Main Output

	//Reconvolve if necessary
	if (config.rawflag == 0) {
		MakeKernel2D(mkernel, size, mzext, zext, config.mzsig, 0, config.psfun, config.zpsfun);
		precompute_fft2D(mkernel, size, mkernel_FFT);

		memcpy(newblur, blur, matsize);
		fftconvolve2D_precomputed(blur, newblur, mkernel_FFT, size, p1, p3, in1, out1);
		printf("Reconvolved with m/z dimension\n");
	}

	//Normalize if necessary
	if (config.datanorm == 1) {
		float blurmax = Max(blur, lines);
		if (dmax != 0) {
			Normalize(lines, blur, blurmax / dmax);
		}
		printf("Maxes: %f %f\n",blurmax, dmax);
	}
	//Write the Output file
	ApplyCutoff(blur, 0, lines);

	// Write to binary
	//char* suffixfit = "";
	//write1D(config.outfile, suffixfit, blur, lines);

	// Write to text
	FILE* out_ptr = NULL;
	char outstring4[510];
	sprintf(outstring4, "%s_decon.txt", config.outfile);
	out_ptr = fopen(outstring4, "w");
	if (out_ptr == 0) { printf("Error Opening %s\n", outstring4); exit(1); }
	for (int i = 0; i < lines; i++)
	{
		fprintf(out_ptr, "%f\n", blur[i]);
	}
	fclose(out_ptr);
	printf("Wrote Deconvolution Output to: %s\n", outstring4);

	//Free memory
	free(mzdat);
	free(zdat);
	free(dataInt);

	free(blur);
	free(newblur);
	free(newblur2);
	free(oldblur);

	free(peakshape);
	fftwf_free(peakshape_FFT);
	free(mkernel);
	fftwf_free(mkernel_FFT);
	fftwf_free(inverse_peakshape_FFT);

	free(zloind);
	free(zupind);
	free(mloind);
	free(mupind);
	free(mzext);
	free(zext);

	fftwf_destroy_plan(p1);
	fftwf_destroy_plan(p3);
	fftwf_free(in1); fftwf_free(out1);

	endtime = time(NULL);
	printf("\nDone in %ds!\n", (int)difftime(endtime, starttime));
	return 0;
}