//
// Created by mm96978 on 7/9/2025.
//

#include "UD_conv.h"

void create_sim_spec(Config config, Input inp, const float* blur, const float mass, float* simint)
{
	for (int i = 0; i < config.numz; i++)
	{
		// Calculate the mz value for the mass and charge state
		float mz = calcmz(mass, inp.nztab[i], config.adductmass);
		// Find the nearest point in the original data
		int nearindex = nearfast(inp.dataMZ, mz, config.lengthmz);
		float nearmz = inp.dataMZ[nearindex];
		// If the nearest point is close enough, set the simulated intensity to the decon intensity at that index and z
		if (fabsf(mz - nearmz) < config.psmzthresh)
		{
			float deconint = blur[index2D(config.numz, nearindex, i)];
			simint[nearindex] = deconint;
		}
		//printf("%f %d %f\n", mz, nearindex, nearmz);
	}
}

// Convolution main function
int conv_main(int argc, char* argv[], Config config)
{
	printf("Running Convolution\n");
	Input inp = SetupInputs();
	ReadInputs(argc, argv, &config, &inp);
	CalcMasses(&config, &inp);

	char pfile[510];
	sprintf(pfile, "%s_peaks.dat", config.outfile);
	int plen = getfilelength(pfile);
	float* peakmasses, * peakints;
	peakmasses = calloc(plen, sizeof(float));
	peakints = calloc(plen, sizeof(float));
	if (peakmasses == NULL || peakints == NULL) {
		printf("Error allocating memory for peak masses and intensities.\n");
		return 1;
	}
	readfile(pfile, plen, peakmasses, peakints);
	printf("Pfile Length: %d\n", plen);

	char gfile[510];
	sprintf(gfile, "%s_grid.bin", config.outfile);
	int glen = getfilelengthbin(gfile, sizeof(float), 1);
	float* blur;
	blur = calloc(glen, sizeof(float));
	if (blur == NULL) {
		printf("Error allocating memory for blur array.\n");
		free(peakmasses);
		free(peakints);
		return 1;
	}
	printf("Gfile Length: %d\n", glen);

	readfile1bin(gfile, glen, blur);


	//...................................................................
	//
	//     Sets the mzdist with the peak shape
	//
	//....................................................................

	Decon decon = SetupDecon();
	float* outarray = calloc(config.lengthmz * config.numz, sizeof(float));
	if (outarray == NULL) {
		printf("Error allocating memory for output array.\n");
		free(peakmasses);
		free(peakints);
		free(blur);
		return 1;
	}

	int maxlength = SetUpPeakShape(config, inp, &decon, 0, 1);

	float* outptr = outarray;

	for (int i = 0; i < plen; i++)
	{
		float mass = peakmasses[i];
		float height = peakints[i];
		printf("%f %f\n", mass, height);

		float* simint, * convint;
		simint = calloc(config.lengthmz, sizeof(float));
		convint = calloc(config.lengthmz, sizeof(float));
		if (simint == NULL || convint == NULL) {
			printf("Error allocating memory for simint or convint.\n");
			free(outarray);
			free(peakmasses);
			free(peakints);
			free(blur);
			return 1;
		}

		create_sim_spec(config, inp, blur, mass, simint);

		convolve_simp(config.lengthmz, maxlength, decon.starttab, decon.endtab, decon.mzdist, simint, convint, config.speedyflag);

		memcpy(outptr, convint, config.lengthmz * sizeof(float));

		free(simint);
		free(convint);
		outptr += config.lengthmz;
	}

	char gfileout[510];
	sprintf(gfileout, "%s_conv.bin", config.outfile);
	writefile1bin(gfileout, config.lengthmz * plen, outarray);

	FreeDecon(decon);
	free(outarray);

	FreeInputs(inp);
	free(blur);
	free(peakints);
	free(peakmasses);

	return 0;
}