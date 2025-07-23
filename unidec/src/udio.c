//
// Created by mm96978 on 7/7/2025.
//

#include "udio.h"


//Reads in x y file.
void readfile(const char *infile, const int lengthmz, float *dataMZ, float *dataInt) {
    char x[500] = {0};
    char y[500] = {0};

    FILE * file_ptr = fopen(infile, "r");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }


    for (int i = 0; i < lengthmz; i++) {
        const int match = fscanf(file_ptr, "%s %s", x, y);
        if (match != 0) {
            // dataMZ[i] = atof(x);
            // dataInt[i] = atof(y);
            char * endptr, * endptr2;
            dataMZ[i] = strtof(x, &endptr);
            dataInt[i] = strtof(y, &endptr2);
            if (*endptr != '\0' || *endptr2 != '\0') {
                printf("Error converting string to float: %s %s\n", x, y);
                exit(2);
            }
        }
    }

    fclose(file_ptr);
}

//Reads in x y z file.
void readfile3(char *infile, int lengthmz, float *array1, float *array2, float *array3) {
    char x[500] = {0};
    char y[500] = {0};
    char z[500] = {0};

    FILE * file_ptr = fopen(infile, "r");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    for (int i = 0; i < lengthmz; i++) {
        const int match = fscanf(file_ptr, "%s %s %s", x, y, z);
        if (match != 0) {
            // array1[i] = atof(x);
            // array2[i] = atof(y);
            // array3[i] = atof(z);
            char * endptr, * endptr2, * endptr3;
            array1[i] = strtof(x, &endptr);
            array2[i] = strtof(y, &endptr2);
            array3[i] = strtof(z, &endptr3);
            if (*endptr != '\0' || *endptr2 != '\0' || *endptr3 != '\0') {
                printf("Error converting string to float: %s %s %s\n", x, y, z);
                exit(2);
            }
        }
    }
    fclose(file_ptr);
}

//Reads in x y z file.
void readfile3bin(const char *infile, const int lengthmz, float *array1, float *array2, float *array3) {
    const int l = lengthmz * 3;
    float * data = (float *) calloc(l, sizeof(float));

    FILE * file_ptr = fopen(infile, "rb");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    if (data) {
        const size_t result = fread(data, sizeof(float), l, file_ptr);
    	if (result != l) {
			printf("Error reading file %s: expected %d floats, got %zu\n", infile, l, result);
			free(data);
			fclose(file_ptr);
			exit(1);
		}
        for (int i = 0; i < lengthmz; i++) {
            array1[i] = data[i * 3];
            array2[i] = data[i * 3 + 1];
            array3[i] = data[i * 3 + 2];
            //printf("%f %f %f\n", array1[i], array2[i], array3[i]);
        }
        free(data);
    }
    fclose(file_ptr);
}

//Reads in list binary file.
void readfile1bin(const char *infile, const int lengthmz, float *data) {
    FILE *file_ptr = fopen(infile, "rb");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }
	const size_t result = fread(data, sizeof(float), lengthmz, file_ptr);
	if (result != lengthmz) {
		printf("Error reading file %s: expected %d floats, got %zu\n", infile, lengthmz, result);
		fclose(file_ptr);
		exit(1);
	}
    fclose(file_ptr);
}

//Writes to binary file
void writefile1bin(const char *outstring, const int lengthmz, const float *data) {
    FILE *out_ptr = fopen(outstring, "wb");
    if (out_ptr == 0) {
        printf("Error Opening %s \n", outstring);
        exit(1);
    }
    fwrite(data, sizeof(float), lengthmz, out_ptr);
    fclose(out_ptr);
}


void write1D(char *outfile, char *suffix, const float *array, const int length)
{
	char outstring[500];
	sprintf(outstring, "%s_%s.bin", outfile, suffix);
	writefile1bin(outstring, length, array);
	printf("File written to: %s\n", outstring);
}

void write2D(char *outfile, char *suffix, const float *array1, const float *array2, const int length)
{
	char outstring[500];
	sprintf(outstring, "%s_%s.txt", outfile, suffix);
	FILE *out_ptr = fopen(outstring, "w");
	if (out_ptr == 0) { printf("Error Opening %s\n", outstring); exit(1); }
	for (int i = 0; i<length; i++)
	{
		fprintf(out_ptr, "%f %f\n", array1[i], array2[i]);
	}
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}

void write3D(char *outfile, char *suffix, const float *array1, const float *array2, const float *array3, const int length)
{
	char outstring[500];
	sprintf(outstring, "%s_%s.txt", outfile, suffix);
	FILE *out_ptr = fopen(outstring, "w");
	if (out_ptr == 0) { printf("Error Opening %s\n", outstring); exit(1); }
	for (int i = 0; i<length; i++)
	{
		fprintf(out_ptr, "%f %f %f\n", array1[i], array2[i], array3[i]);
	}
	fclose(out_ptr);
	printf("File written to: %s\n", outstring);
}


//Reads in single list of values
void readmfile(const char *infile, const int mfilelen, float *testmasses) {
    char input[500] = {0};

    FILE * file_ptr = fopen(infile, "r");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    for (int i = 0; i < mfilelen; i++) {
        const int match = fscanf(file_ptr, "%s", input);
        if (match != 0) {
            // testmasses[i] = atof(input);
            char * endptr;
            testmasses[i] = strtof(input, &endptr);
            if (*endptr != '\0') {
                printf("Error converting string to float: %s\n", input);
                exit(2);
            }
        }
    }

    fclose(file_ptr);
}

// count the number of lines we have in datafile
int getfilelength(const char *infile) {
    int l = 0;
    FILE * file_ptr = fopen(infile, "r");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    int read = -2;
    while (read != EOF) {
        char input[501];
        read = fscanf(file_ptr, "%500[^\n]%*c", input);
        if (read != 0 && read != EOF) {
            l += 1;
        }
    }

    fclose(file_ptr);
    return l;
}

// count the number of lines we have in datafile
int getfilelengthbin(const char *infile, const int size, const int width) {
    int l = 0;
    FILE * file_ptr = fopen(infile, "rb");
    if (file_ptr == 0) {
        printf("Error Opening %s\n", infile);
        exit(1);
    }

    fseek(file_ptr, 0, SEEK_END);
    l = ftell(file_ptr);
    l /= (size * width);

    fclose(file_ptr);
    return l;
}


void ReadInputs(int argc, char* argv[], Config* config, Input* inp)
{

	if (config->filetype == 1) {
		if (config->metamode != -2)
		{
			strcpy(config->dataset, "/ms_dataset");
			char strval[1024];
			sprintf(strval, "/%d", config->metamode);
			strcat(config->dataset, strval);
			if (config->silent == 0) { printf("HDF5: %s\n", config->dataset); }
		}
		else
		{
			strcpy(config->dataset, "/ms_data");
		}

		char outdat[1024];
		strjoin(config->dataset, "/processed_data", outdat);
		//config->file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);

		config->lengthmz = mh5getfilelength(config->file_id, outdat);
		inp->dataMZ = (float*)calloc(config->lengthmz, sizeof(float));
		inp->dataInt = (float*)calloc(config->lengthmz, sizeof(float));
		if (inp->dataMZ == NULL || inp->dataInt == NULL) {
			printf("Error allocating memory for data arrays\n");
			exit(1);
		}
		mh5readfile2d(config->file_id, outdat, config->lengthmz, inp->dataMZ, inp->dataInt);
		if (config->silent == 0) { printf("Length of Data: %d \n", config->lengthmz); }

		//Check the length of the mfile and then read it in.
		if (config->mflag == 1)
		{
			config->mfilelen = mh5getfilelength(config->file_id, "/config/masslist");
			printf("Length of mfile: %d \n", config->mfilelen);
			inp->testmasses = (float*)malloc(sizeof(float) * config->mfilelen);
			mh5readfile1d(config->file_id, "/config/masslist", inp->testmasses);
		}
		else {
			inp->testmasses = (float*)malloc(sizeof(float) * config->mfilelen);
		}

	}
	else {

		//Calculate the length of the data file automatically
		config->lengthmz = getfilelength(config->infile);
		inp->dataMZ = (float*)calloc(config->lengthmz, sizeof(float));
		inp->dataInt = (float*)calloc(config->lengthmz, sizeof(float));
		if (inp->dataMZ == NULL || inp->dataInt == NULL) {
			printf("Error allocating memory for data arrays\n");
			exit(1);
		}

		readfile(config->infile, config->lengthmz, inp->dataMZ, inp->dataInt);//load up the data array
		printf("Length of Data: %d \n", config->lengthmz);

		//Check the length of the mfile and then read it in.

		if (config->mflag == 1)
		{
			config->mfilelen = getfilelength(config->mfile);
			printf("Length of mfile: %d \n", config->mfilelen);
		}
		inp->testmasses = (float*)malloc(sizeof(float) * config->mfilelen);
		if (inp->testmasses == NULL) {
			printf("Error allocating memory for test masses\n");
			exit(1);
		}
		if (config->mflag == 1)
		{
			readmfile(config->mfile, config->mfilelen, inp->testmasses);//read in mass tab
		}
	}

	SetupZtab(*config, inp);
}

void SetupZtab(const Config config, Input *inp) {
	//This for loop creates a list of charge values
	inp->nztab = (int*)calloc(config.numz, sizeof(int));

	if (inp->nztab) {
		for (int i = 0; i < config.numz; i++) { inp->nztab[i] = i + config.startz; }

		//printf("nzstart %d\n",inp.nztab[0]);

		//Test to make sure no charge state is zero
		for (int j = 0; j < config.numz; j++)
		{
			if (inp->nztab[j] == 0) { printf("Error: Charge state cannot be 0"); exit(100); }
		}
		//Test to make sure no two data points has the same x value
		for (int i = 0; i < config.lengthmz - 1; i++)
		{
			if (inp->dataMZ[i] == inp->dataMZ[i + 1]) { printf("Error: Two data points are identical: %f %f \n\n", inp->dataMZ[i], inp->dataMZ[i + 1]); exit(104); }
		}
	}
}

void WritePeaks(const Config config, const Decon* decon) {
	char outdat[1024];
	strjoin(config.dataset, "/peaks", outdat);
	float* ptemp = NULL;
	int l = decon->plen * 4;
	ptemp = (float*)calloc(l, sizeof(float));
	if (ptemp) {
		for (int i = 0; i < decon->plen; i++) {
			ptemp[i * 3] = decon->peakx[i];
			ptemp[i * 3 + 1] = decon->peaky[i];
			ptemp[i * 3 + 2] = decon->dscores[i];
			ptemp[i * 3 + 3] = (float) config.metamode;
		}

		mh5writefile2d_grid(config.file_id, outdat, decon->plen, 4, ptemp);
		free(ptemp);
	}
}




void WriteDecon(const Config config, const Decon* decon, const Input* inp)
{
	hid_t file_id = config.file_id;

	char outdat[1024];
	FILE* out_ptr = NULL;
	//Write the fit data to a file.
	if (config.rawflag >= 0) {
		if (config.filetype == 0) {
			char outstring2[511];
			sprintf(outstring2, "%s_fitdat.bin", config.outfile);
			out_ptr = fopen(outstring2, "wb");
			if (out_ptr == 0) { printf("Error Opening %s \n", outstring2); exit(1); }
			fwrite(decon->fitdat, sizeof(float), config.lengthmz, out_ptr);
			fclose(out_ptr);
			//printf("Fit: %s\t", outstring2);
		}
		else {
			//strjoin(dataset, "/fit_data", outdat);
			//mh5writefile1d(file_id, outdat, config.lengthmz, decon->fitdat);
		}
	}

	//Write the baseline to a file.
	if (config.rawflag >= 0 && config.baselineflag == 1) {
		if (config.filetype == 0) {
			char outstring2[513];
			sprintf(outstring2, "%s_baseline.bin", config.outfile);
			out_ptr = fopen(outstring2, "wb");
			if (out_ptr == 0) { printf("Error Opening %s \n", outstring2); exit(1); }
			fwrite(decon->baseline, sizeof(float), config.lengthmz, out_ptr);
			fclose(out_ptr);
			//printf("Background: %s\t", outstring2);
		}
		else {
			//strjoin(dataset, "/baseline", outdat);
			//mh5writefile1d(file_id, outdat, config.lengthmz, decon->baseline);//
		}
	}

	//Writes the convolved m/z grid in binary format
	if (config.rawflag == 0 || config.rawflag == 1)
	{
		// Note to self
		// rawflag=0 -> Reconvolved/Profile -> newblur
		// rawflag=1 -> Raw/Centroid -> blur
		int l = config.lengthmz * config.numz;
		if (config.filetype == 0) {
			char outstring9[509];
			sprintf(outstring9, "%s_grid.bin", config.outfile);
			out_ptr = fopen(outstring9, "wb");
			if (out_ptr == 0) { printf("Error Opening %s\n", outstring9); exit(1); }
			if (config.rawflag == 0) { fwrite(decon->newblur, sizeof(float), l, out_ptr); }
			if (config.rawflag == 1) { fwrite(decon->blur, sizeof(float), l, out_ptr); }
			fclose(out_ptr);
			//printf("m/z grid: %s\t", outstring9);
		}
		else {
			strjoin(config.dataset, "/mz_grid", outdat);
			if (config.rawflag == 0) { mh5writefile1d(file_id, outdat, l, decon->newblur); }
			if (config.rawflag == 1) { mh5writefile1d(file_id, outdat, l, decon->blur); }

			float* chargedat = NULL;
			chargedat = (float*)calloc(config.numz, sizeof(float));
			float* chargeaxis = NULL;
			chargeaxis = (float*)calloc(config.numz, sizeof(float));

			if (chargedat && chargeaxis) {

				for (int j = 0; j < config.numz; j++) {
					float val = 0;
					chargeaxis[j] = (float)inp->nztab[j];
					for (int i = 0; i < config.lengthmz; i++) {

						val += decon->newblur[index2D(config.numz, i, j)];
					}
					chargedat[j] = val;
				}
				strjoin(config.dataset, "/charge_data", outdat);
				mh5writefile2d(file_id, outdat, config.numz, chargeaxis, chargedat);
				free(chargedat);
				free(chargeaxis);
			}
		}
	}
	else if (config.filetype == 1) {
		// If the rawfile flag not 1 or 2 and this is an HDF5 file, zero out these arrays.
		strjoin(config.dataset, "/mz_grid", outdat);
		delete_group(file_id, outdat);
		strjoin(config.dataset, "/charge_data", outdat);
		delete_group(file_id, outdat);
		strjoin(config.dataset, "/mass_grid", outdat);
		delete_group(file_id, outdat);
	}

	//Writes the convolved mass grid in binary format
	if (config.rawflag == 0 || config.rawflag == 1) {
		int l = decon->mlen * config.numz;
		if (config.filetype == 0) {
			char outstring10[513];
			sprintf(outstring10, "%s_massgrid.bin", config.outfile);
			out_ptr = fopen(outstring10, "wb");
			if (out_ptr == 0) { printf("Error Opening %s\n", outstring10); exit(1); }
			fwrite(decon->massgrid, sizeof(float), l, out_ptr);
			fclose(out_ptr);
			//printf("Mass Grid: %s\t", outstring10);
		}
		else {
			strjoin(config.dataset, "/mass_grid", outdat);
			mh5writefile1d(file_id, outdat, l, decon->massgrid);
		}
	}

	//Writes the mass values convolved with the peak shape
	if (config.rawflag == 0 || config.rawflag == 1 || config.rawflag == 2 || config.rawflag == 3) {
		if (config.filetype == 0) {
			char outstring4[509];
			sprintf(outstring4, "%s_mass.txt", config.outfile);
			out_ptr = fopen(outstring4, "w");
			if (out_ptr == 0) { printf("Error Opening %s\n", outstring4); exit(1); }
			for (int i = 0; i < decon->mlen; i++)
			{
				fprintf(out_ptr, "%f %f\n", decon->massaxis[i], decon->massaxisval[i]);
			}
			fclose(out_ptr);
			//printf("Masses: %s\n", outstring4);
		}
		else {
			//int mlen = remove_middle_zeros(decon->massaxis, decon->massaxisval, decon->mlen); //Doesn't really work
			strjoin(config.dataset, "/mass_data", outdat);
			mh5writefile2d(file_id, outdat, decon->mlen, decon->massaxis, decon->massaxisval);
		}
	}

	if (config.filetype == 1 && decon->plen > 0) {
		WritePeaks(config, decon);
	}

}


void ManualAssign(const float* dataMZ, char* barr, const int* nztab, Config config)
{
	int manlen = 0;
	int lengthmz = config.lengthmz;
	int numz = config.numz;
	if (config.filetype == 1)
	{
		manlen = mh5getfilelength(config.file_id, "/config/manuallist");
	}
	else
	{
		manlen = getfilelength(config.manualfile);
	}
	printf("Length of Manual List: %d \n", manlen);
	float* manualmz = (float*)malloc(sizeof(float) * manlen);
	float* manualwin = (float*)malloc(sizeof(float) * manlen);
	float* manualassign = (float*)malloc(sizeof(float) * manlen);
	if (manualmz == NULL || manualwin == NULL || manualassign == NULL) {
		printf("Error allocating memory for manual assignment arrays\n");
		exit(1);
	}
	if (config.filetype == 1)
	{
		mh5readfile3d(config.file_id, "/config/manuallist", manlen, manualmz, manualwin, manualassign);
	}
	else
	{
		readfile3(config.manualfile, manlen, manualmz, manualwin, manualassign);
	}

	for (int i = 0; i < manlen; i++) {
		if (manualassign[i] < 0) { manualmz[i] = -manualmz[i] * 1000.0f; }
		printf("Manual Assignment: %f %f %f\n", manualmz[i], manualwin[i], manualassign[i]);
		//Cheating a bit...make the manualmz very negative if deassign, so that they aren't discovered by the assign loop
	}

	#pragma omp parallel for schedule(auto)
	for (int i = 0; i < lengthmz; i++)
	{
		for (int j = 0; j < numz; j++)
		{
			//Manual Assign: The nearest value wins in the case of overlap
			float testmz = dataMZ[i];
			int closest = nearunsorted(manualmz, testmz, manlen);
			if (fabsf(manualmz[closest] - testmz) < manualwin[closest] && (float)nztab[j] == manualassign[closest] && manualassign[closest] > 0)
			{
				barr[index2D(numz, i, j)] = 1;
			}
			else if (fabsf(manualmz[closest] - testmz) < manualwin[closest] && (float)nztab[j] != manualassign[closest] && manualassign[closest] > 0)
			{
				barr[index2D(numz, i, j)] = 0;
			}
			//Manual Deassign: Anything within the window is killed
			for (int k = 0; k < manlen; k++) {

				if (fabsf(manualmz[k] - testmz * (-1000.0f)) < manualwin[k] * 1000.0f && (float)nztab[j] == fabsf(manualassign[k]) && manualassign[k] < 0)
				{
					barr[index2D(numz, i, j)] = 0;
				}

			}
		}
	}
	free(manualmz);
	free(manualwin);
	free(manualassign);
	printf("Using Manual Assignments for Some Peaks\n");
}
