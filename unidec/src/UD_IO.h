/*
* UD_IO.h
*
* 
* Author : Michael.Marty
*/

//
//
// Copyright 2023 University of Arizona
//
//

#pragma once

#ifndef UD_IO
#define UD_IO

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void mh5readfile3d(hid_t file_id, char* dataname, int lengthmz, float* dataMZ, float* dataInt, float* data3);
void mh5readfile2d(hid_t file_id, char* dataname, int lengthmz, float* dataMZ, float* dataInt);
void mh5readfile1d(hid_t file_id, char* dataname, float* data);
int mh5getfilelength(hid_t file_id, char* dataname);
void strjoin(const char* s1, const char* s2, char* newstring);
void mh5writefile1d(hid_t file_id, char* dataname, int length, float* data1);
void mh5writefile2d(hid_t file_id, char* dataname, int length, float* data1, float* data2);
void mh5writefile2d_grid(hid_t file_id, char* dataname, int length1, int length2, float* data1);
void delete_group(hid_t file_id, char* dataname);
float calcDtSLIMpoly2(float mass, int z, float ccs, float tcal1, float tcal2, float tcal3, float hmass, float edc);
float calcCCSSLIMpoly2(float mass, int z, float dt, float tcal1, float tcal2, float tcal3, float hmass, float edc);
float calcCCSSLIMpoly3(float mass, int z, float dt, float tcal1, float tcal2, float tcal3, float tcal4, float hmass, float edc);
float calcDtSLIMpoly3(float mass, int z, float ccs, float tcal1, float tcal2, float tcal3, float tcal4, float hmass, float edc);



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

		readfile(config->infile, config->lengthmz, inp->dataMZ, inp->dataInt);//load up the data array
		printf("Length of Data: %d \n", config->lengthmz);

		//Check the length of the mfile and then read it in.

		if (config->mflag == 1)
		{
			config->mfilelen = getfilelength(config->mfile);
			printf("Length of mfile: %d \n", config->mfilelen);
		}
		inp->testmasses = (float*)malloc(sizeof(float) * config->mfilelen);
		if (config->mflag == 1)
		{
			readmfile(config->mfile, config->mfilelen, inp->testmasses);//read in mass tab
		}
	}

	//This for loop creates a list of charge values
	inp->nztab = (int*)calloc(config->numz, sizeof(int));

	if (inp->nztab) {
		for (int i = 0; i < config->numz; i++) { inp->nztab[i] = i + config->startz; }

		//printf("nzstart %d\n",inp.nztab[0]);

		//Test to make sure no charge state is zero
		for (int j = 0; j < config->numz; j++)
		{
			if (inp->nztab[j] == 0) { printf("Error: Charge state cannot be 0"); exit(100); }
		}
		//Test to make sure no two data points has the same x value
		for (int i = 0; i < config->lengthmz - 1; i++)
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

/*
void ReadPeaks(const Config config, const char * dataset, float* peakx, float* peaky, float* dscores)
{
	//Unfinished
	char outdat[1024];
	strjoin(dataset, "/peaks", outdat);
	float* ptemp = NULL;
	//ptemp = calloc(decon->plen * 3, sizeof(float));
}*/


void WriteDecon(const Config config, const Decon* decon, const Input* inp)
{
	hid_t file_id = config.file_id;

	char outdat[1024];
	FILE* out_ptr = NULL;
	//Write the fit data to a file.
	if (config.rawflag >= 0) {
		if (config.filetype == 0) {
			char outstring2[500];
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
			char outstring2[500];
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
			char outstring9[500];
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
			char outstring10[500];
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
			char outstring4[500];
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





void ManualAssign(float* dataMZ, char* barr, int* nztab, Config config)
{
	unsigned int i, j, k;
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
	if (config.filetype == 1)
	{
		mh5readfile3d(config.file_id, "/config/manuallist", manlen, manualmz, manualwin, manualassign);
	}
	else
	{
		readfile3(config.manualfile, manlen, manualmz, manualwin, manualassign);
	}

	for (i = 0; i < manlen; i++) {
		if (manualassign[i] < 0) { manualmz[i] = -manualmz[i] * 1000.0; }
		printf("Manual Assignment: %f %f %f\n", manualmz[i], manualwin[i], manualassign[i]);
		//Cheating a bit...make the manualmz very negative if deassign, so that they aren't discovered by the assign loop
	}

	float testmz;
	int closest;
#pragma omp parallel for private (i,j,testmz,closest), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		for (j = 0; j < numz; j++)
		{
			//Manual Assign: The nearest value wins in the case of overlap
			testmz = dataMZ[i];
			closest = nearunsorted(manualmz, testmz, manlen);
			if (fabs(manualmz[closest] - testmz) < manualwin[closest] && (float)nztab[j] == manualassign[closest] && manualassign[closest] > 0)
			{
				barr[index2D(numz, i, j)] = 1;
			}
			else if (fabs(manualmz[closest] - testmz) < manualwin[closest] && (float)nztab[j] != manualassign[closest] && manualassign[closest] > 0)
			{
				barr[index2D(numz, i, j)] = 0;
			}
			//Manual Deassign: Anything within the window is killed
			for (k = 0; k < manlen; k++) {

				if (fabs(manualmz[k] - testmz * (-1000.0)) < manualwin[k] * 1000.0 && (float)nztab[j] == fabs(manualassign[k]) && manualassign[k] < 0)
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


#endif
