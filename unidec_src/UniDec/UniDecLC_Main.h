/*
* UniDecLC_Main.h
*
*  Created on : 23 December 2020
* Author : Michael.Marty
*/

//
// Copyright 2020 University of Arizona
//
//

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <fftw3.h>

int run_unidec_LC(int argc, char* argv[], Config config) {
	//Initialize
	clock_t starttime;
	starttime = clock();

	FILE* out_ptr = NULL;
	hid_t file_id;

	char dataset[1024];
	char outdat[1024];
	strcpy(dataset, "/ms_dataset");
	char strval[1024];
	sprintf(strval, "/%d", config.metamode);
	strcat(dataset, strval);
	printf("HDF5 Data Set: %s\n", dataset);

	return 0;
}
