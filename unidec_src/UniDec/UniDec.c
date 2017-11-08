// Universal Deconvolution of Mass Spectra
//
// Written by Michael Marty with Erik Marklund and Andrew Baldwin
//
// Copyright University of Oxford 2015
// Copyright University of Arizona 2016
//


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "UniDec_Main.h"
#include "UniDecIM_Main.h"
#include "MetaUniDec_Main.h"



Config ImportConfig(int argc, char * argv[], Config config)
{
	//Check for config file.
	if (argc >= 2)
	{
		if (strstr(argv[1], ".hdf5")) {
			hid_t file_id;
			config.filetype = 1;
			file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
			config = mh5LoadConfig(config, file_id);
			printf("Using HDF5 mode\n");
			H5Fclose(file_id);
		}
		else {
			config.filetype = 0;
			config = LoadConfig(config, argv[1]);
		}
	}
	else {
		//If no config file, print help page
		PrintHelp();
		exit(88);
	}
	return config;
}

int main(int argc, char *argv[])
{
	
	int result = 0;

	if (argc > 2)
	{
		if (strcmp(argv[2], "-nthreads") == 0) { omp_set_num_threads(atoi(argv[3])); printf("Number of Threads: %d\n", atoi(argv[3])); }
	}

	Config config = SetDefaultConfig();

	config=ImportConfig(argc, argv, config);

	if (config.metamode != -2)
	{
		printf("MetaUniDec Run: %d\n", config.metamode);
		result = run_metaunidec(argc, argv, config);
		return result;
	}


	if (config.imflag == 1)
	{
		printf("UniDec Ion Mobility Run\n");
		result = run_unidec_IM(argc, argv, config);
	}
	else
	{
		printf("UniDec Run\n");
		result = run_unidec(argc, argv, config);
	}

	return result;

}
//
