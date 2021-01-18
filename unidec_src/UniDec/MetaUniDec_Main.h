/*
* MetaUniDec_Main.h
*
*  Created on : 16 May 2017
* Author : Michael.Marty
*/

//
// 
// Copyright 2017 University of Arizona
//
//
#include "UD_dataproc.h"
#include "UD_analysis.h"
#include "UD_charge.h"
#include "UniDecLC_Main.h"
//#include "UD_peak_width.h"

int run_metaunidec(int argc, char *argv[], Config config) {
	//Get Length
	int num=0;
	hid_t file_id;
	file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
	num = int_attr(file_id, "/ms_dataset", "num", num);
	H5Fclose(file_id);

	int mode = 0;
	if (argc > 2)
	{
		if (strcmp(argv[2], "-decon") == 0) { mode = 1; }
		else if (strcmp(argv[2], "-proc") == 0) { mode = 2; }
		else if (strcmp(argv[2], "-grids") == 0) { mode = 3; }
		else if (strcmp(argv[2], "-all") == 0) { mode = 4; }
		else if (strcmp(argv[2], "-extract") == 0) { mode = 5; }
		else if (strcmp(argv[2], "-ultraextract") == 0) { mode = 6; }
		else if (strcmp(argv[2], "-charges") == 0) { mode = 7; }
		else if (strcmp(argv[2], "-peaks") == 0) { mode = 8; }
		else if (strcmp(argv[2], "-LC") == 0) { mode = 9; }
		else if (strcmp(argv[2], "-scanpeaks") == 0) { mode = 10; }
	}

	//printf("%d\n\n", mode);
	//Iterate through files
	for (int i = 0; i < num; i++)
	{
		config.metamode = i;
		if(mode==1)
		{			
			//printf("Deconvolving\n");
			run_unidec(argc, argv, config);
		}
		else if (mode == 2)
		{
			//printf("Processing\n");
			process_data(argc, argv, config);
		}
		else if (mode == 3 || mode >= 5)
		{

		}
		else
		{
			process_data(argc, argv, config);
			run_unidec(argc, argv, config);
		}
	}
	
	if ( mode ==3 || mode == 4)
	{
		//printf("Making Merged Grids\n");
		hid_t file_id;
		file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
		if (question_grids(file_id)) { printf("Grids Already Made\n"); H5Fclose(file_id); } //Checks to see if grids are already made
		else{
			H5Fclose(file_id);
			make_grid(argc, argv, config, "/mass_data", "/mass_grid", "/mass_axis", "/mass_sum");
			make_grid(argc, argv, config, "/processed_data", "/mz_grid", "/mz_axis", "/mz_sum");
			file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
			set_got_grids(file_id);
			H5Fclose(file_id);
		}
		get_peaks(argc, argv, config, 0);
		//get_peak_widths(argc, argv, config);
	}

	if (mode == 5)
	{
		//printf("Extracting Data\n");
		get_peaks(argc, argv, config, 0);
	}

	if (mode == 6)
	{
		//clock_t starttime = clock();
		printf("Extracting Data ULTRA\n");
		get_peaks(argc, argv, config, 1);
		charge_peak_extracts(argc, argv, config, 1);
		//clock_t end = clock();
		//float totaltime = (float)(end - starttime) / CLOCKS_PER_SEC;
		//printf("\nFinished in %f seconds\n",totaltime);
	}

	if (mode == 7)
	{
		printf("Extracting Charges\n");
		charge_peak_extracts(argc, argv, config, 0);
	}
	if (mode == 8)
	{
		//printf("Picking Peaks\n");
		get_all_peaks(argc, argv, config);
	}

	if (mode == 9)
	{
		printf("LC Mode. Under Construction...");
		// Merge Data
		make_grid(argc, argv, config, "/raw_data", "/raw_grid", "/raw_axis", "/raw_sum");

		// Process Data
		//process_data(argc, argv, config);

		//Run UniDec LC
		run_unidec_LC(argc, argv, config);
	}

	if (mode == 10)
	{
		printf("Getting Scan Scores\n");
		get_scan_scores(argc, argv, config);
	}

	return 0;
}