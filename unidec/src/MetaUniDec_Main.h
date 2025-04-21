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
//#include "UniDecLC_Main.h"
//#include "UD_peak_width.h"

int run_metaunidec(int argc, char* argv[], Config config) {
	clock_t starttime;
	starttime = clock();
	//Get Length
	int num = 0;
	config.file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
	num = int_attr(config.file_id, "/ms_dataset", "num", num);
	if (num > 20) { config.silent = 1; }

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
		else if (strcmp(argv[2], "-newgrids") == 0) { mode = 9; }
		else if (strcmp(argv[2], "-scanpeaks") == 0) { mode = 10; }
	}

	if (mode == 1 || mode == 2 || mode == 4 || mode==0)
	{
		//Iterate through files
		for (int i = 0; i < num; i++)
		{
			// Stop printing every file for large files
			if (num > 100) {
				if (i % (num / 10) == 0) {
					config.silent = 0;
				}
				else
				{
					config.silent = 1;
				}
			}
			// Run either deconvolution or processing
			config.metamode = i;
			if (mode == 1)
			{
				//printf("Deconvolving\n");
				run_unidec(argc, argv, config);
			}
			else if (mode == 2)
			{
				//printf("Processing\n");
				process_data(argc, argv, config);
			}
			else
			{
				process_data(argc, argv, config);
				run_unidec(argc, argv, config);
			}
		}
	}
	
	
	if (mode == 3 || mode == 4)
	{
		//printf("Making Merged Grids\n");
		if (question_grids(config.file_id)) { printf("Grids Already Made\n"); } //Checks to see if grids are already made
		else {
			make_grid(argc, argv, config, "/mass_data", "/mass_grid", "/mass_axis", "/mass_sum");
			make_grid(argc, argv, config, "/processed_data", "/mz_grid", "/mz_axis", "/mz_sum");
			set_got_grids(config.file_id);
		}
		
		get_peaks(argc, argv, config, 0);

		if (config.exchoice == 6)
		{
			printf("Charge Extraction Avg\n");
			config.exchoicez = 1;
			config.exnormz = config.exnorm;
			charge_peak_extracts(argc, argv, config, 0);
		}
		else if (config.exchoice == 7)
		{
			printf("Charge Extraction Max\n");
			config.exchoicez = 0;
			config.exnormz = config.exnorm;
			charge_peak_extracts(argc, argv, config, 0);
		}
		//get_peak_widths(argc, argv, config);
	}

	if (mode == 9)
	{
		// Forces new grids to be made. Mostly for testing
		printf("Making Merged Grids\n");
		make_grid(argc, argv, config, "/mass_data", "/mass_grid", "/mass_axis", "/mass_sum");
		make_grid(argc, argv, config, "/processed_data", "/mz_grid", "/mz_axis", "/mz_sum");
		set_got_grids(config.file_id);
		get_peaks(argc, argv, config, 0);
		//get_peak_widths(argc, argv, config);
	}
	
	if (mode == 5)
	{
		printf("Extracting Data\n");
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
		printf("Picking Peaks\n");
		get_all_peaks(argc, argv, config);
	}

	if (mode == 10)
	{
		printf("Getting Scan Scores\n");
		get_scan_scores(argc, argv, config);
	}
	
	// Close the file
	H5Fclose(config.file_id);
	
	clock_t end = clock();
	float totaltime = (float)(end - starttime) / CLOCKS_PER_SEC;
	printf("Done in %f s\n", totaltime);
	return 0;
}