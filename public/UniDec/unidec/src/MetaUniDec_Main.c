//
// Created by mm96978 on 7/11/2025.
//

#include "MetaUniDec_Main.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define META_PARALLEL_BATCH_MAX 32

static int min_int(const int a, const int b) { return a < b ? a : b; }

static int can_parallel_fast_decon(const Config config, const int num) {
#ifdef _OPENMP
	return num > 1 && config.rawflag > 1 && config.manualflag == 0 && config.doubledec == 0 &&
		omp_get_max_threads() > 1;
#else
	return 0;
#endif
}

static int run_parallel_fast_decon(Config config, const int num, const int print_spectra) {
#ifndef _OPENMP
	return 0;
#else
	if (!can_parallel_fast_decon(config, num)) { return 0; }
	const int max_workers = min_int(num, min_int(omp_get_max_threads(), META_PARALLEL_BATCH_MAX));
	if (max_workers < 2) { return 0; }
	const int batch_total = (num + max_workers - 1) / max_workers;
	const int batch_size = (num + batch_total - 1) / batch_total;
	const int previous_max_active_levels = omp_get_max_active_levels();
	const int previous_dynamic = omp_get_dynamic();
	omp_set_max_active_levels(1);
	omp_set_dynamic(0);

	for (int batch_start = 0; batch_start < num; batch_start += batch_size) {
		const int batch_count = min_int(batch_size, num - batch_start);
		Config* configs = calloc(batch_count, sizeof(Config));
		Input* inputs = calloc(batch_count, sizeof(Input));
		Decon* decons = calloc(batch_count, sizeof(Decon));
		double* elapsed = calloc(batch_count, sizeof(double));
		if (configs == NULL || inputs == NULL || decons == NULL || elapsed == NULL) {
			fprintf(stderr, "Unable to allocate a MetaUniDec parallel batch.\n");
			free(configs);
			free(inputs);
			free(decons);
			free(elapsed);
			omp_set_max_active_levels(previous_max_active_levels);
			omp_set_dynamic(previous_dynamic);
			exit(11);
		}

		// HDF5 reads stay on the calling thread.
		for (int j = 0; j < batch_count; j++) {
			const int spectrum = batch_start + j;
			if (print_spectra) { printf("Spectrum %d/%d\n", spectrum + 1, num); }
			configs[j] = config;
			configs[j].silent = 1;
			configs[j].metamode = spectrum;
			inputs[j] = InitInputs();
			decons[j] = InitDecon();
			ReadInputs(&configs[j], &inputs[j]);
		}

		if (batch_count >= 4) {
			int* states = calloc(batch_count, sizeof(int));
			if (states == NULL) {
				fprintf(stderr, "Unable to allocate the MetaUniDec writer queue.\n");
				omp_set_max_active_levels(previous_max_active_levels);
				omp_set_dynamic(previous_dynamic);
				exit(11);
			}
			int next_job = 0;

			// Thread zero is the only HDF5 writer. Other threads run independent
			// deconvolution cores and publish completed spectra to the writer.
			#pragma omp parallel num_threads(batch_count) shared(next_job, states)
			{
				if (omp_get_num_threads() == 1) {
					// Resource-constrained runtimes may provide a one-thread team even
					// with dynamic teams disabled. Preserve progress in that case.
					for (int j = 0; j < batch_count; j++) {
						const double start = omp_get_wtime();
						run_unidec_core(configs[j], inputs[j], &decons[j], 0);
						elapsed[j] = omp_get_wtime() - start;
						WriteDecon(configs[j], &decons[j], &inputs[j]);
						WriteGlobalOutputs(configs[j], decons[j], (float)elapsed[j]);
						FreeInputs(inputs[j]);
						FreeDecon(decons[j]);
					}
				}
				else if (omp_get_thread_num() == 0) {
					int written = 0;
					int scan_start = 0;
					while (written < batch_count) {
						int found = 0;
						for (int offset = 0; offset < batch_count; offset++) {
							const int j = (scan_start + offset) % batch_count;
							int state = 0;
							#pragma omp atomic read
							state = states[j];
							if (state == 1) {
								WriteDecon(configs[j], &decons[j], &inputs[j]);
								WriteGlobalOutputs(configs[j], decons[j], (float)elapsed[j]);
								FreeInputs(inputs[j]);
								FreeDecon(decons[j]);
								#pragma omp atomic write
								states[j] = 2;
								written++;
								scan_start = (j + 1) % batch_count;
								found = 1;
								break;
							}
						}
						if (!found) {
							#pragma omp taskyield
						}
					}
				}
				else {
					while (1) {
						int j = 0;
						#pragma omp atomic capture
						j = next_job++;
						if (j >= batch_count) { break; }
						const double start = omp_get_wtime();
						run_unidec_core(configs[j], inputs[j], &decons[j], 0);
						elapsed[j] = omp_get_wtime() - start;
						#pragma omp atomic write
						states[j] = 1;
					}
				}
			}
			free(states);
		}
		else {
			// Small tail batches are faster without reserving a writer thread.
			#pragma omp parallel for schedule(dynamic) num_threads(batch_count)
			for (int j = 0; j < batch_count; j++) {
				const double start = omp_get_wtime();
				run_unidec_core(configs[j], inputs[j], &decons[j], 0);
				elapsed[j] = omp_get_wtime() - start;
			}
			for (int j = 0; j < batch_count; j++) {
				WriteDecon(configs[j], &decons[j], &inputs[j]);
				WriteGlobalOutputs(configs[j], decons[j], (float)elapsed[j]);
				FreeInputs(inputs[j]);
				FreeDecon(decons[j]);
			}
		}
		free(configs);
		free(inputs);
		free(decons);
		free(elapsed);
	}
	omp_set_max_active_levels(previous_max_active_levels);
	omp_set_dynamic(previous_dynamic);
	return 1;
#endif
}

static int full_output_rawflag(const int rawflag) {
	if (rawflag == 2) { return 0; }
	if (rawflag == 3) { return 1; }
	return rawflag;
}

static int has_full_decon_outputs(const hid_t file_id, const int num) {
	char path[1024];
	const char* datasets[] = { "mz_grid", "mass_grid", "charge_data" };

	for (int i = 0; i < num; i++) {
		for (int j = 0; j < 3; j++) {
			sprintf(path, "/ms_dataset/%d/%s", i, datasets[j]);
			if (!H5LTpath_valid(file_id, path, 1)) { return 0; }
		}
	}
	return 1;
}

static int ensure_full_decon_outputs(int argc, char* argv[], Config config, const int num,
	const char* reason) {
	if (has_full_decon_outputs(config.file_id, num)) { return 1; }

	config.rawflag = full_output_rawflag(config.rawflag);
	config.silent = 1;
	printf("Full per-spectrum outputs are required for %s; generating them now.\n", reason);
	for (int i = 0; i < num; i++) {
		printf("Spectrum %d/%d\n", i + 1, num);
		config.metamode = i;
		run_unidec(argc, argv, config);
	}

	if (!has_full_decon_outputs(config.file_id, num)) {
		fprintf(stderr, "Unable to generate the per-spectrum mz, mass, and charge grids required for %s.\n",
			reason);
		return 0;
	}
	return 1;
}


int run_metaunidec(int argc, char* argv[], Config config) {
	clock_t starttime;
	starttime = clock();
	//Get Length
	int num = 0;
	int result = 0;
	config.file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
	num = int_attr(config.file_id, "/ms_dataset", "num", num);
	if (num > 20) { config.silent = 1; }
	const int meta_output_silent = config.silent;

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

	// Charge extraction during -all needs the grids omitted by Fast Profile and
	// Fast Centroid. Generate them in the first pass instead of deconvolving twice.
	if (mode == 4 && (config.exchoice == 6 || config.exchoice == 7) && config.rawflag > 1) {
		printf("Full per-spectrum outputs are required for charge extraction; generating them now.\n");
		config.rawflag = full_output_rawflag(config.rawflag);
	}

	if (mode == 1 || mode == 2 || mode == 4 || mode==0)
	{
		int used_parallel_decon = 0;
		if ((mode == 0 || mode == 4) && can_parallel_fast_decon(config, num)) {
			// Processing and all HDF5 input traffic remain serial.
			for (int i = 0; i < num; i++) {
				printf("Spectrum %d/%d\n", i + 1, num);
				config.silent = 1;
				config.metamode = i;
				process_data(argc, argv, config);
			}
			used_parallel_decon = run_parallel_fast_decon(config, num, 0);
		}
		else if (mode == 1) {
			used_parallel_decon = run_parallel_fast_decon(config, num, 1);
		}

		if (!used_parallel_decon) {
			//Iterate through files
			for (int i = 0; i < num; i++)
			{
				printf("Spectrum %d/%d\n", i + 1, num);
				// Per-spectrum peak detection/scoring is not needed for MetaUniDec
				// grids or merged peaks. Dedicated -peaks and -scanpeaks modes can
				// generate per-spectrum peak products later when requested.
				config.silent = 1;
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
		config.silent = meta_output_silent;
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
			if (!ensure_full_decon_outputs(argc, argv, config, num, "charge extraction")) {
				result = 12;
				goto cleanup;
			}
			printf("Charge Extraction Avg\n");
			config.exchoicez = 1;
			config.exnormz = config.exnorm;
			charge_peak_extracts(argc, argv, config, 0);
		}
		else if (config.exchoice == 7)
		{
			if (!ensure_full_decon_outputs(argc, argv, config, num, "charge extraction")) {
				result = 12;
				goto cleanup;
			}
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
		if (!ensure_full_decon_outputs(argc, argv, config, num, "ultra extraction")) {
			result = 12;
			goto cleanup;
		}
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
		if (!ensure_full_decon_outputs(argc, argv, config, num, "charge extraction")) {
			result = 12;
			goto cleanup;
		}
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
		if (!ensure_full_decon_outputs(argc, argv, config, num, "scan peak scoring")) {
			result = 12;
			goto cleanup;
		}
		printf("Getting Scan Scores\n");
		get_scan_scores(argc, argv, config);
	}

	cleanup:
	// Close the file
	H5Fclose(config.file_id);

	clock_t end = clock();
	float totaltime = (float)(end - starttime) / CLOCKS_PER_SEC;
	printf("Done in %f s\n", totaltime);
	return result;
}

int run_chromatogram(int argc, char* argv[], Config config) {
	clock_t starttime;
	starttime = clock();
	//Get Length
	int num = 0;
	config.file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
	num = int_attr(config.file_id, "/ms_dataset", "num", num);

	// Create an array of spectra
	Spectrum * spectra = (Spectrum*)malloc(sizeof(Spectrum) * num);
	if (spectra == NULL) {
		printf("Error allocating memory for spectra array\n");
		exit(1);
	}

	//Iterate through files
	for (int i = 0; i < num; i++) {
		// Run either deconvolution or processing
		config.metamode = i;
		process_data(argc, argv, config);
		run_unidec(argc, argv, config);
	}

	// Close the file
	H5Fclose(config.file_id);

	// Free memory
	free(spectra);

	clock_t end = clock();
	float totaltime = (float)(end - starttime) / CLOCKS_PER_SEC;
	printf("Done in %f s\n", totaltime);
	return 0;
}
