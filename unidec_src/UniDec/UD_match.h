/*
* UD_match.h
*
*  Created on : 9/28/22
* Author : Michael.Marty
*/

//
// 
// Copyright 2022 University of Arizona
//
//

void lindex_to_nindex(const int olen, const unsigned long long lindex, const unsigned long long * mults, int* nindex)
{
	unsigned long long running_index = lindex;
	for (int i = 0; i < olen; i++)
	{
		unsigned long long ni = running_index/mults[i];
		nindex[i] = (int) ni;
		running_index -= ni * mults[i];
	}
}

int nindex_to_lindex(const int olen, const int* nindex, const int* dims)
{
	int output = 0;
	for (int i = 0; i < olen; i++)
	{
		int mult = 1;
		for (int j = olen-1; j > i; j--)
		{
			mult *= dims[j];
		}
		output += nindex[i] * mult;
	}
	return output;
}

void read_ofile(const char* infile, const int olength, float* basemasses, float* monomasses, int* minn, int* maxn, char* onames)
{
	FILE* file_ptr;
	int i;
	char a[500] = { 0 };
	char b[500] = { 0 };
	char c[500] = { 0 };
	char d[500] = { 0 };
	char e[1000] = { 0 };

	file_ptr = fopen(infile, "r");
	if (file_ptr == 0) { printf("Error Opening %s\n", infile); exit(1); }

	for (i = 0; i < olength; i++)
	{
		int match = fscanf(file_ptr, "%s %s %s %s %999[^\n]%*c\n", a, b, c, d, e);
		if (match != 0) {
			printf("%s %s %s %s %s\n", a, b, c, d, e);
			basemasses[i] = atof(a);
			monomasses[i] = atof(b);
			minn[i] = atoi(c);
			maxn[i] = atoi(d);
			strcpy(&onames[i], e);
		}
	}

	fclose(file_ptr);
}


int unidec_match(int argc, char* argv[]) 
{
	printf("Running Match\n");

	//char buffer[1024];
	//setvbuf(stdout, buffer, _IOLBF, sizeof(buffer));

	clock_t starttime;
	starttime = clock();

	char* ofile = argv[2];
	char* pfile = argv[3];

	int olen = getfilelength(ofile);
	int plen = getfilelength(pfile);

	printf("Ofile Length: %d\n", olen);
	printf("Pfile Length: %d\n", plen);

	float* peakmasses, *peakints;
	peakmasses = calloc(plen, sizeof(float));
	peakints = calloc(plen, sizeof(float));

	readfile(pfile, plen, peakmasses, peakints);
	free(peakints);

	//Sort Peak Masses
	/*
	float* psorted = calloc(plen, sizeof(float));
	memcpy(psorted, peakmasses, sizeof(float) * plen);
	qsort(psorted, plen, sizeof(float), compare_highest_first);*/

	for (int i = 0; i < plen; i++)
	{
		printf("Peak: %f\n", peakmasses[i]);
	}

	float* basemasses, * monomasses;
	int* mins, * maxs;
	char * onames;
	basemasses = calloc(olen, sizeof(float));
	monomasses = calloc(olen, sizeof(float));
	mins = calloc(olen, sizeof(int));
	maxs = calloc(olen, sizeof(int));
	onames = calloc(olen, olen*1000);

	read_ofile(ofile, olen, basemasses, monomasses, mins, maxs, onames);

	//
	// Maybe this just would slow things down...
	// Put in the sorted pfile
	// Need to sort the ofile
	// Need to make it skip ahead when it's already too big
	// Need to take a break first...
	//

	int* dims;
	int* nindex;

	dims = calloc(olen, sizeof(int));
	nindex = calloc(olen, sizeof(int));

	for (int i = 0; i < olen; i++)
	{
		dims[i] = maxs[i] - mins[i] + 1;
	}

	unsigned long long maxn = 1;
	for (int i = 0; i < olen; i++)
	{
		maxn *= (unsigned long long) dims[i];
		//printf("Size: %llu\n", maxn);
	}
	printf("Max Number: %llu\n", maxn);

	unsigned long long *mults;
	mults = calloc(olen, sizeof(unsigned long long));
	
	for (int i = 0; i < olen; i++)
	{
		unsigned long long mult = 1;
		for (int j = olen - 1; j > i; j--)
		{
			mult *= dims[j];
		}
		mults[i] = mult;
	}

	float tolerance = 1000;

	int matchindex = 0;
	int maxmatches = plen * 1000;

	unsigned long long* matchindexes;
	matchindexes = calloc(maxmatches, sizeof(unsigned long long));

	float* errors;
	errors = calloc(maxmatches, sizeof(float));

	int* matches;
	matches = calloc(maxmatches, sizeof(int));
	

	//#pragma omp parallel for
	for (unsigned long long i = 0; i < maxn; i++)
	{		
		if (i%20000000==0 && i!=0){
			printf("Iteration: %llu\t",i);
			float percent = (float) i / (float) maxn;
			printf("Percent: %f\t", percent*100);
			clock_t inttime = clock();
			float intsec = (float)(inttime - starttime) / CLOCKS_PER_SEC;
			if (percent != 0) {
				printf("Time Remaining: %f (hrs)\t", intsec / percent / 3600);
				printf("Time per iteration: %f (us)\n", 1e6 * intsec / (float)i);
			}
		}

		lindex_to_nindex(olen, i, mults, nindex);
		
		float value = 0;
		for (int j = 0; j < olen; j++)
		{
			value += basemasses[j] + (float) nindex[j] * monomasses[j];
		}
		
		for (int j=0; j<plen; j++)
		{
			float abserror = fabsf(value - peakmasses[j]);
			//if (error < tolerance)
			//{continue;}
			//if (error > tolerance){break;}
			if (abserror < tolerance) {
				//printf("%f %d\n", value, j);
				matches[matchindex] = j;
				matchindexes[matchindex] = i;
				errors[matchindex] = abserror;
				matchindex++;
				if (matchindex == maxmatches)
				{
					matchindex = maxmatches-1;
				}
			}
		}
	}
	
	clock_t end = clock();
	float totaltime = (float)(end - starttime) / CLOCKS_PER_SEC;
	printf("Done in %f seconds!\n", totaltime);
	
	char str1[10000];

	//sprintf(str1, "Matchindex: %d\n", matchindex);
	for(int i=0; i < matchindex; i++)
	{	
		sprintf(str1, "%d %llu \n", matches[i], matchindexes[i]);
	}
	puts(str1);

	free(peakmasses);
	free(monomasses);
	free(basemasses);
	free(onames);
	free(mins);
	free(maxs);
	free(dims);
	free(nindex);
	free(mults);
	free(matches);
	free(matchindexes);
	free(errors);

	return 0;
}