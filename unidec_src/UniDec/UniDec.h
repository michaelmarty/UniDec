/*
 * UniDec.h
 *
 *  Created on: 1 Jul 2014
 *      Author: Michael.Marty
 */

/*
#ifndef UNIDEC_H_
#define UNIDEC_H_
#endif*/ /* UNIDEC_H_ */

#ifndef UNIDEC_HEADER
#define UNIDEC_HEADER

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "hdf5.h"
#include "hdf5_hl.h"
//#include "UD_sg.h"

void mh5readfile3d(hid_t file_id, char *dataname, int lengthmz, float *dataMZ, float *dataInt, float *data3);
void mh5readfile2d(hid_t file_id, char* dataname, int lengthmz, float* dataMZ, float* dataInt);
void mh5readfile1d(hid_t file_id, char* dataname, float* data);
int mh5getfilelength(hid_t file_id, char *dataname);
void strjoin(const char* s1, const char* s2, char* newstring);
void mh5writefile1d(hid_t file_id, char* dataname, int length, float* data1);
void mh5writefile2d(hid_t file_id, char* dataname, int length, float* data1, float* data2);
void mh5writefile2d_grid(hid_t file_id, char* dataname, int length1, int length2, float* data1);
void delete_group(hid_t file_id, char* dataname);

typedef struct Input Input;

struct Input {
	float* dataMZ;
	float* dataInt;
	float* testmasses;
	int* nztab;
	float* mtab;
	char* barr;
	int* isotopepos;
	float* isotopeval;
	float isoparams[10];
};

Input SetupInputs()
{
	Input inp;
	inp.dataMZ = NULL;
	inp.dataInt = NULL;
	inp.testmasses = NULL;
	inp.nztab = NULL;
	inp.mtab = NULL;
	inp.barr = NULL;
	inp.isotopepos = NULL;
	inp.isotopeval = NULL;
	//Set default isotope parameters. These can be overwritten by config file.
	float temp[10] = { 1.00840852e+00, 1.25318718e-03, 2.37226341e+00, 8.19178000e-04, -4.37741951e-01, 6.64992972e-04, 9.94230511e-01, 4.64975237e-01, 1.00529041e-02, 5.81240792e-01 };
	memcpy(inp.isoparams, temp, sizeof(temp));
	
	return inp;
}

void FreeInputs(Input inp)
{
	free(inp.dataMZ);
	free(inp.dataInt);
	free(inp.nztab);
	free(inp.mtab);
	free(inp.testmasses);
	free(inp.barr);
	free(inp.isotopepos);
	free(inp.isotopeval);
}

typedef struct Decon Decon;

struct Decon {
	float* fitdat;
	float* baseline;
	float* noise;
	float* massgrid;
	float* massaxis;
	float* massaxisval;
	float* blur;
	float* newblur;
	float* peakx;
	float* peaky;
	float* dscores;
	float error;
	float rsquared;
	int iterations;
	float uniscore;
	float conv;
	float threshold;
	int mlen;
	int plen;
};

Decon SetupDecon() {
	Decon decon;
	decon.fitdat = NULL;
	decon.baseline = NULL;
	decon.noise = NULL;
	decon.massgrid = NULL;
	decon.massaxis = NULL;
	decon.massaxisval = NULL;
	decon.blur = NULL;
	decon.newblur = NULL;
	decon.peakx = NULL;
	decon.peaky = NULL;
	decon.dscores = NULL;
	decon.error = 0;
	decon.rsquared = 0;
	decon.iterations = 0;
	decon.uniscore = 0;
	decon.conv = 0;
	decon.threshold = 0;
	decon.mlen = 0;
	decon.plen = 0;
	return decon;
}

void FreeDecon(Decon decon)
{
	free(decon.fitdat);
	free(decon.baseline);
	free(decon.noise);
	free(decon.massgrid);
	free(decon.massaxis);
	free(decon.massaxisval);
	free(decon.blur);
	free(decon.newblur);
	free(decon.peakx);
	free(decon.peaky);
	free(decon.dscores);
}

typedef struct Config Config;

struct Config
{
	char infile[500];
	char outfile[500];
	int numit;
	int numz;
	int endz;
	int startz;
	float zsig;
	float psig;
	float beta;
	float mzsig;
	float msig;
	float molig;
	float massub;
	float masslb;
	int psfun;
	float mtabsig;
	char mfile[500];
	char manualfile[500];
	int mflag;
	float massbins;
	int limitflag;
	float psthresh;
	int speedyflag;
	int linflag;
	int aggressiveflag;
	float adductmass;
	int rawflag;
	float nativezub;
	float nativezlb;
	int poolflag;
	int manualflag;
	float intthresh;
	float peakshapeinflate;
	float killmass;
	int fixedmassaxis;
	int isotopemode;
	int filetype;
	int imflag;
	//IM Parameters
	float dtsig;
	float csig;
	float ccsub;
	float ccslb;
	float ccsbins;
	float temp;
	float press;
	float volt;
	float tcal1;
	float tcal2;
	float twaveflag;
	float hmass;
	float to;
	float len;
	float edc;
	float nativeccsub;
	float nativeccslb;
	int baselineflag;
	int noiseflag;
	int zout;
	int metamode;
	float minmz;
	float maxmz;
	int mzbins;
	float bsub;
	float datareduction;
	float peakwin;
	float peakthresh;
	float exwindow;
	int exchoice;
	int exchoicez;
	float exthresh;
	int exnorm;
	int exnormz;
	int peaknorm;
	int orbimode;
	int datanorm;
	//Experimental Parameters
	int filterwidth;
	float zerolog;
	int lengthmz;
	int mfilelen;
	int isolength;
	hid_t file_id;
	char dataset[1024];
};

Config SetDefaultConfig()
{
	Config config;
	strcpy(config.infile, "default_data.txt");
	strcpy(config.outfile, "default_output");
	strcpy(config.mfile, "default_mfile.txt");
	config.numit = 50;
	config.endz = 100;
	config.startz = 1;
	config.zsig = 1;
	config.psig = 1;
	config.beta = 0;
	config.mzsig = 15;
	config.msig = 0;
	config.molig = 0;
	config.massub = 5000000;
	config.masslb = 100;
	config.psfun = 0;
	config.mtabsig = 0;
	config.mflag = 0;
	config.massbins = 100;
	config.limitflag = 0;
	config.psthresh = 6;
	config.speedyflag = 0;
	config.aggressiveflag = 0;
	config.adductmass = 1.007276467;
	config.rawflag = 1;
	config.nativezub = 100;
	config.nativezlb = -200;
	config.poolflag = 1;
	config.manualflag = 0;
	config.intthresh = 0;
	config.peakshapeinflate = 1;
	config.killmass = 0;
	config.fixedmassaxis = 0;
	config.isotopemode = 0;
	config.filetype = 0;
	config.imflag = 0;
	config.linflag = -1;
	//IM Parameters
	config.dtsig=0.2;
	config.csig=1;
	config.ccsub=20000;
	config.ccslb=-20000;
	config.ccsbins=100;
	config.temp=0;
	config.press=2;
	config.volt=50;
	config.tcal1=0.3293;
	config.tcal2 = 6.3597;
	config.twaveflag = -1;
	config.hmass = 4.002602;
	config.to = 0.97;
	config.len = 0.18202;
	config.edc = 1.57;
	config.nativeccsub = 20000;
	config.nativeccslb = -20000;
	config.zout = 0;
	config.baselineflag = 1;
	config.noiseflag = 0;
	config.metamode = -2;
	config.minmz = -1;
	config.maxmz = -1;
	config.mzbins = 0;
	config.bsub = 0;
	config.datareduction = 0;
	config.peakwin = 500;
	config.peakthresh = 0.1;
	config.exchoice = 0;
	config.exchoicez = 1;
	config.exthresh = 10;
	config.exnorm = 0;
	config.exnormz = 0;
	config.peaknorm = 1;
	config.exwindow = 0;
	config.orbimode = 0;
	config.datanorm = 1;
	//Experimental
	config.filterwidth = 20;
	config.zerolog = -12;
	config.lengthmz = 0;
	config.mfilelen = 0;
	config.isolength = 0;
	config.file_id = 0;
	return config;
}

Config PostImport(Config config)
{
	//Convert gaussian FWHM to sigma
	if (config.psfun == 0) { config.mzsig = config.mzsig / 2.35482; }
	config.dtsig = config.dtsig / 2.35482;
	//Check whether to turn off or on the config.limitflag. Limit flag on means it will select only the mass values closest to the mass values from mfile.
	if (config.mflag == 1 && config.mtabsig == 0) { config.limitflag = 1; }

	config.numz = config.endz - config.startz + 1;

	//If linflag is active, overwrite speedy}
	if (config.linflag !=-1){
		if (config.linflag != 2) { config.speedyflag = 1; }
		else { config.speedyflag = 0; }
	}
	if (config.speedyflag == 1) { printf("Speedy mode: Assuming Linearized Data\n"); }

	if(config.filetype==0)
	{
		//Print inputs for check
		printf("infile = %s\n", config.infile);
		//printf("outfile = %s\n", config.outfile);
		//printf("\n");
	}

	//Check to see if the mass axis should be fixed
	if (config.massub<0 || config.masslb<0) {
		config.fixedmassaxis = 1;
		config.massub = fabs(config.massub);
		config.masslb = fabs(config.masslb);
	}

	if (config.twaveflag == -1 && config.imflag == 1) { printf("\n\nNeed to define twaveflag for CCS calculation\n\n"); }

	//bug proofing so we don't get a 1/0 problem
	if (config.msig == 0) { config.msig = 0.00001; }
	if (config.zsig == 0) { config.zsig = 0.00001; }
	if (config.massbins == 0) { config.massbins = 1; }

	//Convert aggressiveflag to baselineflag
	if (config.aggressiveflag == 1 || config.aggressiveflag == 2) { config.baselineflag = 1; }
	else { config.baselineflag = 0; }

	//Experimental correction. Not sure why this is necessary.
	if (config.psig < 0) { config.mzsig /= 3; }

	return config;
}

Config LoadConfig(Config config,const char *filename)
{
	// We assume argv[1] is a filename to open
	FILE *file = fopen(filename, "r");

	if (file == 0)
	{
		printf("\nCould not open configuration file: \n \nSomething is wrong...\n\n");
	}
	//Read parameters from configuration file
	//Configuration file should be formatted: name value
	else
	{
		char x[500];
		char y[500];
		//printf("\nRead from file:");
		while (fscanf(file, "%s %500[^\n]", x, y) != EOF)
		{
			//printf( "read in: %s %s \n", x,y );
			if (strstr(x, "input") != NULL) { strcpy(config.infile, y); }// printf(" input");
			if (strstr(x, "output") != NULL){ strcpy(config.outfile, y);}//  printf(" output"); }
			if (strstr(x, "mfile") != NULL){ strcpy(config.mfile, y); config.mflag = 1;  }// printf(" mfile"); }
			if (strstr(x, "numit") != NULL){ config.numit = atoi(y); }// printf(" numit"); }
			//if (strstr(x, "numz") != NULL){ config.numz = atoi(y); printf(" numz"); }
			if (strstr(x, "startz") != NULL){ config.startz = atoi(y); }// printf(" startz"); }
			if (strstr(x, "endz") != NULL) { config.endz = atoi(y); }// printf(" endz"); }
			if (strstr(x, "zzsig") != NULL){ config.zsig = atof(y); }// printf(" zzsig"); }
			if (strstr(x, "psig") != NULL) { config.psig = atof(y); }// printf(" psig"); }
			if (strstr(x, "beta") != NULL) { config.beta = atof(y); }// printf(" beta"); }
			if (strstr(x, "mzsig") != NULL){ config.mzsig = atof(y); }// printf(" mzsig"); }
			if (strstr(x, "msig") != NULL){ config.msig = atof(y); }// printf(" msig"); }
			if (strstr(x, "molig") != NULL){ config.molig = atof(y); }// printf(" molig"); }
			if (strstr(x, "massub") != NULL){ config.massub = atof(y); }// printf(" massub"); }
			if (strstr(x, "masslb") != NULL){ config.masslb = atof(y); }// printf(" masslb"); }
			if (strstr(x, "psfun") != NULL){ config.psfun = atoi(y); }// printf(" psfun"); }
			if (strstr(x, "mtabsig") != NULL){ config.mtabsig = atof(y); }// printf(" mtabsig"); }
			if (strstr(x, "massbins") != NULL){ config.massbins = atof(y); }// printf(" massbins"); }
			if (strstr(x, "psthresh") != NULL){ config.psthresh = atof(y); }// printf(" psthresh"); }
			if (strstr(x, "speedy") != NULL){ config.speedyflag = atoi(y);}//  printf(" speedy"); }
			if (strstr(x, "aggressive") != NULL){ config.aggressiveflag = atoi(y);}//  printf(" aggressive"); }
			if (strstr(x, "adductmass") != NULL){ config.adductmass = atof(y);}//  printf(" adductmass"); }
			if (strstr(x, "rawflag") != NULL){ config.rawflag = atoi(y); }// printf(" rawflag"); }
			if (strstr(x, "nativezub") != NULL){ config.nativezub = atof(y); }// printf(" nativezub"); }
			if (strstr(x, "nativezlb") != NULL){ config.nativezlb = atof(y); }// printf(" nativezlb"); }
			if (strstr(x, "poolflag") != NULL){ config.poolflag = atoi(y); }// printf(" poolflag"); }
			if (strstr(x, "manualfile") != NULL){ config.manualflag = 1; strcpy(config.manualfile, y);}//  printf(" manualfile"); }
			if (strstr(x, "intthresh") != NULL){ config.intthresh = atof(y); }// printf(" intthresh"); }
			if (strstr(x, "peakshapeinflate") != NULL){ config.peakshapeinflate = atof(y); }// printf(" peakshapeinflate"); }
			if (strstr(x, "killmass") != NULL){ config.killmass = atof(y); }// printf(" killmass"); }
			if (strstr(x, "isotopemode") != NULL){ config.isotopemode = atoi(y); }// printf(" isotopemode"); }
			if (strstr(x, "orbimode") != NULL) { config.orbimode = atoi(y); }// printf(" orbimode"); }
			if (strstr(x, "imflag") != NULL) { config.imflag = atoi(y); }// printf(" imflag"); }
			if (strstr(x, "linflag") != NULL) { config.linflag = atoi(y); }// printf(" linflag"); }
			//IM Parameters
			if (strstr(x, "csig") != NULL) { config.csig = atof(y); }// printf(" csig"); }
			if (strstr(x, "dtsig") != NULL) { config.dtsig = atof(y); }// printf(" dtsig"); }
			if (strstr(x, "ccsub") != NULL) { config.ccsub = atof(y); }// printf(" ccsub"); }
			if (strstr(x, "ccslb") != NULL) { config.ccslb = atof(y); }// printf(" ccslb"); }
			if (strstr(x, "ccsbins") != NULL) { config.ccsbins = atof(y); }// printf(" ccsbins"); }
			if (strstr(x, "temp") != NULL) { config.temp = atof(y); }// printf(" temp"); }
			if (strstr(x, "pressure") != NULL) { config.press = atof(y); }// printf(" pressure"); }
			if (strstr(x, "volt") != NULL) { config.volt = atof(y); }// printf(" volt"); }
			if (strstr(x, "gasmass") != NULL) { config.hmass = atof(y); }// printf(" gasmass"); }
			if (strstr(x, "tnaught") != NULL) { config.to = atof(y); }// printf(" to"); }
			if (strstr(x, "tcal1") != NULL) { config.tcal1 = atof(y); }// printf(" tcal1"); }
			if (strstr(x, "tcal2") != NULL) { config.tcal2 = atof(y); }// printf(" tcal2"); }
			if (strstr(x, "edc") != NULL) { config.edc = atof(y); }// printf(" edc"); }
			if (strstr(x, "zout") != NULL) { config.zout = atoi(y); }// printf(" zout"); }
			if (strstr(x, "twaveflag") != NULL) { config.twaveflag = atoi(y); }// printf(" twaveflag"); }
			if (strstr(x, "ubnativeccs") != NULL) { config.nativeccsub = atof(y); }// printf(" ubnativeccs"); }
			if (strstr(x, "lbnativeccs") != NULL) { config.nativeccslb = atof(y); }// printf(" lbnativeccs"); }
			if (strstr(x, "driftlength") != NULL) { config.len = atof(y); }// printf(" driftlength"); }
			if (strstr(x, "baselineflag") != NULL) { config.baselineflag = atoi(y); }// printf(" baselineflag"); }
			if (strstr(x, "noiseflag") != NULL) { config.noiseflag = atoi(y); }// printf(" noiseflag"); }
			//Experimental
			if (strstr(x, "filterwidth") != NULL) { config.filterwidth = atoi(y); }// printf(" filterwidth"); }
			if (strstr(x, "zerolog") != NULL) { config.zerolog = atof(y); }// printf(" zerolog"); }
			//Peak Parameters
			if (strstr(x, "peakwindow") != NULL) { config.peakwin = atof(y); }// printf(" peakwindow"); }
			if (strstr(x, "peakthresh") != NULL) { config.peakthresh = atof(y); }// printf(" peakthresh"); }
			if (strstr(x, "peaknorm") != NULL) { config.peaknorm = atoi(y); }// printf(" peaknorm"); }
		}
		//printf("\n\n");
	}
	fclose(file);

	config=PostImport(config);

	return config;

}


//Gives command line help options
void PrintHelp()
{
	printf("\nUniDec: Universal Deconvolution of Mass and Ion Mobility Spectra\n\nWritten by Michael Marty\n\twith contributions from Erik Marklund and Andrew Baldwin\n");
	printf("Copyright University of Oxford 2016 and University of Arizona 2017\n");
	printf("\nUniDec runs off a configuration text or hdf5 file.\nExample usage: UniDec.exe conf.txt\n");
	printf("Text configuration file should be written for each line: keyword argument\n");
	printf("HDF5 files should consult the Python API source code for a model of directory and metadata construction.\n");
	printf("File names should include the path if not in the current directory\n");
	printf("\nPossible configuration keywords include:\n");
	printf("\t\"imflag\" \tSelects between MS and IM-MS\n");
	printf("\t\t\t\t0=MS\t1=IM-MS\n");
	printf("\t\"input\" \tInput text file name in x y form (MS) or x y z (IM-MS)\n");
	printf("\t\"output\" \tHeader for output file names\n");
	printf("\t\"numit\" \tNumber of iterations\n");
	printf("\t\"startz\" \tMinimum charge value\n");
	printf("\t\"endz\"   \tMaximum charge value\n");
	printf("\t\"mzsig\" \tPeak Full-Width at Half Max\n");
	printf("\t\"psfun\" \tPeak shape function\n\t\t\t\t0=Gaussian\n\t\t\t\t1=Lorentzian\n\t\t\t\t2=Split Gaussian\\Lorentzian\n");
	printf("\t\"zzsig\" \tCharge state smoothing width parameter\n");
	printf("\t\"psig\" \tPoint smoothing width parameter\n");
	printf("\t\"beta\" \tBoltzman/Softmax factor for charge state distributions.\n");
	printf("\t\"molig\" \tMass to be included in smooth\n");
	printf("\t\"msig\"   \tWidth of mass smoothing\n");
	printf("\t\"mfile\" \tText file list of limiting masses\n");
	printf("\t\"mtabsig\" \tWindow on mass limitations from mfile in Da\n");
	printf("\t\"massub\" \tUpper bound on mass\n");
	printf("\t\"masslb\" \tLower bound on mass\n");
	printf("\t\"massbins\" \tSize of output mass bin in Da\n");
	printf("\t\"adductmass\" \tMass of electrospray adduct in Da, usually 1\n");
	printf("\t\"rawflag\" \tTurns on convolution before transformation\n");
	printf("\t\t\t\t0=Deconvolved outputs\n");
	printf("\t\t\t\t1=Outputs reconvolved with peak shape\n");
	printf("\t\"nativezub\" \tLimits to native charge + this upper bound.\n");
	printf("\t\"nativezlb\" \tLimits to native charge - this lower bound.\n");
	printf("\t\"manualfile\" \tManual assignments file: M/Zvalue window charge\n");
	printf("\t\"intthresh\" \tLimits to m/z values with intensities > threshold \n");
	printf("\t\"isotopemode\" \t0=off 1=monoisotopic 2=average: Uses isotope distributions in deconvolution (MS Only)\n");
	printf("\t\"orbimode\" \t0=off 1=on: Uses charge scaling by 1/z to scale intensities based on charge (MS Only)\n");
	printf("\t\"poolflag\" \tSpecifies how to transform from m/z to mass axis\n");
	printf("\t\t\t\t0=Integration\n");
	printf("\t\t\t\t1=Interpolation\n");
	printf("\t\"speedy\" \tAllows faster convolution on linearized data (MS Only)\n");
	printf("\t\t\t\t0=No speedup. Works for nonlinear data\n");
	printf("\t\t\t\t1=Speedup! Works only for linearized data\n");
	printf("\tIM-MS mode assumes linearized data and includes the following parameters:\n");
	printf("\t\"csig\"   \tWidth of CCS smoothing\n");
	printf("\t\"dtsig\"   \tFWHM or peak in drift time\n");
	printf("\t\"ccsub\"   \tUpper bound on CCS\n");
	printf("\t\"ccslb\"   \tLower bound on CCS\n");
	printf("\t\"ccsbins\"   \tDensity of CCS sampling\n");
	printf("\t\"ubnativeccs\"   \tLimits to native CCS + this upper bound.\n");
	printf("\t\"lbnativeccs\"   \tLimits to native CCS - this lower bound.\n");
	printf("\t\"zout\"   \tSpecific charge state to extract data from 3D final grid (optional).\n");
	printf("\t\"twaveflag\"   \tType of equation for calculating CCS from drift time\n");
	printf("\t\t\t\t0=Linear Cell\n");
	printf("\t\t\t\t\t\"temp\"=Temperature in K\n");
	printf("\t\t\t\t\t\"pressure\"=Pressure in Torr\n");
	printf("\t\t\t\t\t\"volt\"=Drift voltage\n");
	printf("\t\t\t\t\t\"driftlength\"=Length of drift cell\n");
	printf("\t\t\t\t\t\"tnaught\"=Dead time between IM cell and detector\n");
	printf("\t\t\t\t1=T-Wave Log Calibration\n");
	printf("\t\t\t\t\t\"tcal1\"=Calibration paramter 1 (P1)\n");
	printf("\t\t\t\t\t\"tcal2\"=Calibration parmater 2 (P2)\n");
	printf("\t\t\t\t\t\"EDC\"=Instrumental constant\n");
	printf("\t\t\t\t\tReduced CCS = Exp(P1 * log(Reduced Drift Time) + P2)\n");
	printf("\t\t\t\t2=T-Wave Linear Calibration\n");
	printf("\t\t\t\t\t Reduced CCS = P1 * Reduced Drift Time + P2\n");
	printf("\t\t\t\t3=T-Wave Power Law Calibration\n");
	printf("\t\t\t\t\tReduced CCS =P1 * (Reduced Drift Time ^ P2)\n");
	printf("\nEnjoy! Please report bugs to Michael Marty (mtmarty@email.arizona.edu) v.1141.\n");
	//printf("\nsize of: %d",sizeof(char));
}

//Print a float array
void floatPrint(const float* array, const int length)
{
	for (int i = 0; i < length; i++)
	{
		printf("%f\n", array[i]);
	}
}

//Print an int array
void IntPrint(const int* array, const int length)
{
	for (int i = 0; i < length; i++)
	{
		printf("%d\n", array[i]);
	}
}

//Calculate Average
float Average(const int length, const float* xarray)
{
	float temp1 = 0;
	float temp2 = (float) length;
	for (int i = 0; i < length; i++)
	{
		temp1 += xarray[i];
	}
	if (temp2 == 0) { return 0; }
	return temp1 / temp2;
}

float ndis(float x, float y, float sig) 
{
	if (sig == 0) { return 0; }
	return 1 / (sig*2.50663)*exp(-(pow(x - y, 2.)) / (2. * sig*sig));
}

//Actual Modulus operator rather than Remainder operator %
int mod(int a, int b) { int r = a % b; return r < 0 ? r + b : r; }

//Reads in x y file.
void readfile(char *infile,int lengthmz,float *dataMZ,float *dataInt)
{
  FILE *file_ptr;
  int i;
  char x[500];
  char y[500];

  file_ptr=fopen(infile,"r");

  if ( file_ptr == 0 )
    {
    printf( "Could not open %s file\n",infile);
    exit(10);
    }
  else{

  for(i=0;i<lengthmz;i++)
    {
      fscanf(file_ptr,"%s %s",x,y);
      dataMZ[i]=atof(x);
      dataInt[i]=atof(y);
    }
  }
  fclose(file_ptr);

}

//Reads in x y z file.
void readfile3(char *infile,int lengthmz,float *array1,float *array2,float *array3)
{
  FILE *file_ptr;
  int i;
  char x[500];
  char y[500];
  char z[500];

  file_ptr=fopen(infile,"r");

  if ( file_ptr == 0 )
    {
    printf( "Could not open %s file\n",infile);
    exit(10);
    }
  else{

  for(i=0;i<lengthmz;i++)
    {
      fscanf(file_ptr,"%s %s %s",x,y,z);
      array1[i]=atof(x);
      array2[i]=atof(y);
      array3[i]=atof(z);
    }
  }
  fclose(file_ptr);

}

//Reads in single list of values
void readmfile(char *infile,int mfilelen,float *testmasses)
{
  FILE *file_ptr;
  int i;
  char input[500];

  file_ptr=fopen(infile,"r");

  if ( file_ptr == 0 )
    {
    printf( "Could not open %s file\n",infile);
    exit(20);
    }
  else{

  for(i=0;i<mfilelen;i++)
    {
	fscanf(file_ptr,"%s",input);
	testmasses[i]=atof(input);
    }
  }
  fclose(file_ptr);
}

// count the number of lines we have in datafile
int getfilelength(char *infile)
{
	FILE *file_ptr;
	int l = 0;
	char input[501];
	file_ptr = fopen(infile, "r");

	if (file_ptr == 0)
	{
		printf("Could not open %s file\n", infile);
		exit(9);
	}
	else
	{
		while (fscanf(file_ptr, "%500[^\n]%*c", input) != EOF)
		{
			l += 1;
		}
	}
	fclose(file_ptr);
	return l;
}

//Slow nearest unsorted.
int nearunsorted(float *testmasses, float point, int lengthtest)
{
	float minval = fabs(point - testmasses[0]);
	float val = testmasses[0];
	float difftest;
	int pos=0;
	for (int i = 1; i<lengthtest; i++)
	{
		difftest = fabs(point - testmasses[i]);
		if (difftest<minval)
		{
			minval = difftest;
			val = testmasses[i];
			pos = i;
		}
	}
	return pos;
}

//Slow way to test if two points are within a cutoff distance of each other. Works on unsorted list.
int neartest(float *testmasses,float point,int lengthtest,float cutoff)
{
    float minval=fabs(point-testmasses[0]);
    float val=testmasses[0];
	for (int i = 0; i<lengthtest; i++)
    {
        float difftest=fabs(point-testmasses[i]);
		if (difftest<minval)
        {
			minval = difftest;
            val=testmasses[i];
        }
    }
    int test=0;
    if(fabs(point-val)<cutoff)
    {
        test=1;
    }
    return test;
}

//Fast way of finding the nearest data point in an ordered list.
int nearfast(const float *dataMZ,const float point,const int numdat)
{
    int start=0;
    int length=numdat-1;
    int end=0;
    int diff=length-start;
    while(diff>1)
        {
        if(point<dataMZ[start+(length-start)/2])
            {
            length=start+(length-start)/2;
            }
        else if(point==dataMZ[start+(length-start)/2])
            {
            end=start+(length-start)/2;
            length=start+(length-start)/2;
            start=start+(length-start)/2;
            return end;
            }
        else if(point>dataMZ[start+(length-start)/2])
            {
            start=start+(length-start)/2;
            }
        diff=length-start;
        }
    if(fabs(point-dataMZ[start])>=fabs(point-dataMZ[length]))
    {
        end=length;
    }
    else
    {
        end=start;
    }
    return end;
}

int nearfast_test(const float *dataMZ, const float point, const int numdat, float cutoff)
{
	int index = nearfast(dataMZ, point, numdat);
	float val = dataMZ[index];
	if (fabs(val - point) < cutoff)
	{
		return index;
	}
	else { return -1; }
}

//Perform a linear interpolation. I've left the code for other interpolation functions below but they don't seem to matter.
float LinearInterpolate(float y1,float y2,float mu)
{
   return(y1*(1-mu)+y2*mu);
}

float LinearInterpolatePosition(float x1, float x2, float x)
{
	if (x2 - x1 == 0) { return 0; }
	return (x - x1) / (x2 - x1);
}

float CubicInterpolate(float y0,float y1,float y2,float y3,float mu)
{float a0,a1,a2,a3,mu2;
   mu2 = mu*mu;
   a0 = y3 - y2 - y0 + y1;
   a1 = y0 - y1 - a0;
   a2 = y2 - y0;
   a3 = y1;
   return(a0*mu*mu2+a1*mu2+a2*mu+a3);}

float CRSplineInterpolate(float y0,float y1,float y2,float y3,float mu)
{float a0,a1,a2,a3,mu2;
   mu2 = mu*mu;
   a0 = -0.5*y0 + 1.5*y1 - 1.5*y2 + 0.5*y3;
   a1 = y0 - 2.5*y1 + 2*y2 - 0.5*y3;
   a2 = -0.5*y0 + 0.5*y2;
   a3 = y1;
   return(a0*mu*mu2+a1*mu2+a2*mu+a3);}

float clip(float x,float cutoff)
{
	if(x>cutoff){return x;}
	else{return 0;}
}

//Function for defining m/z peak shape. Other peak shape functions could be easily added her.
float mzpeakshape(float x,float y,float sig,int psfun)
{
	if (sig == 0) { printf("Error: mzpeakshape sigma is 0\n"); exit(103); }
    float result;
    if(psfun==0)
    {
        result=exp(-(pow(x-y,2))/(2*sig*sig)); 
    }
    else if(psfun==1)
    {
        result=pow((sig/2),2)/(pow((x-y),2)+pow((sig/2),2));
    }
    else if(psfun==2)
    {
        if(y<x)
        {
        result=exp(-(pow(x-y,2))/(2*sig*sig*0.180337));
        }
        else
        {
        result=(sig/2)*(sig/2)/(pow((x-y),2)+pow((sig/2),2));
        }
    }
    else
    {
        printf("Invalid Peak Function");
        exit(14);
    }
    return result;
}

inline int index2D(const int ncols, const int r, const int c)
{ return r*ncols + c; }

inline int indexmod(int length, int r, int c)
{
	return mod((c - r), length);
}

inline int index3D(const int ncols, const int nrows, const int r, const int c, const int d)
{
	return r*ncols*nrows + c*nrows+d;
}


float Max(const float *blur, const int length){
	float blurmax = 0;
	for (int i = 0; i<length; i++)
	{
		if (blur[i]>blurmax)
		{
			blurmax = blur[i];
		}
	}
	return blurmax;
}

float Sum(const float *blur, int length) {
	float sum = 0;
	for (int i = 0; i<length; i++)
	{
		sum += blur[i];
	}
	return sum;
}

float max1d(float *blur, int lengthmz) {
	float blurmax = blur[0];
	unsigned int i;
	for (i = 0; i<lengthmz; i++)
	{
		if (blur[i]>blurmax)
		{
			blurmax = blur[i];
		}
	}
	return blurmax;
}

float min1d(float *blur, int lengthmz) {
	float blurmin = blur[0];
	unsigned int i;
	for (i = 0; i<lengthmz; i++)
	{
		if (blur[i]<blurmin)
		{
			blurmin = blur[i];
		}
	}
	return blurmin;
}

int argmax(float *blur, int lengthmz)
{
	float max = blur[0];
	int pos = 0;
	unsigned int i;
	for (i = 0; i<lengthmz; i++)
	{
		if (blur[i]>max)
		{
			max = blur[i];
			pos = i;
		}
	}
	return pos;
}

//Convolution of neighborhood function with gaussian filter.
void blur_it(const int lengthmz,
                    const int numz,
                    const int numclose,
                    const int * __restrict closeind,
                    const float * __restrict closearray,
                    float * __restrict newblur,
                    const float * __restrict blur,
					const char * __restrict barr)
{
  if (numclose == 1)
  {
	  memcpy(newblur, blur, lengthmz*numz * sizeof(float));
  }
  else {
	#pragma omp parallel for schedule(auto)
	  for (int i = 0; i < lengthmz* numz; i++)
	  {
		float temp = 0;
		if (barr[i] == 1)
		{
			for (int k = 0; k < numclose; k++)
			{
				if (closeind[index2D(numclose, i, k)] != -1)
				{
					temp += closearray[index2D(numclose, i, k)] * blur[closeind[index2D(numclose, i, k)]];
				}
			}
		}
		newblur[i] = temp;
	  }
  }
}

//Charge state smooth using a mean filter of the log
void blur_it_mean(const int lengthmz,
                    const int numz,
                    const int numclose,
                    const int * __restrict closeind,
                    float * __restrict newblur,
                    const float * __restrict blur,
					const char * __restrict barr,
					const float * __restrict closearray,
					const float zerolog)
{
  if (numclose == 1)
  {
	  memcpy(newblur, blur, lengthmz*numz*sizeof(float));
  }
  else{
	#pragma omp parallel for schedule(auto)
	  for (int i = 0; i < lengthmz*numz; i++)
	  {
		float temp = 0;
		if (barr[i] == 1)
		{
			for (int k = 0; k < numclose; k++)
			{
				float temp2 = 0;
				if (closeind[index2D(numclose, i, k)] != -1)
				{
					temp2 = blur[closeind[index2D(numclose, i, k)]] *closearray[index2D(numclose, i, k)];
				}
				if (temp2 > 0) { temp += log(temp2); }
				else { temp += zerolog; }
			}
			temp = exp(temp / (float)numclose);
		}
		newblur[i] = temp;
	  }
  }
}

/*
//Charge state smooth using a mean filter of the log
void blur_it_geometric_mean(const int lengthmz,
	const int numz,
	const int numclose,
	const int* __restrict closeind,
	float* __restrict newblur,
	const float* __restrict blur,
	const char* __restrict barr,
	const float zerolog)
{	
	float mult = 10;
	if (numclose == 1)
	{
		memcpy(newblur, blur, lengthmz * numz * sizeof(float));
	}
	else {
#		pragma omp parallel for schedule(auto)
		for (int i = 0; i < lengthmz * numz; i++)
		{
			float temp = 0;
			if (barr[i] == 1)
			{
				for (int k = 0; k < numclose; k++)
				{
					float temp2 = 0;
					if (closeind[index2D(numclose, i, k)] != -1)
					{
						temp2 = pow(blur[closeind[index2D(numclose, i, k)]], mult);
					}
					temp += temp2;
				}
				temp = pow(temp, 1/mult)/ (float)numclose;
			}
			newblur[i] = temp;
		}
	}
}*/
/*
//Convolution of neighborhood function with gaussian filter.
void blur_it_hybrid1alt(const int lengthmz,
	const int numz,
	const int zlength,
	const int mlength,
	const int * __restrict closeind,
	const int * __restrict closemind,
	const int * __restrict closezind,
	const float * __restrict mdist,
	const float * __restrict zdist,
	float * __restrict newblur,
	const float * __restrict blur,
	const char * __restrict barr,
	const float* __restrict closearray,
	const float zerolog)
{
	int i, j, k, n;
	int numclose = zlength * mlength;
	if (numclose == 1)
	{
		memcpy(newblur, blur, lengthmz*numz * sizeof(float));
	}
	else {
	#pragma omp parallel for private (i,k,j, n), schedule(auto)
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				float temp = 0;
				if (barr[index2D(numz, i, j)] == 1)
				{
					for (k = 0; k < zlength; k++)
					{
						float temp2 = 0;
						for (n = 0; n < mlength; n++)
						{
							int m = index2D(mlength, k, n);
							if (closeind[index3D(numz, numclose, i, j, m)] != -1)
							{
								temp2 += blur[closeind[index3D(numz, numclose, i, j, m)]] * mdist[n] *closearray[index3D(numz, numclose, i, j, m)];
							}
						}
						if (temp2 > 0) { temp += log(temp2); }
						else { temp += zerolog; }
					}
					temp = exp(temp / (float)zlength);
				}
				newblur[index2D(numz, i, j)] = temp;
			}
		}
	}
}*/

//Convolution of neighborhood function with gaussian filter.
void blur_it_hybrid1(const int lengthmz,
	const int numz,
	const int zlength,
	const int mlength,
	const int * __restrict closeind,
	const int * __restrict closemind,
	const int * __restrict closezind,
	const float * __restrict mdist,
	const float * __restrict zdist,
	float * __restrict newblur,
	const float * __restrict blur,
	const char * __restrict barr,
	const float* __restrict closearray,
	const float zerolog)
{
	int i, j, k, n;
	int numclose = zlength * mlength;
	if (numclose == 1)
	{
		memcpy(newblur, blur, lengthmz*numz * sizeof(float));
	}
	else {
		#pragma omp parallel for private (i,k,j, n), schedule(auto)
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				float temp = 0;
				if (barr[index2D(numz, i, j)] == 1)
				{
					for (n = 0; n < mlength; n++)
					{
						float temp2 = 0;
						for (k = 0; k < zlength; k++)
						{
							int m = index2D(mlength, k, n);
							float temp3 = 0;
							if (closeind[index3D(numz, numclose, i, j, m)] != -1)
							{
								temp3 = blur[closeind[index3D(numz, numclose, i, j, m)]] *closearray[index3D(numz, numclose, i, j, m)];
							}
							if (temp3 > 0) { temp2 += log(temp3); }
							else { temp2 += zerolog; }
						}
						temp += exp(temp2 / (float)zlength) * mdist[n];
					}
				}
				newblur[index2D(numz, i, j)] = temp;
			}
		}
	}
}


//Convolution of neighborhood function with gaussian filter.
void blur_it_hybrid2(const int lengthmz,
	const int numz,
	const int zlength,
	const int mlength,
	const int * __restrict closeind,
	const int * __restrict closemind,
	const int * __restrict closezind,
	const float * __restrict mdist,
	const float * __restrict zdist,
	float * __restrict newblur,
	const float * __restrict blur,
	const char * __restrict barr,
	const float* __restrict closearray,
	const float zerolog)
{
	int i, j, k, n;
	int numclose = zlength * mlength;
	if (numclose == 1)
	{
		memcpy(newblur, blur, lengthmz*numz * sizeof(float));
	}
	else {
		#pragma omp parallel for private (i,k,j, n), schedule(auto)
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				float temp = 0;
				if (barr[index2D(numz, i, j)] == 1)
				{
					for (n = 0; n < mlength; n++)
					{
						float temp2 = 0;
						for (k = 0; k < zlength; k++)
						{
							int m = index2D(mlength, k, n);
							if (closeind[index3D(numz, numclose, i, j, m)] != -1)
							{
								temp2 += blur[closeind[index3D(numz, numclose, i, j, m)]]*zdist[k] * closearray[index3D(numz, numclose, i, j, m)];
							}
						}
					if (temp2 > 0) { temp += log(temp2); }// / (float)mlength);}
					else { temp += zerolog; }
					}
					temp = exp(temp / (float)mlength);
				}
				newblur[index2D(numz, i, j)] = temp;
			}
		}
	}
}

//Second Derivative of a Gaussian
float secderndis(float m, float s,float x)
{
	if (s == 0) { return 0; }
	return s/2*(2*exp(-pow(m-x,2)/s))/s-(4*exp(-pow(m-x,2)/s)*pow(m-x,2))/pow(s,2);
}


void blur_baseline(float *baseline, const int lengthmz, const float *dataMZ, const float mzsig, int mult, int filterwidth)
{
	int mulin = mult;
	float *temp = NULL;
	temp = calloc(lengthmz, sizeof(float));
	memcpy(temp, baseline, sizeof(float)*lengthmz);
	int i,j;
	#pragma omp parallel for private (i,j), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		float mzdiff = 0;
		if (i > 0) {
			mzdiff = dataMZ[i] - dataMZ[i - 1];
		}
		else {
			mzdiff = dataMZ[i + 1] - dataMZ[i];
		}
		if (mulin == 0 && mzdiff>0)
		{
			//mult = lengthmz / 500;
			mult = (int)(2 * mzsig / mzdiff);
		}
		if (mult < 1) { mult = 1; }

		float val = 0;
		int window = filterwidth;
		for (j = -window; j < window; j++)
		{
			int k = i + j*(mult);
			float newval = 0;
			if(k>=0&&k<lengthmz){
				newval = temp[k];
			}
			if (k < 0)
			{
				newval = temp[-k];
			}
			if (k >= lengthmz)
			{
				newval = temp[2 * lengthmz - k];
			}
			//if(newval>0){
			val += newval;
			//}
		}		
		baseline[i] = val/ ((float)window*2+1);
	}
	free(temp);
}

int compare_function(const void *a, const void *b)
{
	float *x = (float *)a;
	float *y = (float *)b;
	if (*x < *y) { return -1; }
	else if (*x > *y) { return 1; }
	return 0;
}

void midblur_baseline(float *baseline, const int lengthmz,const float *dataMZ, const float mzsig, int mult)
{
	if (mult == 0)
	{
		mult = lengthmz / 400;
	}
	int window = 25;

	float *temp = NULL;
	temp = calloc(lengthmz, sizeof(float));
	memcpy(temp, baseline, sizeof(float)*lengthmz);
	int i, j;
	#pragma omp parallel for private(i,j), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		float val = 0;
		float len = window * 2;
		float *med = NULL;
		med = calloc(len, sizeof(float));
		int index = 0;
		for (j = -window; j < window; j++)
		{
			int k = i + j * mult;
			float newval = 0;
			if (k >= 0 && k<lengthmz) {
				newval = temp[k];
			}
			if (k < 0)
			{
				newval = temp[-k];
			}
			if (k >= lengthmz)
			{
				newval = temp[2 * lengthmz - k];
			}
			med[index] = newval;
			index++;
		}
		qsort(med, len, sizeof(float), compare_function);
		index = 0;
		for (j = 0; j < window; j++)
		{
			val += med[j];
			index++;
		}
		if(index!=0){
			baseline[i] = val/((float)index);
		}
		else { baseline[i] = 0; }
		free(med);
	}
	free(temp);
}

void blur_noise(float *noise, int lengthmz)
{
	float *temp = NULL;
	temp = calloc(lengthmz, sizeof(float));
	memcpy(temp, noise, sizeof(float)*lengthmz);
	int i, j;
	float filter[5] = { -0.1,-0.4,1,-0.4,-0.1 };
	#pragma omp parallel for private (i,j), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		float val = 0;
		for (j = -2; j <=2 ; j++)
		{
			int k = i + j;
			float newval=0;
			if (k >= 0 && k<lengthmz) {
				newval = temp[k];
			}
			if (k < 0)
			{
				newval = temp[-k];
			}
			if (k >= lengthmz)
			{
				newval = temp[2 * lengthmz - k];
			}
			val += newval*filter[j+2];
		}
		noise[i] = val;
	}
	free(temp);
}

inline int fixk(int k, int lengthmz)
{
	k = abs(k);
	if (k >= lengthmz){ k = 2 * lengthmz - k-2; }
	//if (k < 0) { k = 0; }
	return k;
}

void convolve_simp(const int lengthmz, const int maxlength, const int*starttab, const int *endtab, const float *mzdist, const float *deltas, float *denom, const int speedyflag)
{
	if (speedyflag == 0){
		#pragma omp parallel for schedule(auto)
		for (int i = 0; i < lengthmz; i++)
		{
			float cv = 0;
			for (int k = starttab[i]; k <= endtab[i]; k++)
			{
				int k2 = fixk(k, lengthmz);
				int start = starttab[k2];
				cv += deltas[k2] * mzdist[index2D(maxlength, k2, i - start)];
			}
			denom[i] = cv;
		}
	}
	else{
		#pragma omp parallel for schedule(auto)
		for (int i = 0; i < lengthmz; i++)
		{
			float cv = 0;
			for (int k = starttab[i]; k <= endtab[i]; k++)
			{
				cv += deltas[k] * mzdist[indexmod(lengthmz,k, i)];
			}
			denom[i] = cv;
		}
	}
}


void sum_deltas(const int lengthmz, const int numz, const float * __restrict blur, const char * __restrict barr,
	const int isolength, const int *__restrict isotopepos, const float *__restrict isotopeval, float *deltas)
{
	int i,j,k;
	float temp;
	//Collapse the grid into a 1D array of delta function values
	if (isolength == 0) {
	#pragma omp parallel for private (i,j), schedule(auto)
		for (i = 0; i < lengthmz; i++)
		{
			temp = 0;
			for (j = 0; j < numz; j++)
			{
				if (barr[index2D(numz, i, j)] == 1)
				{
					temp += blur[index2D(numz, i, j)];
				}
			}
			deltas[i] = temp;
		}
	}
	else {
		//#pragma omp parallel for private (i,j,k), schedule(auto)
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				if (barr[index2D(numz, i, j)] == 1) {
					float topval = blur[index2D(numz, i, j)];
					for (k = 0; k < isolength; k++)
					{
						int pos = isotopepos[index3D(numz, isolength, i, j, k)];
						float val = isotopeval[index3D(numz, isolength, i, j, k)];
						deltas[pos] += topval* (float)val;
					}
				}
			}
		}
	}
}

void apply_ratios(const int lengthmz, const int numz, const float * __restrict blur, const char * __restrict barr,
	const int isolength, const int *__restrict isotopepos, const float *__restrict isotopeval, const float * __restrict denom, float *blur2)
{
	int i, j, k;
	if (isolength == 0)
	{
		#pragma omp parallel for private (i,j), schedule(auto)
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				if (barr[index2D(numz, i, j)] == 1)
				{
					blur2[index2D(numz, i, j)] = denom[i] * blur[index2D(numz, i, j)];
				}
				else
				{
					blur2[index2D(numz, i, j)] = 0;
				}
			}
		}
	}
	else
	{
		#pragma omp parallel for private (i,j,k), schedule(auto)
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				if (barr[index2D(numz, i, j)] == 1)
				{
					float temp = 0;
					for (k = 0; k < isolength; k++)
					{
						int pos = isotopepos[index3D(numz, isolength, i, j, k)];
						float val = isotopeval[index3D(numz, isolength, i, j, k)];
						temp += val*(denom[pos]);
					}
					blur2[index2D(numz, i, j)] = (temp)*blur[index2D(numz, i, j)];
				}
				else
				{
					blur2[index2D(numz, i, j)] = 0;
				}
			}
		}
	}
}

void deconvolve_baseline(const int lengthmz, const float *dataMZ, const float * dataInt, float * baseline, const float mzsig)
{
	float *denom = NULL;
	denom = calloc(lengthmz, sizeof(float));

	midblur_baseline(baseline, lengthmz, dataMZ,mzsig, 0);
	midblur_baseline(baseline, lengthmz, dataMZ, mzsig, 5);

	memcpy(denom, baseline, sizeof(float)*lengthmz);
	int i;
	#pragma omp parallel for private(i), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		if (denom[i] != 0 && dataInt[i] >= 0) { denom[i] = dataInt[i] / denom[i]; }
	}

	midblur_baseline(denom, lengthmz, dataMZ, mzsig, 0);
	midblur_baseline(denom, lengthmz, dataMZ, mzsig, 5);
	#pragma omp parallel for private(i), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		baseline[i] = baseline[i] * (denom[i]);
	}
	free(denom);
}

float deconvolve_iteration_speedy(const int lengthmz, const int numz,const int maxlength, const float * __restrict blur, float * __restrict blur2,
						  const char * __restrict barr, const int aggressiveflag,const float * __restrict  dataInt,
							const int isolength,const int *__restrict isotopepos,const float *__restrict isotopeval
						  ,const int*starttab, const int *endtab, const float *mzdist, const float* rmzdist, const int speedyflag,const int baselineflag, float *baseline,
							float *noise, const float mzsig, const float *dataMZ, const float filterwidth, const float psig)
{
	unsigned int i,j,k;
	float *deltas=NULL,*denom=NULL;
	deltas=calloc(lengthmz,sizeof(float));
	denom=calloc(lengthmz,sizeof(float));

	if (aggressiveflag==1 && mzsig!=0) {
		blur_baseline(baseline, lengthmz, dataMZ, fabs(mzsig), 0, filterwidth);
		//blur_baseline(baseline, lengthmz, 10);
		//blur_noise(noise, lengthmz); 
	}
	//printf("1\n");
	//Sum deltas
	sum_deltas(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, deltas);
	//printf("2\n");
	if (mzsig != 0 && psig>=0)
	{
		//Convolve with peak shape
		convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, denom, speedyflag);
	}
	else
	{
		memcpy(denom, deltas, sizeof(float)*lengthmz);
	}
	//printf("3\n");
	if (aggressiveflag == 1)
	{
		#pragma omp parallel for private(i), schedule(auto)
		for(i=0;i<lengthmz;i++){
			denom[i] += baseline[i];// +noise[i]);
		}
	}
	//printf("4\n");
	//Calculate Ratio
	#pragma omp parallel for private(i), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		if (denom[i] != 0 && dataInt[i] >= 0) { denom[i] = dataInt[i] / denom[i]; }
	}
	//printf("5\n");
	if ( mzsig < 0)
	{
		//Real Richardson-Lucy Second Convolution
		convolve_simp(lengthmz, maxlength, starttab, endtab, rmzdist, denom, deltas, speedyflag);
		memcpy(denom, deltas, sizeof(float)*lengthmz);
		/*
		for (i = 0; i<lengthmz; i++)
		{
			deltas[lengthmz - 1 - i] = denom[i];
		}
		convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, denom, speedyflag);
		for (i = 0; i<lengthmz; i++)
		{
			deltas[lengthmz - 1 - i] = denom[i];
		}
		memcpy(denom, deltas, sizeof(float)*lengthmz);*/
	}

	//printf("6\n");
	//Multiply Ratio by prior
	apply_ratios(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, denom, blur2);

	//printf("7\n");
	if (aggressiveflag ==1)
	{
		//memcpy(deltas, denom, sizeof(float)*lengthmz);
		blur_baseline(denom, lengthmz, dataMZ, fabs(mzsig), 0, filterwidth);
		//blur_baseline(denom, lengthmz, 10);
		//blur_noise(deltas, lengthmz);
		#pragma omp parallel for private(i), schedule(auto)
		for (i = 0; i < lengthmz; i++)
		{
			baseline[i] = baseline[i] * (denom[i]);
			//noise[i] = noise[i]*(deltas[i]);
		}
	}
	//printf("8\n");
 free(deltas);
 free(denom);
 return 0;
}

void ApplyCutoff1D(float *array, float cutoff, int lengthmz)
{
	unsigned int i;
	//#pragma omp parallel for private (i,j), schedule(auto)
	for (i = 0; i<lengthmz; i++)
	{
		if (array[i]<cutoff) {array[i] = 0;}
	}
}

float getfitdatspeedy(float *fitdat, const float *blur, const char *barr, const int lengthmz,const int numz,const int maxlength, const float maxint,
	const int isolength, const int *isotopepos, const float *isotopeval, const int*starttab, const int *endtab, const float *mzdist, const int speedyflag)
{
	unsigned int i,j,k;
	float *deltas=NULL;
	deltas=calloc(lengthmz,sizeof(float));
	if (isolength == 0){
		#pragma omp parallel for private (i,j), schedule(auto)
		for (i = 0; i < lengthmz; i++) //Collapse the grid into a 1D array of delta function values
		{
			float temp = 0;
			for (j = 0; j < numz; j++)
			{
				temp += blur[index2D(numz, i, j)];
			}
			deltas[i] = temp;
		}
	}
	else{
		for (i = 0; i < lengthmz; i++) //Collapse the grid into a 1D array of delta function values
		{
			for (j = 0; j < numz; j++)
			{
				float topval = blur[index2D(numz, i, j)];
				for (k = 0; k < isolength; k++)
				{
					int pos = isotopepos[index3D(numz, isolength, i, j, k)];
					float val = isotopeval[index3D(numz, isolength, i, j, k)];
					deltas[pos] += topval* (float)val;
				}
			}
		}

	}
	if (maxlength != 0)
	{
		convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, fitdat, speedyflag);
	}
	else
	{
		memcpy(fitdat, deltas, sizeof(float)*lengthmz);
	}

	free(deltas);
	float fitmax=0;
	//#pragma omp parallel for private (i), schedule(auto)
	for(i=0;i<lengthmz;i++)
	{
		if(fitdat[i]>fitmax){fitmax=fitdat[i];}
	}
	//#pragma omp parallel for private (i), schedule(auto)
	if (fitmax != 0)
	{
		for (i = 0; i<lengthmz; i++)
		{
			if (fitdat[i]<0) { fitdat[i] = 0; }
			else { fitdat[i] = fitdat[i] * maxint / fitmax; }
		}
	}
	return fitmax;
}

float errfunspeedy(Config config, Decon decon, const char *barr, const float *dataInt, const int maxlength,
	const int *isotopepos, const float *isotopeval, const int*starttab, const int *endtab, const float *mzdist, float *rsquared)
{
	//Get max intensity
	float maxint = 0;
	for (int i = 0; i < config.lengthmz; i++)
	{
		if (dataInt[i] > maxint) { maxint = dataInt[i]; }
	}

	getfitdatspeedy(decon.fitdat, decon.blur, barr, config.lengthmz, config.numz,maxlength,
		maxint, config.isolength, isotopepos, isotopeval,starttab,endtab,mzdist,config.speedyflag);
    
	if (config.baselineflag == 1)
	{
	#pragma omp parallel for schedule(auto)
		for (int i = 0; i < config.lengthmz; i++) {
			decon.fitdat[i] += decon.baseline[i];// +decon.noise[i];
			//decon.fitdat[i] = decon.noise[i]+0.1;
		}
	}
	ApplyCutoff1D(decon.fitdat, 0, config.lengthmz);

	float fitmean = Average(config.lengthmz, dataInt);

	float error=0;
	float sstot = 0;
    for(int i=0;i< config.lengthmz;i++)
    {
        error+=pow((decon.fitdat[i]-dataInt[i]),2);
		sstot += pow((dataInt[i] - fitmean), 2);
    }

	//Calculate R-squared
	if (sstot != 0) { *rsquared = 1 - (error / sstot); }

    return error;
}

//Finds nearest factor of two for optimizing FFTs
int twopow(int num){
	int n,val;
	n=0;
	val=1;
	while(n<100&&val<num){
		n++;
		val=pow(2,n);
	}
	return val;
}

//Average native charge state from Champ
float nativecharge(float mass,float fudge)
{
	return 0.0467*pow(mass,0.533)+fudge;
}



//Calculate Standard Deviation
float StdDev(int length,float *xarray,float wmean)
{
	float temp1=0;
	float temp2=0;
	int i;
	for(i=0;i<length;i++)
	{
		temp1+=pow(xarray[i]-wmean,2);
		temp2+=1;
	}
	if (temp2 == 0) { return 0; }
	return sqrt(temp1/temp2);
}

//Confidence score
float Confidence(float fit, float data)
{
	if (data == 0) { return 0; }
	return clip(1-fabs(fit-data)/data,0);
}

void KillMass(float killmass, int lengthmz, int numz, char *barr, int *nztab, float adductmass, float *dataMZ, float psfun, float mzsig)
{
	float thresh;
	if (psfun == 0){ thresh = mzsig*2.35482; }
	else{ thresh = mzsig; }
	int j, i1, i2, k;
	float testmz;
	//#pragma omp parallel for private (j,testmz,i1,i2), schedule(auto)
	for (j = 0; j < numz; j++)
	{
		testmz = (killmass + adductmass*nztab[j]) / (float)nztab[j];
		i1 = nearfast(dataMZ, testmz - thresh, lengthmz);
		i2 = nearfast(dataMZ, testmz + thresh, lengthmz);
		for (k = i1; k <= i2; k++){ barr[index2D(numz, k, j)] = 0; }
	}
}

void TestMassListWindowed(int lengthmz, int numz, char *barr, float *mtab, float nativezub, float nativezlb, float massub, float masslb, int *nztab, float *testmasses, int mfilelen, float mtabsig)
{
	unsigned int i, j;
	float testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		for (j = 0; j < numz; j++)
		{
			barr[index2D(numz, i, j)] = 0;
			testmass = mtab[index2D(numz, i, j)];
			nativelimit = nativecharge(testmass, 0);
			if (testmass<massub&&testmass>masslb&&nztab[j]<nativelimit + nativezub&&nztab[j]>nativelimit + nativezlb)
			{
				if (neartest(testmasses, testmass, mfilelen, mtabsig) == 1)
				{
					barr[index2D(numz, i, j)] = 1;
				}
			}
		}
	}
}

void TestMassListLimit(int lengthmz, int numz, char *barr, float *mtab, float nativezub, float nativezlb, float massub, float masslb, int *nztab, int *testmasspos, int mfilelen)
{
	unsigned int i, j, k;
	float testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		for (j = 0; j < numz; j++)
		{
			barr[index2D(numz, i, j)] = 0;
			testmass = mtab[index2D(numz, i, j)];
			nativelimit = nativecharge(testmass, 0);
			if (testmass<massub&&testmass>masslb&&nztab[j]<nativelimit + nativezub&&nztab[j]>nativelimit + nativezlb)
			{
				for (k = 0; k < mfilelen; k++)
				{
					if (testmasspos[index2D(numz, k, j)] == i)
					{
						barr[index2D(numz, i, j)] = 1;
					}
				}
			}
		}
	}
}

void TestMass(int lengthmz, int numz, char *barr, float *mtab, float nativezub, float nativezlb, float massub, float masslb, int *nztab)
{
	unsigned int i, j;
	float testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(auto)
	for (i = 0; i < lengthmz; i++)
	{
		for (j = 0; j < numz; j++)
		{
			testmass = mtab[index2D(numz, i, j)];
			nativelimit = nativecharge(testmass, 0);
			if (testmass<massub&&testmass>masslb&&nztab[j]<nativelimit + nativezub&&nztab[j]>nativelimit + nativezlb)
			{
				barr[index2D(numz, i, j)] = 1;
			}
			else
			{
				barr[index2D(numz, i, j)] = 0;
			}
		}
	}
}

void ManualAssign(float *dataMZ, char *barr, int *nztab, Config config)
{
	unsigned int i, j,k;
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
	float *manualmz = malloc(sizeof(float)*manlen);
	float *manualwin = malloc(sizeof(float)*manlen);
	float *manualassign = malloc(sizeof(float)*manlen);
	if (config.filetype == 1)
	{
		mh5readfile3d(config.file_id, "/config/manuallist", manlen, manualmz, manualwin, manualassign);
	}
	else
	{
		readfile3(config.manualfile, manlen, manualmz, manualwin, manualassign);
	}
	
	for( i = 0; i < manlen; i++){
		if (manualassign[i] < 0){ manualmz[i] = -manualmz[i]*1000.0; }
		//Cheating a bit...make the manualmz very negative if deassign, so that they aren't discovered by the assign loop
		}

	float testmz;
	int closest;
	#pragma omp parallel for private (i,j,testmz,closest), schedule(auto)
	for(i=0;i<lengthmz;i++)
	{
		for(j=0;j<numz;j++)
		{
			//Manual Assign: The nearest value wins in the case of overlap
			testmz=dataMZ[i];
			closest=nearunsorted(manualmz,testmz,manlen);
			if (fabs(manualmz[closest] - testmz)<manualwin[closest] && (float)nztab[j] == manualassign[closest] && manualassign[closest]>0)
				{
				barr[index2D(numz, i, j)]=1;
				}
			else if (fabs(manualmz[closest] - testmz)<manualwin[closest] && (float)nztab[j] != manualassign[closest] && manualassign[closest]>0)
				{
				barr[index2D(numz, i, j)]=0;
				}
			//Manual Deassign: Anything within the window is killed
			for (k = 0; k < manlen; k++){
				
				if (fabs(manualmz[k] - testmz*(-1000.0))<manualwin[k]*1000.0 && (float)nztab[j] == fabs(manualassign[k]) && manualassign[k]<0)
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

/*
void MakeBlur(const int lengthmz, const int numz, const int numclose, char *barr, const int *closezind,
	const int *closemind, const float *mtab, const float molig, const float adductmass, const int *nztab, 
	const float *dataMZ,int *closeind, const float threshold, const Config config)
{
	
	  #pragma omp parallel for schedule(auto)
	  for(int i=0;i<lengthmz;i++)
	    {
	      for(int j=0;j<numz;j++)
	        {
	            if(barr[index2D(numz,i,j)]==1)
	            {
					//Reset the threshold if it is zero
					float newthreshold = threshold;
					if (newthreshold == 0) { newthreshold = config.massbins * 3 / (float)nztab[j];}
				
					int num = 0;
					for(int k=0;k<numclose;k++)
					{
						//Find the z index and test if it is within the charge range
						int indz=(int)(j+closezind[k]);
						if (indz < 0 || indz >= numz || (nztab[j] + closezind[k])==0) { closeind[index3D(numz, numclose, i, j, k)] = -1; }
						else
						{
							//Find the nearest m/z value and test if it is close enough to the predicted one and within appropriate ranges
							float point=(float)((mtab[index2D(numz, i, j)]+closemind[k]*molig+adductmass*(float)(nztab[j]+closezind[k]))/(float)(nztab[j]+closezind[k]));
							if(point<dataMZ[0]- newthreshold || point>dataMZ[lengthmz-1]+ newthreshold){ closeind[index3D(numz, numclose, i, j, k)] = -1; }
							else {
								int ind = nearfast(dataMZ, point, lengthmz);
								float closepoint = dataMZ[ind];
								int newind= index2D(numz, ind, indz);
								if (barr[newind] == 1 && fabs(point- closepoint)< newthreshold) {
									closeind[index3D(numz, numclose, i, j, k)] = newind;
									num += 1;
								}
								else { closeind[index3D(numz, numclose, i, j, k)] = -1; }
								//printf("%d %d %d %f %f %d\n", i, j, k, point, closepoint, closeind[index3D(numz, numclose, i, j, k)]);
							}
						  
						}
					}
					if (num < 2 && config.isotopemode == 0) { barr[index2D(numz, i, j)] = 0; }// printf("%d %d \n", i, j);}
				}
				else
				{
					for (int k=0;k<numclose;k++)
					{
						closeind[index3D(numz, numclose, i, j, k)] = -1;
					}
				}
	        }
	    }
}*/

void MakeSparseBlur(const int numclose, char* barr, const int* closezind,
	const int* closemind, const float* mtab, const int* nztab,
	const float* dataMZ, int* closeind, float *closeval, float *closearray, const Config config)
{
	int lengthmz = config.lengthmz;
	int numz = config.numz;
	float molig = config.molig;
	float adductmass = config.adductmass;

	#pragma omp parallel for schedule(auto)
	for (int i = 0; i < lengthmz; i++)
	{
		for (int j = 0; j < numz; j++)
		{
			if (barr[index2D(numz, i, j)] == 1)
			{
				int num = 0;

				//Reset the threshold if it is zero
				float mzsig = config.mzsig;
				if (mzsig == 0) { 
					int i1 = i - 1;
					int i2 = i + 1;
					if (i >= lengthmz - 1) { i2 = i; }
					if (i == 0) { i1 = i; }
					mzsig = 2*fabs(dataMZ[i2] - dataMZ[i1]);
					if (mzsig > config.massbins || mzsig == 0) { mzsig = config.massbins*2; }
				}
				float newthreshold = mzsig * 2;

				for (int k = 0; k < numclose; k++)
				{
					//Find the z index and test if it is within the charge range
					int indz = (int)(j + closezind[k]);
					if (indz < 0 || indz >= numz || (nztab[j] + closezind[k]) == 0) 
					{ 
						closeind[index3D(numz, numclose, i, j, k)] = -1; 
						closearray[index3D(numz, numclose, i, j, k)] = 0;
					}
					else
					{
						//Find the nearest m/z value and test if it is close enough to the predicted one and within appropriate ranges
						float point = (float)((mtab[index2D(numz, i, j)] + closemind[k] * molig + adductmass * (float)(nztab[j] + closezind[k])) / (float)(nztab[j] + closezind[k]));
						if (point<dataMZ[0] - newthreshold || point>dataMZ[lengthmz - 1] + newthreshold) 
						{ 
							closeind[index3D(numz, numclose, i, j, k)] = -1;
							closearray[index3D(numz, numclose, i, j, k)] = 0;
						}
						else {
							int ind = nearfast(dataMZ, point, lengthmz);
							float closepoint = dataMZ[ind];
							int newind = index2D(numz, ind, indz);
							if (barr[newind] == 1 && fabs(point - closepoint) < newthreshold) {
								closeind[index3D(numz, numclose, i, j, k)] = newind;
								closearray[index3D(numz, numclose, i, j, k)] = closeval[k] *mzpeakshape(point, closepoint, mzsig, config.psfun);
								num += 1;
							}
							else { 
								closeind[index3D(numz, numclose, i, j, k)] = -1;
								closearray[index3D(numz, numclose, i, j, k)] = 0;
							}
							//printf("%d %d %d %f %f %d\n", i, j, k, point, closepoint, closeind[index3D(numz, numclose, i, j, k)]);
						}

					}
				}
				if (num < 2 && config.isotopemode == 0) { barr[index2D(numz, i, j)] = 0; }// printf("%d %d \n", i, j);}
			}
			else
			{
				for (int k = 0; k < numclose; k++)
				{
					closeind[index3D(numz, numclose, i, j, k)] = -1;
					closearray[index3D(numz, numclose, i, j, k)] = 0;
				}
			}
		}
	}
}

/*
void MakeBlurZ(int lengthmz, int numz, int numclose, char *barr, int *closezind, float *mtab, float adductmass, int *nztab, float *dataMZ, int *closeind)
{
	unsigned int i, j, k;
	#pragma omp parallel for private (i,j,k), schedule(auto)
	for (i = 0; i<lengthmz; i++)
	{
		for (j = 0; j<numz; j++)
		{
			if (barr[index2D(numz, i, j)] == 1)
			{
				for (k = 0; k<numclose; k++)
				{
					int indz = (int)(j + closezind[k]);
					int ind;
					if (indz < 0 || indz >= numz || (nztab[j] + closezind[k]) == 0) { closeind[index3D(numz, numclose, i, j, k)] = -1; }
					else
					{
						float point = (float)((mtab[index2D(numz, i, j)] + adductmass*(float)(nztab[j] + closezind[k])) / (float)(nztab[j] + closezind[k]));
						if (point<dataMZ[0] || point>dataMZ[lengthmz - 1]) { closeind[index3D(numz, numclose, i, j, k)] = -1; }
						else {
							ind = nearfast(dataMZ, point, lengthmz);
							closeind[index3D(numz, numclose, i, j, k)] = index2D(numz, ind, indz);
						}
					}
				}
			}
			else
			{
				for (k = 0; k<numclose; k++)
				{
					closeind[index3D(numz, numclose, i, j, k)] = -1;
				}
			}
		}
	}
}


void MakeBlurM(int lengthmz, int numz, int numclose, char *barr, int *closemind, float *mtab, float molig, float adductmass, int *nztab, float *dataMZ, int *closeind)
{
	unsigned int i, j, k;
#pragma omp parallel for private (i,j,k), schedule(auto)
	for (i = 0; i<lengthmz; i++)
	{
		for (j = 0; j<numz; j++)
		{
			if (barr[index2D(numz, i, j)] == 1)
			{
				for (k = 0; k<numclose; k++)
				{
					int ind;
					float point = (float)((mtab[index2D(numz, i, j)] + closemind[k] * molig + adductmass*(float)(nztab[j])) / (float)(nztab[j]));
					if (point<dataMZ[0] || point>dataMZ[lengthmz - 1]) { closeind[index3D(numz, numclose, i, j, k)] = -1; }
					else {
						ind = nearfast(dataMZ, point, lengthmz);
						closeind[index3D(numz, numclose, i, j, k)] = index2D(numz, ind, j);
					}
				}
			}
			else
			{
				for (k = 0; k<numclose; k++)
				{
					closeind[index3D(numz, numclose, i, j, k)] = -1;
				}
			}
		}
	}
}*/


void MakePeakShape2D(int lengthmz,int maxlength,int *starttab,int *endtab,float *dataMZ,float mzsig,int psfun,int speedyflag,float *mzdist, float *rmzdist, int makereverse)
{
	#pragma omp parallel for schedule(auto)
	for(int i=0;i<lengthmz;i++)
	  {
		  int start = starttab[i];
		  int end = endtab[i];
		  for (int j = start; j <= end; j++)
		  {	
			  int j2 = fixk(j, lengthmz);
			  mzdist[index2D(maxlength, i, j2 - start)] = mzpeakshape(dataMZ[i], dataMZ[j2], mzsig, psfun);
			  if (makereverse == 1) { rmzdist[index2D(maxlength, i, j2 - start)] = mzpeakshape(dataMZ[j2], dataMZ[i], mzsig, psfun); }
		  }
	  }
}

void MakePeakShape1D(float *dataMZ,float threshold,int lengthmz,int speedyflag,float mzsig,int psfun,float *mzdist, float *rmzdist, int makereverse)
{
	float binsize=dataMZ[1]-dataMZ[0];
	float newrange=threshold/binsize;
	int n;
	for(n=(int)-newrange;n<(int)newrange;n++)
	{
		mzdist[indexmod(lengthmz,0,n)]=mzpeakshape(0,n*binsize,mzsig,psfun);
		if (makereverse == 1) { rmzdist[indexmod(lengthmz, 0, n)] = mzpeakshape(n * binsize, 0, mzsig, psfun); }
	}
	printf("\nNotice: Assuming linearized data. \n\n");
}


void monotopic_to_average(const int lengthmz, const int numz, float *blur, const char *barr, int isolength, const int *__restrict isotopepos, const float *__restrict isotopeval)
{
	float *newblur = NULL;
	newblur = calloc(lengthmz*numz, sizeof(float));
	unsigned int i, j, k;
	for (i = 0; i < lengthmz; i++)
	{
		for (j = 0; j < numz; j++)
		{
			if (barr[index2D(numz, i, j)] == 1) {
				float topval = blur[index2D(numz, i, j)];
				for (k = 0; k < isolength; k++)
				{
					int pos = isotopepos[index3D(numz, isolength, i, j, k)];
					float val = isotopeval[index3D(numz, isolength, i, j, k)];
					newblur[index2D(numz, pos, j)] += topval* (float)val;
				}
			}
		}
	}
	memcpy(blur, newblur, sizeof(float)*lengthmz*numz);
	free(newblur);
}


float Reconvolve(const int lengthmz, const int numz, const int maxlength,const int *starttab,const int *endtab,const float *mzdist, const float *blur, float *newblur, const int speedyflag, const char *barr)
{
	float newblurmax=0;
	if (speedyflag == 0){
		#pragma omp parallel for schedule(auto)
		for (int i = 0; i < lengthmz; i++)
		{
			for (int j = 0; j < numz; j++)
			{

				float cv = 0;
				if (barr[index2D(numz, i, j)] == 1) {
					for (int k = starttab[i]; k <= endtab[i]; k++)
					{
						int k2 = fixk(k, lengthmz);
						if (blur[index2D(numz, k2, j)] != 0)
						{
							int start = starttab[k2];
							cv += blur[index2D(numz, k2, j)] * mzdist[index2D(maxlength, k2, i - start)];
						}
					}
				}
				newblur[index2D(numz, i, j)] = cv;
				if (cv > newblurmax)
				{
					newblurmax = cv;
				}
			}
		}
	}
	else{
		#pragma omp parallel for schedule(auto)
		for (int i = 0; i<lengthmz; i++)
		{
			for (int j = 0; j<numz; j++)
			{

				float cv = 0;
				if (barr[index2D(numz, i, j)] == 1) {
					for (int k = starttab[i]; k <= endtab[i]; k++)
					{
						if (blur[index2D(numz, k, j)] != 0)
						{
							cv += blur[index2D(numz, k, j)] * mzdist[indexmod(lengthmz, k, i)];
						}
					}
				}
			newblur[index2D(numz, i, j)] = cv;
			if (cv>newblurmax)
			{
				newblurmax = cv;
			}
			}
		}

	}
	return newblurmax;
}


void IntegrateTransform(const int lengthmz, const int numz, const float *mtab, float massmax, float massmin,
	const int maaxle, float *massaxis,float *massaxisval, const float*blur,float *massgrid)
{
	for(int i=0;i<lengthmz;i++)
		{
			for(int j=0;j<numz;j++)
			{
			float testmass=mtab[index2D(numz,i,j)];
			if(testmass<massmax&&testmass>massmin){
				int index=nearfast(massaxis,testmass,maaxle);
				float newval=blur[index2D(numz,i,j)];
				if (massaxis[index] == testmass){
					massaxisval[index]+=newval;
					massgrid[index2D(numz,index,j)]+=newval;
				}
				
				if (massaxis[index] < testmass &&index<maaxle-2)
				{
					int index2 = index + 1;
					float interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
					massaxisval[index] += (1.0 - interpos)*newval;
					massgrid[index2D(numz, index, j)] += (1.0 - interpos)*newval;
					massaxisval[index2] += (interpos)*newval;
					massgrid[index2D(numz, index2, j)] += (interpos)*newval;
				}
				
				if (massaxis[index] > testmass&&index>0)
				{
					int index2 = index - 1;
					float interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
					massaxisval[index] += (1 - interpos)*newval;
					massgrid[index2D(numz, index, j)] += (1 - interpos)*newval;
					massaxisval[index2] += (interpos)*newval;
					massgrid[index2D(numz, index2, j)] += (interpos)*newval;
				}			
			}
			}
		}
}

void InterpolateTransform(const int maaxle,const int numz,const int lengthmz, const int *nztab,float *massaxis,
	const float adductmass,const float *dataMZ, float *massgrid, float *massaxisval, const float *blur)
{
	float startmzval = dataMZ[0];
	float endmzval = dataMZ[lengthmz - 1];
	//#pragma omp parallel for schedule(auto)
	for(int i=0;i<maaxle;i++)
	{
		float val=0;
		
		for(int j=0;j<numz;j++)
		{
			float newval = 0;
			float mztest=(massaxis[i]+nztab[j]*adductmass)/(float)nztab[j];

			if(mztest>startmzval&&mztest<endmzval)
			{
				int index=nearfast(dataMZ,mztest,lengthmz);
				int index2=index;
				if(dataMZ[index]==mztest)
				{
					newval=blur[index2D(numz,index,j)];
					val+=newval;
					massgrid[index2D(numz,i,j)]=newval;
				}
				else
				{
					if(dataMZ[index]>mztest&&index>1&&index<lengthmz-1)
					{
						index2=index;
						index=index-1;
					}
					else if(dataMZ[index]<mztest&&index<lengthmz-2&&index>0)
					{
						index2=index+1;
					}
					if(index2>index&&(dataMZ[index2]-dataMZ[index])!=0)
					{
						float mu=(mztest-dataMZ[index])/(dataMZ[index2]-dataMZ[index]);
						float y0=blur[index2D(numz,index-1,j)];
						float y1=blur[index2D(numz,index,j)];
						float y2=blur[index2D(numz,index2,j)];
						float y3=blur[index2D(numz,index2+1,j)];
						newval=clip(CubicInterpolate(y0,y1,y2,y3,mu),0);
						//newval=CRSplineInterpolate(y0,y1,y2,y3,mu);
						//newval=LinearInterpolate(y1,y2,mu);
						val+=newval;
						massgrid[index2D(numz,i,j)]=newval;
						}
					}
				}
			}
		massaxisval[i]=val;
	}
}

void SmartTransform(const int maaxle, const int numz, const int lengthmz, const int* nztab, float* massaxis,
	const float adductmass, const float* dataMZ, float* massgrid, float* massaxisval, const float* blur)
{
	float startmzval = dataMZ[0];
	float endmzval = dataMZ[lengthmz - 1];
	#pragma omp parallel for schedule(auto)
	for (int i = 0; i < maaxle; i++)
	{
		float val = 0;
		for (int j = 0; j < numz; j++)
		{
			float z = (float) nztab[j];
			float mtest = massaxis[i];
			float mztest = (mtest + z * adductmass) / z;

			float mzlower;
			float mlower;
			if (i>0){
				mlower = massaxis[i - 1];
				mzlower = (mlower + z * adductmass) / z;
				//mzlower = (mzlower + mztest) / 2;
			}
			else { mzlower = mztest; mlower = mtest; }

			float mzupper;
			float mupper;
			if (i < maaxle-1) {
				mupper = massaxis[i + 1];
				mzupper = (mupper + z * adductmass) / z; 
				//mzupper = (mzupper + mztest) / 2;
			}
			else { mzupper = mztest; mupper = mtest; }

			float newval = 0;
			if ( mzupper > startmzval&& mzlower < endmzval)
			{
				int index = nearfast(dataMZ, mztest, lengthmz);
				int index1 = nearfast(dataMZ, mzlower, lengthmz);
				int index2 = nearfast(dataMZ, mzupper, lengthmz);
				float imz = dataMZ[index];

				if (index2-index1 < 5) { 
					
					if (imz == mztest)
					{
						newval = clip(blur[index2D(numz, index, j)],0);
						val += newval;
						massgrid[index2D(numz, i, j)] = newval;
						
					}
					else
					{
						int edge = 0;
						index2 = index;
						if (imz > mztest)
						{
							index = index - 1;
						}
						else if (imz < mztest)
						{
							index2 = index + 1;
						}
						
						if (index < 1 || index2 >= lengthmz - 1)
						{
							edge = 1;
						}
						if (index < 0 || index2 >= lengthmz)
						{
							edge = 2;
						}

						
						if (index2 > index && (dataMZ[index2] - dataMZ[index]) != 0 && edge ==0)
						{
							float mu = (mztest - dataMZ[index]) / (dataMZ[index2] - dataMZ[index]);
							float y0 = blur[index2D(numz, index - 1, j)];
							float y1 = blur[index2D(numz, index, j)];
							float y2 = blur[index2D(numz, index2, j)];
							float y3 = blur[index2D(numz, index2 + 1, j)];
							newval = clip(CubicInterpolate(y0, y1, y2, y3, mu), 0);
							//newval=clip(CRSplineInterpolate(y0,y1,y2,y3,mu), 0);
							//newval=clip(LinearInterpolate(y1,y2,mu),0);
							val += newval;
							massgrid[index2D(numz, i, j)] = newval;
							//printf("0\n");
						}
						else if (edge == 1 && (dataMZ[index2] - dataMZ[index]) != 0)
						{
							float mu = (mztest - dataMZ[index]) / (dataMZ[index2] - dataMZ[index]);
							float y1 = blur[index2D(numz, index, j)];
							float y2 = blur[index2D(numz, index2, j)];
							newval=clip(LinearInterpolate(y1,y2,mu),0);
							val += newval;
							massgrid[index2D(numz, i, j)] = newval;
							//printf("1 %d %d %f %f %f\n", index, index2, mztest, imz, newval, massaxis[i]);
						}
						else if (edge == 2 && (dataMZ[index2] - dataMZ[index]) != 0)
						{
							float factor = 1;
							if (index2 == 0) { index = 0; index2 = 1; }
							if (index == lengthmz - 1) { index = lengthmz - 1; index2 = lengthmz - 2; factor = -1; }
							float mu = (mztest - dataMZ[index]) / (dataMZ[index] - dataMZ[index2]);
							float y1 = blur[index2D(numz, index, j)];
							float y2 = 0;
							newval = clip(LinearInterpolate(y1, y2, mu), 0);
							val += newval;
							massgrid[index2D(numz, i, j)] = newval;
							//printf("2\n");
						}
					}
				}

				else{
					//printf("Integrating\n");
					float num = 0;
					for (int k = index1; k < index2+1; k++)
					{
						float kmz = dataMZ[k];
						float km = (kmz - adductmass) * z;
						float scale = 1;
						
						if (mztest < kmz && km < mupper)
						{
							//scale = LinearInterpolatePosition(mzupper, mztest, kmz);
							scale = LinearInterpolatePosition(mupper, mtest, km);
						}
						else if (kmz < mztest && km > mlower)
						{
							//scale = LinearInterpolatePosition(mzlower, mztest, kmz);
							scale = LinearInterpolatePosition(mlower, mtest, km);
						}
						else if (kmz == mztest) { scale = 1; }
						else {scale = 0;}

						newval += scale * blur[index2D(numz, k, j)];
						num += scale;

					}
					if (num != 0) { newval /= num; }
					newval = clip(newval, 0);
					val += newval;
					massgrid[index2D(numz, i, j)] = newval;
				}
			}
		}
		massaxisval[i] = val;
	}
}




void ignorezeros(char *barr, const float * dataInt, const int lengthmz, const int numz)
{
#pragma omp parallel for schedule(auto)
	for (int i = 0; i < lengthmz; i++)
	{
		float val = dataInt[i];
		if (val == 0) {
			for (int j = 0; j < numz; j++)
			{
				barr[index2D(numz, i, j)] = 0;
			}
		}
	}
}


void KillB(float *I, char *B, float intthresh, int lengthmz, int numz, const int isolength, int *isotopepos, float *isotopeval)
{
	unsigned int i, j, k;
	if (isolength == 0){
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				if (I[i] <= intthresh){ B[index2D(numz, i, j)] = 0; }
			}
		}
	}
	else
	{
		float cutoff = 0.5;
		printf("Removing species where isotope fall below %f\n", cutoff * 100);
		#pragma omp parallel for private (i,j,k), schedule(auto)
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				float max = 0;
				for (k = 0; k < isolength; k++)
				{
					float val = isotopeval[index3D(numz, isolength, i, j, k)];
					if (val > max){ max = val; }
				}
				for (k = 0; k < isolength; k++)
				{
					float val = isotopeval[index3D(numz, isolength, i, j, k)];
					if (val>cutoff*max){
						int pos = isotopepos[index3D(numz, isolength, i, j, k)];
						if (I[pos] <= intthresh){ B[index2D(numz, i, j)] = 0; }
					}
				}

			}
		}
	}
}


float isotopemid(float mass,float *isoparams)
{
float a,b,c;
a=isoparams[4];
b=isoparams[5];
c=isoparams[6];
return a+b*pow(mass,c);
}

float isotopesig(float mass,float *isoparams)
{
float a,b,c;
a=isoparams[7];
b=isoparams[8];
c=isoparams[9];
return a+b*pow(mass,c);
}

float isotopealpha(float mass,float *isoparams)
{
float a,b;
a=isoparams[0];
b=isoparams[1];
return a*exp(-mass*b);
}

float isotopebeta(float mass,float *isoparams)
{
float a,b;
a=isoparams[2];
b=isoparams[3];
return a*exp(-mass*b);
}

int setup_isotopes(float *isoparams,int *isotopepos,float *isotopeval,float *mtab,int *ztab,char *barr,float *dataMZ,int lengthmz,int numz)
{
    float minmass=100000000;
    float maxmass=1;
    int i;
    for(i=0;i<lengthmz*numz;i++)
    {
        if(barr[i]==1)
        {
			float mass = mtab[i];
            if(mass<minmass){minmass=mass;}
            if(mass>maxmass){maxmass=mass;}
        }
    }

    float minmid=isotopemid(minmass,isoparams);
    float minsig=isotopesig(minmass,isoparams);
    float maxmid=isotopemid(maxmass,isoparams);
    float maxsig=isotopesig(maxmass,isoparams);

    int isostart=(int)(minmid-4*minsig);
    int isoend=(int)(maxmid+4*maxsig);
    if(isostart<0){isostart=0;}
	if (isoend < 4){ isoend = 4; }
    int isolength=isoend-isostart;
	return isolength;
}

void make_isotopes(float *isoparams, int *isotopepos, float *isotopeval, float *mtab, int *ztab, char *barr, float *dataMZ, int lengthmz, int numz)
{
	float minmass = 100000000;
	float maxmass = 1;
	int i, j, k;
	for (i = 0; i<lengthmz*numz; i++)
	{
		if (barr[i] == 1)
		{
			float mass = mtab[i];
			if (mass<minmass){ minmass = mass; }
			if (mass>maxmass){ maxmass = mass; }
		}
	}
	float massdiff = 1.0026;

	float minmid = isotopemid(minmass, isoparams);
	float minsig = isotopesig(minmass, isoparams);
	float maxmid = isotopemid(maxmass, isoparams);
	float maxsig = isotopesig(maxmass, isoparams);

	//int isostart = (int)(minmid - 4 * minsig);
	int isostart = 0;
	int isoend = (int)(maxmid + 4 * maxsig);
	//if (isostart<0){ isostart = 0; }
	if (isoend < 4){ isoend = 4; }
	int isolength = isoend - isostart;
	float *isorange = NULL;
	int *isoindex = NULL;
	isorange = calloc(isolength, sizeof(float));
	isoindex = calloc(isolength, sizeof(int));
	for (i = 0; i<isolength; i++)
	{
		isorange[i] = (isostart + i)*massdiff;
		isoindex[i] = (isostart + i);
	}
	#pragma omp parallel for private (i,j,k), schedule(auto)
	for (i = 0; i<lengthmz; i++)
	{
		for (j = 0; j<numz; j++)
		{
			if (barr[index2D(numz, i, j)] == 1)
			{
				float mz = dataMZ[i];
				int z = ztab[j];
				for (k = 0; k<isolength; k++)
				{
					float newmz = mz + (isorange[k] / ((float)z));
					int pos = nearfast(dataMZ, newmz, lengthmz);
					isotopepos[index3D(numz, isolength, i, j, k)] = pos;
				}
			}
		}
	}
	#pragma omp parallel for private (i,j,k), schedule(auto)
	for (i = 0; i<lengthmz; i++)
	{
		for (j = 0; j<numz; j++)
		{
			if (barr[index2D(numz, i, j)] == 1)
			{
				float mass = mtab[index2D(numz, i, j)];
				float mid = isotopemid(mass, isoparams);
				float sig = isotopesig(mass, isoparams);
				if (sig == 0) { printf("Error: Sigma Isotope Parameter is 0"); exit(102); }
				float alpha = isotopealpha(mass, isoparams);
				float amp = (1.0 - alpha) / (sig*2.50662827);
				float beta = isotopebeta(mass, isoparams);
				float tot = 0;
				for (k = 0; k<isolength; k++)
				{
					float e = alpha*exp(-isoindex[k]*beta);
					float g = amp *exp(-pow(isoindex[k] - mid, 2) / (2 * pow(sig, 2)));
					float temp = e + g;
					tot+=temp;
					isotopeval[index3D(numz, isolength, i, j, k)] =temp;
				}
				for (k = 0; k<isolength; k++)
				{
					if (tot > 0){ isotopeval[index3D(numz, isolength, i, j, k)] = isotopeval[index3D(numz, isolength, i, j, k)] / tot; }
				}
			}
		}
	}
	free(isorange);
	free(isoindex);
}

void isotope_dist(float mass,int isolength, int *isoindex, float *isovals, float *isoparams)
{
	float mid = isotopemid(mass, isoparams);
	float sig = isotopesig(mass, isoparams);
	if (sig == 0) { printf("Error: Sigma Isotope Parameter is 0"); exit(102); }
	float alpha = isotopealpha(mass, isoparams);
	float amp = 1.0 - alpha;
	float beta = isotopebeta(mass, isoparams);
	float tot = 0;
	int k;
	for (k = 0; k<isolength; k++)
	{
		float e = alpha*exp(-isoindex[k] * beta);
		float g = amp / (sig*2.50662827)*exp(-pow(isoindex[k] - mid, 2) / (2 * pow(sig, 2)));
		tot += e + g;
		isovals[k] = e + g;
	}
	for (k = 0; k<isolength; k++)
	{
		if (tot > 0){ isovals[k] = isovals[k] / tot; }
	}
}

void textvectorprint(float *arr, int length)
{
	int levels = 20;
	float grad = 1.0 / ((float)levels);
	int i, j;
	printf("\n");
	float max = 0;
	for (i = 0; i < length; i++){
		if (arr[i]>max){ max = arr[i]; }
	}
	if (max != 0){
		for (i = 0; i < length; i++){
			arr[i]=arr[i] / max;
		}
	}
	for (i = 0; i < levels; i++){
		for (j = 0; j < length; j++)
		{
			if (arr[j]>grad*(levels - i )){ printf("| "); }
			else{ printf("  "); }
		}
		printf("\n");
	}
}

void test_isotopes(float mass, float *isoparams)
{
	int i;
	float maxmid = isotopemid(mass, isoparams);
	float maxsig = isotopesig(mass, isoparams);

	int isostart =0;
	int isoend = (int)(maxmid + 4 * maxsig);
	if (isoend < 4){ isoend = 4; }
	int isolength = isoend - isostart;
	float *isorange = NULL;
	int *isoindex = NULL;
	isorange = calloc(isolength, sizeof(float));
	isoindex = calloc(isolength, sizeof(int));
	for (i = 0; i<isolength; i++)
	{
		isoindex[i] = (isostart + i);
	}
	isotope_dist(mass, isolength, isoindex, isorange, isoparams);
	printf("Distribution: %f", mass);
	for (i = 0; i < isolength; i++)
	{
		printf("%f ", isorange[i]);
	}
	printf("\n");
	textvectorprint(isorange, isolength);
	free(isorange);
	free(isoindex);
}

void setup_and_make_isotopes(Config *config, Input *inp) {
	printf("Isotope Mode: %d\n", config->isotopemode);

	config->isolength = setup_isotopes(inp->isoparams, inp->isotopepos, inp->isotopeval, inp->mtab, inp->nztab, inp->barr, inp->dataMZ, config->lengthmz, config->numz);

	inp->isotopepos = calloc(config->isolength * config->lengthmz * config->numz, sizeof(int));
	inp->isotopeval = calloc(config->isolength * config->lengthmz * config->numz, sizeof(float));

	make_isotopes(inp->isoparams, inp->isotopepos, inp->isotopeval, inp->mtab, inp->nztab, inp->barr, inp->dataMZ, config->lengthmz, config->numz);

	printf("Isotopes set up, Length: %d\n", config->isolength);
}

void deepcopy(float *out, const float *in, int len)
{
	for (int i = 0; i < len; i++)
	{
		out[i] = in[i];
	}
}

void simp_norm(const int length, float *data)
{
	float norm = Max(data, length);
	//printf("\nMax: %f %d\n", norm, length);
	if(norm>0){
		for (int i = 0; i<length; i++)
		{
			data[i] = data[i] / norm;
		}
	}
	return;
}

void simp_norm_sum(const int length, float *data)
{
	float norm = Sum(data, length);
	//printf("\nMax: %f %d\n", norm, length);
	if (norm>0) {
		for (int i = 0; i<length; i++)
		{
			data[i] = data[i] / norm;
		}
	}
	return;
}



void charge_scaling(float *blur, const int *nztab, const int lengthmz, const int numz)
{
	for (int i = 0; i < lengthmz; i++)
	{
		for (int j = 0; j < numz; j++)
		{
			int charge = nztab[j];
			float z = (float)charge;
			if (z != 0) { blur[index2D(numz, i, j)] /= z; }
		}
	}
	return; 
}


void point_smoothing(float *blur, const char * barr, const int lengthmz, const int numz, const int width)
{
	float *newblur;
	newblur = calloc(lengthmz*numz, sizeof(float));
	memcpy(newblur, blur, lengthmz*numz * sizeof(float));
	#pragma omp parallel for schedule(auto)
	for (int i = 0; i < lengthmz; i++)
	{
		for (int j = 0; j < numz; j++)
		{
			if( barr[index2D(numz, i, j)]==1){
				int low = i - width;
				if (low < 0) { low = 0; }
				int high = i + width+1;
				if (high > lengthmz) { high = lengthmz; }

				float sum = 0;
				for (int k = low; k < high; k++)
				{
					sum += newblur[index2D(numz, k, j)];
				}

				blur[index2D(numz, i, j)] = sum / ((float)1 + 2 * width);
			}
		}
	}
	free(newblur);
	return;
}

void point_smoothing_peak_width(const int lengthmz, const int numz, const int maxlength, const int* starttab, const int* endtab, const float* mzdist, float* blur, const int speedyflag, const char *barr)
{
	float* newblur;
	newblur = calloc(lengthmz * numz, sizeof(float));
	memcpy(newblur, blur, lengthmz * numz * sizeof(float));
	Reconvolve(lengthmz, numz, maxlength, starttab, endtab, mzdist, newblur, blur, speedyflag, barr);
	return;
}

/*
void point_smoothing_sg(float *blur, const int lengthmz, const int numz, const int width)
{
	float *newblur;
	newblur = calloc(lengthmz*numz, sizeof(float));
	memcpy(newblur, blur, lengthmz*numz * sizeof(float));

	float *c;
	int np = width * 2 + 1;
	c = calloc(np, sizeof(float));
	savgol(c, np, width, width, 0, 2);

	//printf("Savgol: ");
	//for (int i = 0; i < np; i++){printf("%f ", c[i]);}
	//printf("\n");

	#pragma omp parallel for schedule(auto)
	for (int i = 0; i < lengthmz; i++)
	{
		for (int j = 0; j < numz; j++)
		{
			int low = i - width;
			if (low < 0) { low = 0; }
			int high = i + width + 1;
			if (high > lengthmz) { high = lengthmz; }

			float sum = 0;
			for (int k = low; k < high; k++)
			{
				int cind = mod((i - k) , np);
				sum += c[cind] * newblur[index2D(numz, k, j)];
			}

			blur[index2D(numz, i, j)] = sum;
		}
	}
	free(newblur);
	return;
}
*/

void softargmax_transposed(float *blur, const int lengthmz, const int numz, const float beta, const char *barr, const int maxlength, 
	const int isolength, const int* __restrict isotopepos, const float* __restrict isotopeval, const int speedyflag
	, int* starttab, int* endtab, float* mzdist, const float mzsig)
{
	float *newblur, *deltas, *deltas2, *denom, *denom2;
	newblur = calloc(lengthmz*numz, sizeof(float));
	deltas = calloc(lengthmz, sizeof(float));
	deltas2 = calloc(lengthmz, sizeof(float));
	denom = calloc(lengthmz, sizeof(float));
	denom2 = calloc(lengthmz, sizeof(float));
	memcpy(newblur, blur, lengthmz*numz * sizeof(float));

	//Sum deltas
	sum_deltas(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, deltas);

	if (mzsig != 0)//Convolve with peak shape
	{convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, denom, speedyflag);}
	else{memcpy(denom, deltas, sizeof(float) * lengthmz);}

	#pragma omp parallel for schedule(auto)
	for (int i = 0; i < lengthmz * numz; i++)
	{
		blur[i] = exp(beta * newblur[i])-1;
	}

	//Sum deltas
	sum_deltas(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, deltas2);

	if (mzsig != 0)//Convolve with peak shape
	{	convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas2, denom2, speedyflag);}
	else	{memcpy(denom2, deltas, sizeof(float) * lengthmz);}

	#pragma omp parallel for schedule(auto)
	for (int i = 0; i < lengthmz; i++)
	{
		float factor = 0;
		if (denom2[i] != 0) { factor = denom[i] / denom2[i]; };
		for (int j = 0; j < numz; j++)
		{
			//blur[index2D(numz, i, j)] -= 1;
			blur[index2D(numz, i, j)] *= factor;
		}
		//printf("%f %f\n", sum1, sum3);
	}
	free(newblur);
	free(deltas);
	free(deltas2);
	free(denom);
	free(denom2);
	return;
}

void softargmax_everything(float *blur, const int lengthmz, const int numz, const float beta)
{
	float *newblur;
	newblur = calloc(lengthmz*numz, sizeof(float));
	memcpy(newblur, blur, lengthmz*numz * sizeof(float));
	//float max1 = 0;
	//float max2 = 0;
	float sum2 = 0;
	float sum1 = 0;
	float min2 = 100000000000000000;
	//#pragma omp parallel for schedule(auto), reduction(min:min2), reduction(+:sum1), reduction(+:sum2)
	for (int i = 0; i < lengthmz*numz; i++)
	{
		float d = newblur[i];
		//if (d > max1) { max1 = d; }
		float e = exp(beta * d);
		//float e = pow(d, -beta);
		//if (e > max2) { max2 = e; }
		blur[i] = e;
		if (e < min2) { min2 = e; }
		sum1 += d;
		sum2 += e;
	}
	float factor = 0;
	//float denom = (max2 - min2);
	//if (denom != 0) { factor = max1 / denom; };
	float denom = (sum2 -min2*numz*lengthmz);
	if (denom != 0) { factor = sum1 / denom; };
	//if (sum2 != 0) { factor = sum1 / sum2; };
	if (factor > 0) {
		#pragma omp parallel for schedule(auto)
		for (int i = 0; i < lengthmz*numz; i++)
		{
			blur[i] -= min2;
			blur[i] *= factor;
		}
	}
	free(newblur);
	return;
}

void softargmax(float *blur, const int lengthmz, const int numz, const float beta)
{
	if (beta < 0) {
		//softargmax_transposed(blur, lengthmz, numz, fabs(beta));
		softargmax_everything(blur, lengthmz, numz, fabs(beta));
		return;
	}
	
	float *newblur;
	newblur = calloc(lengthmz*numz, sizeof(float));
	memcpy(newblur, blur, lengthmz*numz * sizeof(float));
	#pragma omp parallel for schedule(auto)
	for (int i = 0; i < lengthmz; i++)
	{
		float sum2 = 0;
		float sum1 = 0;
		float factor = 0;
		//float max1 = 0;
		//float max2 = 0;
		float min2 = 1000000000000;
		//float min1 = 1000000000000;
		for (int j = 0; j < numz; j++)
		{
			float d = newblur[index2D(numz, i, j)];
			sum1 += d;
			//if (d < min1) { min1 = d; }
			//if (d > max1) { max1 = d; }
			float e = exp(beta * d);
			//if (beta > 0) { e = exp(beta * d); }
			//else { e = pow(d, -beta); }
			if (e < min2) { min2 = e; }
			//if (e > max2) { max2 = e; }
			blur[index2D(numz, i, j)] = e;
			sum2 += e;
		}
		//float min = min2 - min1;
		//if (beta < 0) { min = min2; } 
		float denom = (sum2 - min2*numz);
		if (denom != 0) { factor = sum1 / denom; };
		//float factor = max1 / (max2 - min);
		//float sum3 = 0;
		if (factor > 0) {
			for (int j = 0; j < numz; j++)
			{
				blur[index2D(numz, i, j)] -= min2;
				blur[index2D(numz, i, j)] *= factor;
				//sum3 += blur[index2D(numz, i, j)];
			}
		}
		else {
			for (int j = 0; j < numz; j++)
			{
				blur[index2D(numz, i, j)] = 0;
			}
		}
		//printf("%f %f\n", sum1, sum3);
	}
	free(newblur);
	return;
}


/*

void softargmax_transposed(float *blur, const int lengthmz, const int numz, const float beta)
{
	float *newblur;
	newblur = calloc(lengthmz*numz, sizeof(float));
	memcpy(newblur, blur, lengthmz*numz * sizeof(float));
#pragma omp parallel for schedule(auto)
	for (int j = 0; j < numz; j++)
	{
		//float sum2 = 0;
		//float sum1 = 0;
		float factor = 0;
		float max1 = 0;
		float max2 = 0;
		float min2 = 1000000000000;
		//float min1 = 1000000000000;
		for (int i = 0; i < lengthmz; i++)
		{
			float d = newblur[index2D(numz, i, j)];
			//sum1 += d;
			//if (d < min1) { min1 = d; }
			if (d > max1) { max1 = d; }
			float e = exp(beta * d);
			//if (beta > 0) { e = exp(beta * d); }
			//else { e = pow(d, -beta); }
			if (e < min2) { min2 = e; }
			if (e > max2) { max2 = e; }
			blur[index2D(numz, i, j)] = e;
			//sum2 += e;
		}
		//float min = min2 - min1;
		//if (beta < 0) { min = min2; }
		//float denom = (sum2 - min2*numz);
		//if (denom != 0) { factor = sum1 / denom; };
		float denom = (max2 - min2);
		if (denom != 0) { factor = max1 / denom; };
		//float sum3 = 0;
		if (factor > 0) {
			for (int i = 0; i < lengthmz; i++)
			{
				blur[index2D(numz, i, j)] -= min2;
				blur[index2D(numz, i, j)] *= factor;
				//sum3 += blur[index2D(numz, i, j)];
			}
		}
		else {
			for (int i = 0; i < lengthmz; i++)
			{
				blur[index2D(numz, i, j)] = 0;
			}
		}
		//printf("%f %f\n", sum1, sum3);
	}
	free(newblur);
	return;
}

void point_smoothing_iso(float *blur, const int lengthmz, const int numz, const int width)
{
	float *newblur;
	newblur = calloc(lengthmz*numz, sizeof(float));
	memcpy(newblur, blur, lengthmz*numz * sizeof(float));
	#pragma omp parallel for schedule(auto)
	for (int i = 0; i < lengthmz; i++)
	{
		for (int j = 0; j < numz; j++)
		{
			int low = i - width;
			if (low < 0) { low = 0; }
			int high = i + width + 1;
			if (high > lengthmz) { high = lengthmz; }

			float sum = 0;
			float val = 0;
			for (int k = low; k < high; k++)
			{
				float e= exp(500*newblur[index2D(numz, k, j)]);
				sum += e;
				if (k == i) { val = e; }
			}
			if(sum!=0){ blur[index2D(numz, i, j)] = val / sum; }
			else { blur[index2D(numz, i, j)] = 0; }
		}
	}
	free(newblur);
	return;
}*/

void SetLimits(const Config config, Input *inp)
{
	//Determines the indexes of each test mass from mfile in m/z space
	int* testmasspos = malloc(sizeof(float) * config.mfilelen * config.numz);
	if (config.mflag == 1 && config.limitflag == 1) {
		for (int i = 0; i < config.mfilelen; i++)
		{
			for (int j = 0; j < config.numz; j++)
			{
				float mztest = (inp->testmasses[i] + config.adductmass * inp->nztab[j]) / (float)inp->nztab[j];
				testmasspos[index2D(config.numz, i, j)] = nearfast(inp->dataMZ, mztest, config.lengthmz);
			}
		}
	}
	
	//If there is a mass file read, it will only allow masses close to those masses within some config.mtabsig window.
	if (config.mflag == 1 && config.limitflag == 0)
	{
		TestMassListWindowed(config.lengthmz, config.numz, inp->barr, inp->mtab, config.nativezub, config.nativezlb, config.massub, config.masslb, inp->nztab, inp->testmasses, config.mfilelen, config.mtabsig);
	}
	//If there is a mass file read and the mass table window (config.mtabsig) is 0, it will only write intensities at the m/z values closest to the m/z values read in from the mfile.
	else if (config.mflag == 1 && config.limitflag == 1)
	{
		TestMassListLimit(config.lengthmz, config.numz, inp->barr, inp->mtab, config.nativezub, config.nativezlb, config.massub, config.masslb, inp->nztab, testmasspos, config.mfilelen);
	}
	//Normally, write the intensity values if the values fall within the mass upperbound and lower bound
	else
	{
		TestMass(config.lengthmz, config.numz, inp->barr, inp->mtab, config.nativezub, config.nativezlb, config.massub, config.masslb, inp->nztab);
	}
	free(testmasspos);
}

//Sets the maxlength parameter and the start and end values for the m/z peak shape convolution
//Convolution uses a reflection for the edges, so some care needs to be taken when things are over the edge.
int SetStartsEnds(const Config config, const Input * inp, int *starttab, int *endtab, const float threshold) {
	int maxlength = 1;
	for (int i = 0; i < config.lengthmz; i++)
	{
		float point = inp->dataMZ[i] - threshold;
		int start, end;
		if (point < inp->dataMZ[0] && config.speedyflag == 0) {
			//start = (int)((point - inp->dataMZ[0]) / (inp->dataMZ[1] - inp->dataMZ[0]));
			start = 0 - nearfast(inp->dataMZ, 2 * inp->dataMZ[0] - point, config.lengthmz);
		}
		else {
			start = nearfast(inp->dataMZ, point, config.lengthmz);
		}
		starttab[i] = start;

		point = inp->dataMZ[i] + threshold;
		if (point > inp->dataMZ[config.lengthmz - 1] && config.speedyflag == 0) {
			//end = config.lengthmz - 1 + (int)((point - inp->dataMZ[config.lengthmz - 1]) / (inp->dataMZ[config.lengthmz - 1] - inp->dataMZ[config.lengthmz - 2]));
			end = config.lengthmz - 1 + nearfast(inp->dataMZ, 2 * inp->dataMZ[0] - point, config.lengthmz);
		}
		else {
			end = nearfast(inp->dataMZ, point, config.lengthmz);
		}
		endtab[i] = end;
		if (end - start > maxlength) { maxlength = end - start; }
		//printf("%d %d\n", start, end);
	}
	//printf("Max Length: %d\t", maxlength);
	return maxlength;
}

void WritePeaks(const Config config, const Decon* decon) {
	char outdat[1024];
	strjoin(config.dataset, "/peaks", outdat);
	float* ptemp = NULL;
	ptemp = calloc(decon->plen * 3, sizeof(float));

	for (int i = 0; i < decon->plen; i++) {
		ptemp[i * 3] = decon->peakx[i];
		ptemp[i * 3 + 1] = decon->peaky[i];
		ptemp[i * 3 + 2] = decon->dscores[i];
	}

	mh5writefile2d_grid(config.file_id, outdat, decon->plen, 3, ptemp);
	free(ptemp);
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

void WriteDecon(const Config config, const Decon * decon, const Input * inp)
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

		if (config.filetype == 0) {
			char outstring9[500];
			sprintf(outstring9, "%s_grid.bin", config.outfile);
			out_ptr = fopen(outstring9, "wb");
			if (config.rawflag == 0) { fwrite(decon->newblur, sizeof(float), config.lengthmz * config.numz, out_ptr); }
			if (config.rawflag == 1) { fwrite(decon->blur, sizeof(float), config.lengthmz * config.numz, out_ptr); }
			fclose(out_ptr);
			//printf("m/z grid: %s\t", outstring9);
		}
		else {
			strjoin(config.dataset, "/mz_grid", outdat);
			if (config.rawflag == 0) { mh5writefile1d(file_id, outdat, config.lengthmz * config.numz, decon->newblur); }
			if (config.rawflag == 1) { mh5writefile1d(file_id, outdat, config.lengthmz * config.numz, decon->blur); }
			
			float* chargedat = NULL;
			chargedat = calloc(config.numz, sizeof(float));
			float* chargeaxis = NULL;
			chargeaxis = calloc(config.numz, sizeof(float));

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
		if (config.filetype == 0) {
			char outstring10[500];
			sprintf(outstring10, "%s_massgrid.bin", config.outfile);
			out_ptr = fopen(outstring10, "wb");
			fwrite(decon->massgrid, sizeof(float), decon->mlen * config.numz, out_ptr);
			fclose(out_ptr);
			//printf("Mass Grid: %s\t", outstring10);
		}
		else {
			strjoin(config.dataset, "/mass_grid", outdat);
			mh5writefile1d(file_id, outdat, decon->mlen * config.numz, decon->massgrid);
		}
	}

	//Writes the mass values convolved with the peak shape
	if (config.rawflag == 0 || config.rawflag == 1 || config.rawflag == 2 || config.rawflag == 3) {
		if (config.filetype == 0) {
			char outstring4[500];
			sprintf(outstring4, "%s_mass.txt", config.outfile);
			out_ptr = fopen(outstring4, "w");
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

	if (config.filetype == 1 && decon->plen>0) {
		WritePeaks(config, decon);
	}

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
			printf("HDF5 Data Set: %s\n", config->dataset);
		}
		else
		{
			strcpy(config->dataset, "/ms_data");
		}

		char outdat[1024];
		config->file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
		strjoin(config->dataset, "/processed_data", outdat);
		config->lengthmz = mh5getfilelength(config->file_id, outdat);
		inp->dataMZ = calloc(config->lengthmz, sizeof(float));
		inp->dataInt = calloc(config->lengthmz, sizeof(float));
		mh5readfile2d(config->file_id, outdat, config->lengthmz, inp->dataMZ, inp->dataInt);
		printf("Length of Data: %d \n", config->lengthmz);

		//Check the length of the mfile and then read it in.
		if (config->mflag == 1)
		{
			config->mfilelen = mh5getfilelength(config->file_id, "/config/masslist");
			printf("Length of mfile: %d \n", config->mfilelen);
			inp->testmasses = malloc(sizeof(float) * config->mfilelen);
			mh5readfile1d(config->file_id, "/config/masslist", inp->testmasses);
		}
		else {
			inp->testmasses = malloc(sizeof(float) * config->mfilelen);
		}
	}
	else {

		//Calculate the length of the data file automatically
		config->lengthmz = getfilelength(config->infile);
		inp->dataMZ = calloc(config->lengthmz, sizeof(float));
		inp->dataInt = calloc(config->lengthmz, sizeof(float));

		readfile(config->infile, config->lengthmz, inp->dataMZ, inp->dataInt);//load up the data array
		printf("Length of Data: %d \n", config->lengthmz);

		//Check the length of the mfile and then read it in.

		if (config->mflag == 1)
		{
			config->mfilelen = getfilelength(config->mfile);
			printf("Length of mfile: %d \n", config->mfilelen);
		}
		inp->testmasses = malloc(sizeof(float) * config->mfilelen);
		if (config->mflag == 1)
		{
			readmfile(config->mfile, config->mfilelen, inp->testmasses);//read in mass tab
		}
	}


	//This for loop creates a list of charge values
	inp->nztab = calloc(config->numz, sizeof(int));
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
		if (inp->dataMZ[i] == inp->dataMZ[i + 1]) { printf("Error: Two data points are identical"); exit(104); }
	}

}

#endif