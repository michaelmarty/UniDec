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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "hdf5.h"
#include "hdf5_hl.h"

void mh5readfile3d(hid_t file_id, char *dataname, int lengthmz, double *dataMZ, double *dataInt, double *data3);
int mh5getfilelength(hid_t file_id, char *dataname);

typedef struct Config Config;

struct Config
{
	char infile[500];
	char outfile[500];
	int numit;
	int numz;
	int endz;
	int startz;
	double zsig;
	double mzsig;
	double msig;
	double molig;
	double massub;
	double masslb;
	int psfun;
	double mtabsig;
	char mfile[500];
	char manualfile[500];
	int mflag;
	double massbins;
	int limitflag;
	double cutoff;
	double psthresh;
	int speedyflag;
	int linflag;
	int aggressiveflag;
	double adductmass;
	int rawflag;
	double nativezub;
	double nativezlb;
	int poolflag;
	int manualflag;
	double intthresh;
	double peakshapeinflate;
	double killmass;
	int fixedmassaxis;
	int isotopemode;
	int filetype;
	int imflag;
	//IM Parameters
	double dtsig;
	double csig;
	double ccsub;
	double ccslb;
	double ccsbins;
	double temp;
	double press;
	double volt;
	double tcal1;
	double tcal2;
	double twaveflag;
	double hmass;
	double to;
	double len;
	double edc;
	double nativeccsub;
	double nativeccslb;
	int baselineflag;
	int noiseflag;
	int zout;
	int metamode;
	double minmz;
	double maxmz;
	int mzbins;
	double bsub;
	double peakwin;
	double peakthresh;
	double exwindow;
	int exchoice;
	int exnorm;
	int peaknorm;
	int orbimode;
	int datanorm;
	//Experimental Parameters
	int filterwidth;
	int zerolog;
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
	config.cutoff = 0;
	config.psthresh = 6;
	config.speedyflag = 1;
	config.aggressiveflag = 0;
	config.adductmass = 0;
	config.rawflag = 0;
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
	config.peakwin = 500;
	config.peakthresh = 0.1;
	config.exchoice = 0;
	config.exnorm = 0;
	config.peaknorm = 1;
	config.exwindow = 0;
	config.orbimode = 0;
	config.datanorm = 1;
	//Experimental
	config.filterwidth = 20;
	config.zerolog = -8;
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
		printf("outfile = %s\n", config.outfile);
		printf("\n");
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
		printf("\nRead from file:");
		while (fscanf(file, "%s %500[^\n]", x, y) != EOF)
		{
			//printf( "read in: %s %s \n", x,y );
			if (strstr(x, "input") != NULL){ strcpy(config.infile, y); printf(" input"); }
			if (strstr(x, "output") != NULL){ strcpy(config.outfile, y); printf(" output"); }
			if (strstr(x, "mfile") != NULL){ strcpy(config.mfile, y); config.mflag = 1;  printf(" mfile"); }
			if (strstr(x, "numit") != NULL){ config.numit = atoi(y); printf(" numit"); }
			//if (strstr(x, "numz") != NULL){ config.numz = atoi(y); printf(" numz"); }
			if (strstr(x, "startz") != NULL){ config.startz = atoi(y); printf(" startz"); }
			if (strstr(x, "endz") != NULL) { config.endz = atoi(y); printf(" endz"); }
			if (strstr(x, "zzsig") != NULL){ config.zsig = atof(y); printf(" zzsig"); }
			if (strstr(x, "mzsig") != NULL){ config.mzsig = atof(y); printf(" mzsig"); }
			if (strstr(x, "msig") != NULL){ config.msig = atof(y); printf(" msig"); }
			if (strstr(x, "molig") != NULL){ config.molig = atof(y); printf(" molig"); }
			if (strstr(x, "massub") != NULL){ config.massub = atof(y); printf(" massub"); }
			if (strstr(x, "masslb") != NULL){ config.masslb = atof(y); printf(" masslb"); }
			if (strstr(x, "psfun") != NULL){ config.psfun = atoi(y); printf(" psfun"); }
			if (strstr(x, "mtabsig") != NULL){ config.mtabsig = atof(y); printf(" mtabsig"); }
			if (strstr(x, "massbins") != NULL){ config.massbins = atof(y); printf(" massbins"); }
			if (strstr(x, "psthresh") != NULL){ config.psthresh = atof(y); printf(" psthresh"); }
			if (strstr(x, "speedy") != NULL){ config.speedyflag = atoi(y); printf(" speedy"); }
			if (strstr(x, "aggressive") != NULL){ config.aggressiveflag = atoi(y); printf(" aggressive"); }
			if (strstr(x, "adductmass") != NULL){ config.adductmass = atof(y); printf(" adductmass"); }
			if (strstr(x, "rawflag") != NULL){ config.rawflag = atoi(y); printf(" rawflag"); }
			if (strstr(x, "nativezub") != NULL){ config.nativezub = atof(y); printf(" nativezub"); }
			if (strstr(x, "nativezlb") != NULL){ config.nativezlb = atof(y); printf(" nativezlb"); }
			if (strstr(x, "poolflag") != NULL){ config.poolflag = atoi(y); printf(" poolflag"); }
			if (strstr(x, "manualfile") != NULL){ config.manualflag = 1; strcpy(config.manualfile, y); printf(" manualfile"); }
			if (strstr(x, "intthresh") != NULL){ config.intthresh = atof(y); printf(" intthresh"); }
			if (strstr(x, "peakshapeinflate") != NULL){ config.peakshapeinflate = atof(y); printf(" peakshapeinflate"); }
			if (strstr(x, "killmass") != NULL){ config.killmass = atof(y); printf(" killmass"); }
			if (strstr(x, "isotopemode") != NULL){ config.isotopemode = atoi(y); printf(" isotopemode"); }
			if (strstr(x, "orbimode") != NULL) { config.orbimode = atoi(y); printf(" orbimode"); }
			if (strstr(x, "imflag") != NULL) { config.imflag = atoi(y); printf(" imflag"); }
			if (strstr(x, "linflag") != NULL) { config.linflag = atoi(y); printf(" linflag"); }
			//IM Parameters
			if (strstr(x, "csig") != NULL) { config.csig = atof(y); printf(" csig"); }
			if (strstr(x, "dtsig") != NULL) { config.dtsig = atof(y); printf(" dtsig"); }
			if (strstr(x, "ccsub") != NULL) { config.ccsub = atof(y); printf(" ccsub"); }
			if (strstr(x, "ccslb") != NULL) { config.ccslb = atof(y); printf(" ccslb"); }
			if (strstr(x, "ccsbins") != NULL) { config.ccsbins = atof(y); printf(" ccsbins"); }
			if (strstr(x, "temp") != NULL) { config.temp = atof(y); printf(" temp"); }
			if (strstr(x, "pressure") != NULL) { config.press = atof(y); printf(" pressure"); }
			if (strstr(x, "volt") != NULL) { config.volt = atof(y); printf(" volt"); }
			if (strstr(x, "gasmass") != NULL) { config.hmass = atof(y); printf(" gasmass"); }
			if (strstr(x, "tnaught") != NULL) { config.to = atof(y); printf(" to"); }
			if (strstr(x, "tcal1") != NULL) { config.tcal1 = atof(y); printf(" tcal1"); }
			if (strstr(x, "tcal2") != NULL) { config.tcal2 = atof(y); printf(" tcal2"); }
			if (strstr(x, "edc") != NULL) { config.edc = atof(y); printf(" edc"); }
			if (strstr(x, "zout") != NULL) { config.zout = atoi(y); printf(" zout"); }
			if (strstr(x, "twaveflag") != NULL) { config.twaveflag = atoi(y); printf(" twaveflag"); }
			if (strstr(x, "ubnativeccs") != NULL) { config.nativeccsub = atof(y); printf(" ubnativeccs"); }
			if (strstr(x, "lbnativeccs") != NULL) { config.nativeccslb = atof(y); printf(" lbnativeccs"); }
			if (strstr(x, "driftlength") != NULL) { config.len = atof(y); printf(" driftlength"); }
			if (strstr(x, "baselineflag") != NULL) { config.baselineflag = atoi(y); printf(" baselineflag"); }
			if (strstr(x, "noiseflag") != NULL) { config.noiseflag = atoi(y); printf(" noiseflag"); }
			//Experimental
			if (strstr(x, "filterwidth") != NULL) { config.filterwidth = atoi(y); printf(" filterwidth"); }
			if (strstr(x, "zerolog") != NULL) { config.zerolog = atoi(y); printf(" zerolog"); }
		}
		printf("\n\n");
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
	printf("\t\"isotopemode\" \t0=off 1=on: Uses isotope distributions in deconvolution (MS Only)\n");
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
	printf("\nEnjoy! Please report bugs to Michael Marty (mtmarty@email.arizona.edu) v.622.\n");
	//printf("\nsize of: %d",sizeof(char));
}

double ndis(double x, double y, double sig)
{
	return 1 / (sig*sqrt(2 * 3.14159265359))*exp(-(pow(x - y, 2)) / (2 * sig*sig));
}

//Actual Modulus operator rather than Remainder operator %
int mod(int a, int b) { int r = a % b; return r < 0 ? r + b : r; }

//Reads in x y file.
void readfile(char *infile,int lengthmz,double *dataMZ,double *dataInt)
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
void readfile3(char *infile,int lengthmz,double *array1,double *array2,double *array3)
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
void readmfile(char *infile,int mfilelen,double *testmasses)
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
		printf("TESTCould not open %s file\n", infile);
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
int nearunsorted(double *testmasses, double point, int lengthtest)
{
	double minval = fabs(point - testmasses[0]);
	double val = testmasses[0];
	double difftest;
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
int neartest(double *testmasses,double point,int lengthtest,double cutoff)
{
    double minval=fabs(point-testmasses[0]);
    double val=testmasses[0];
    double difftest;
	for (int i = 0; i<lengthtest; i++)
    {
        difftest=fabs(point-testmasses[i]);
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
int nearfast(const double *dataMZ,const double point,const int numdat)
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

//Perform a linear interpolation. I've left the code for other interpolation functions below but they don't seem to matter.
double LinearInterpolate(double y1,double y2,double mu)
{
   return(y1*(1-mu)+y2*mu);
}

double LinearInterpolatePosition(double x1, double x2, double x)
{
	return (x - x1) / (x2 - x1);
}

double CubicInterpolate(double y0,double y1,double y2,double y3,double mu)
{double a0,a1,a2,a3,mu2;
   mu2 = mu*mu;
   a0 = y3 - y2 - y0 + y1;
   a1 = y0 - y1 - a0;
   a2 = y2 - y0;
   a3 = y1;
   return(a0*mu*mu2+a1*mu2+a2*mu+a3);}

double CRSplineInterpolate(double y0,double y1,double y2,double y3,double mu)
{double a0,a1,a2,a3,mu2;
   mu2 = mu*mu;
   a0 = -0.5*y0 + 1.5*y1 - 1.5*y2 + 0.5*y3;
   a1 = y0 - 2.5*y1 + 2*y2 - 0.5*y3;
   a2 = -0.5*y0 + 0.5*y2;
   a3 = y1;
   return(a0*mu*mu2+a1*mu2+a2*mu+a3);}

double clip(double x,double cutoff)
{
	if(x>cutoff){return x;}
	else{return 0;}
}

//Function for defining m/z peak shape. Other peak shape functions could be easily added her.
double mzpeakshape(double x,double y,double sig,int psfun)
{
    double result;
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

inline int index2D(int ncols, int r, int c)
{ return r*ncols + c; }

inline int indexmod(int length, int r, int c)
{
	return mod((c - r), length);
}

inline int index3D(int ncols,int nrows, int r, int c, int d)
{
	return r*ncols*nrows + c*nrows+d;
}


double Max(const double *blur, const int length){
	double blurmax = 0;
	for (int i = 0; i<length; i++)
	{
		if (blur[i]>blurmax)
		{
			blurmax = blur[i];
		}
	}
	return blurmax;
}

double Sum(const double *blur, int length) {
	double sum = 0;
	for (int i = 0; i<length; i++)
	{
		sum += blur[i];
	}
	return sum;
}

double max1d(double *blur, int lengthmz) {
	double blurmax = blur[0];
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

double min1d(double *blur, int lengthmz) {
	double blurmin = blur[0];
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

//Convolution of neighborhood function with gaussian filter.
void blur_it(const int lengthmz,
                    const int numz,
                    const int numclose,
                    const int * __restrict closeind,
                    const double * __restrict closeval,
                    double * __restrict newblur,
                    const double * __restrict blur,
					const char * __restrict barr)
{
  int i, j, k;
  if (numclose == 1)
  {
	  memcpy(newblur, blur, lengthmz*numz * sizeof(double));
  }
  else {
	#pragma omp parallel for private (i,k,j), schedule(dynamic)
	  for (i = 0; i < lengthmz; i++)
	  {
		  for (j = 0; j < numz; j++)
		  {
			  double temp = 0;
			  if (barr[index2D(numz, i, j)] == 1)
			  {
				  for (k = 0; k < numclose; k++)
				  {
					  if (closeind[index3D(numz, numclose, i, j, k)] != -1)
					  {
						  temp += closeval[k] * blur[closeind[index3D(numz, numclose, i, j, k)]];
					  }
				  }
			  }
			  newblur[index2D(numz, i, j)] = temp;
		  }
	  }
  }
}

//Charge state smooth using a mean filter of the log
void blur_it_mean(const int lengthmz,
                    const int numz,
                    const int numclose,
                    const int * __restrict closeind,
                    double * __restrict newblur,
                    const double * __restrict blur,
					const char * __restrict barr,
					const int zerolog)
{
  unsigned  int i, j, k;
  double zeroval = zerolog;
  if (numclose == 1)
  {
	  memcpy(newblur, blur, lengthmz*numz*sizeof(double));
  }
  else{
	#pragma omp parallel for private (i,j, k), schedule(dynamic)
	  for (i = 0; i < lengthmz; i++)
	  {
		  for (j = 0; j < numz; j++)
		  {
			  double temp = 0;
			  if (barr[index2D(numz, i, j)] == 1)
			  {
				  for (k = 0; k < numclose; k++)
				  {
					  if (closeind[index3D(numz, numclose, i, j, k)] != -1)
					  {
						  double temp2 = blur[closeind[index3D(numz, numclose, i, j, k)]];
						  if (temp2 > 0){ temp += log(temp2); }
						  else{ temp += zeroval; }
					  }
					  else{ temp += zeroval; }
				  }
				  temp = exp(temp / (double)numclose);
			  }
			  newblur[index2D(numz, i, j)] = temp;
		  }
	  }
  }
}


//Second Derivative of a Gaussian
double secderndis(double m, double s,double x){
	return s/2*(2*exp(-pow(m-x,2)/s))/s-(4*exp(-pow(m-x,2)/s)*pow(m-x,2))/pow(s,2);
}


void blur_baseline(double *baseline, const int lengthmz, const double *dataMZ, const double mzsig, int mult, int filterwidth)
{
	int mulin = mult;
	double *temp = NULL;
	temp = calloc(lengthmz, sizeof(double));
	memcpy(temp, baseline, sizeof(double)*lengthmz);
	int i,j;
	#pragma omp parallel for private (i,j), schedule(dynamic)
	for (i = 0; i < lengthmz; i++)
	{
		double mzdiff = 0;
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

		double val = 0;
		int window = filterwidth;
		for (j = -window; j < window; j++)
		{
			int k = i + j*(mult);
			double newval = 0;
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
	double *x = (double *)a;
	double *y = (double *)b;
	if (*x < *y) { return -1; }
	else if (*x > *y) { return 1; }
	return 0;
}

void midblur_baseline(double *baseline, const int lengthmz,const double *dataMZ, const double mzsig, int mult)
{
	if (mult == 0)
	{
		mult = lengthmz / 400;
	}
	int window = 25;

	double *temp = NULL;
	temp = calloc(lengthmz, sizeof(double));
	memcpy(temp, baseline, sizeof(double)*lengthmz);
	int i, j;
	#pragma omp parallel for private(i,j), schedule(dynamic)
	for (i = 0; i < lengthmz; i++)
	{
		double val = 0;
		double len = window * 2;
		double *med = NULL;
		med = calloc(len, sizeof(double));
		int index = 0;
		for (j = -window; j < window; j++)
		{
			int k = i + j * mult;
			double newval = 0;
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
		qsort(med, len, sizeof(double), compare_function);
		index = 0;
		for (j = 0; j < window; j++)
		{
			val += med[j];
			index++;
		}
		baseline[i] = val/((float)index);
		free(med);
	}
	free(temp);
}

void blur_noise(double *noise, int lengthmz)
{
	double *temp = NULL;
	temp = calloc(lengthmz, sizeof(double));
	memcpy(temp, noise, sizeof(double)*lengthmz);
	int i, j;
	double filter[5] = { -0.1,-0.4,1,-0.4,-0.1 };
	#pragma omp parallel for private (i,j), schedule(dynamic)
	for (i = 0; i < lengthmz; i++)
	{
		double val = 0;
		for (j = -2; j <=2 ; j++)
		{
			int k = i + j;
			double newval=0;
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
}

inline int fixk(int k, int lengthmz)
{
	k = abs(k);
	if (k >= lengthmz){ k = 2 * lengthmz - k-2; }
	return k;
}

void convolve_simp(int lengthmz, int maxlength, int*starttab, int *endtab, double *mzdist, double *deltas, double *denom, int speedyflag)
{

	unsigned int i, k;
	if (speedyflag == 0){
		#pragma omp parallel for private (i,k), schedule(dynamic)
		for (i = 0; i < lengthmz; i++)
		{
			double cv = 0;
			for (k = starttab[i]; k <= endtab[i]; k++)
			{
				int k2 = fixk(k, lengthmz);
				int start = starttab[k2];
					cv += deltas[k2] * mzdist[index2D(maxlength, k2, i - start)];
			}
			denom[i] = cv;
		}
	}
	else{
		#pragma omp parallel for private (i,k), schedule(dynamic)
		for (i = 0; i < lengthmz; i++)
		{
			double cv = 0;
			for (k = starttab[i]; k <= endtab[i]; k++)
			{
				cv += deltas[k] * mzdist[indexmod(lengthmz,k, i)];
			}
			denom[i] = cv;
		}
	}
}


void sum_deltas(const int lengthmz, const int numz, const double * __restrict blur, const char * __restrict barr,
	const int isolength, const int *__restrict isotopepos, const float *__restrict isotopeval, double *deltas)
{
	int i,j,k;
	//Collapse the grid into a 1D array of delta function values
	if (isolength == 0) {
	#pragma omp parallel for private (i,j), schedule(dynamic)
		for (i = 0; i < lengthmz; i++)
		{
			double temp = 0;
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
		//#pragma omp parallel for private (i,j,k), schedule(dynamic)
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				if (barr[index2D(numz, i, j)] == 1) {
					double topval = blur[index2D(numz, i, j)];
					for (k = 0; k < isolength; k++)
					{
						int pos = isotopepos[index3D(numz, isolength, i, j, k)];
						float val = isotopeval[index3D(numz, isolength, i, j, k)];
						deltas[pos] += topval* (double)val;
					}
				}
			}
		}
	}
}

void apply_ratios(const int lengthmz, const int numz, const double * __restrict blur, const char * __restrict barr,
	const int isolength, const int *__restrict isotopepos, const float *__restrict isotopeval, const double * __restrict denom, double *blur2)
{
	int i, j, k;
	if (isolength == 0)
	{
		#pragma omp parallel for private (i,j), schedule(dynamic)
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
		#pragma omp parallel for private (i,j,k), schedule(dynamic)
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				if (barr[index2D(numz, i, j)] == 1)
				{
					double temp = 0;
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

void deconvolve_baseline(const int lengthmz, const double *dataMZ, const double * dataInt, double * baseline, const double mzsig)
{
	double *denom = NULL;
	denom = calloc(lengthmz, sizeof(double));

	midblur_baseline(baseline, lengthmz, dataMZ,mzsig, 0);
	midblur_baseline(baseline, lengthmz, dataMZ, mzsig, 5);

	memcpy(denom, baseline, sizeof(double)*lengthmz);
	int i;
	#pragma omp parallel for private(i), schedule(dynamic)
	for (i = 0; i < lengthmz; i++)
	{
		if (denom[i] != 0 && dataInt[i] >= 0) { denom[i] = dataInt[i] / denom[i]; }
	}

	midblur_baseline(denom, lengthmz, dataMZ, mzsig, 0);
	midblur_baseline(denom, lengthmz, dataMZ, mzsig, 5);
	#pragma omp parallel for private(i), schedule(dynamic)
	for (i = 0; i < lengthmz; i++)
	{
		baseline[i] = baseline[i] * (denom[i]);
	}
	free(denom);
}

double deconvolve_iteration_speedy(const int lengthmz, const int numz,const int maxlength, const double * __restrict blur, double * __restrict blur2,
						  const char * __restrict barr, const int aggressiveflag,const double * __restrict  dataInt,
							const int isolength,const int *__restrict isotopepos,const float *__restrict isotopeval
						  ,int*starttab, int *endtab, double *mzdist, int speedyflag,int baselineflag, double *baseline,
							double *noise, const double mzsig, const double *dataMZ, const double filterwidth)
{
	unsigned int i,j,k;
	double *deltas=NULL,*denom=NULL;
	deltas=calloc(lengthmz,sizeof(double));
	denom=calloc(lengthmz,sizeof(double));

	if (aggressiveflag==1) {
		blur_baseline(baseline, lengthmz, dataMZ, mzsig, 0, filterwidth);
		//blur_baseline(baseline, lengthmz, 10);
		//blur_noise(noise, lengthmz); 
	}

	//Sum deltas
	sum_deltas(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, deltas);

	//Convolve with peak shape
	convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, denom, speedyflag);	
	
	if (aggressiveflag == 1)
	{
		#pragma omp parallel for private(i), schedule(dynamic)
		for(i=0;i<lengthmz;i++){
			denom[i] += baseline[i];// +noise[i]);
		}
	}

	//Calculate Ratio
	#pragma omp parallel for private(i), schedule(dynamic)
	for (i = 0; i < lengthmz; i++)
	{
		if (denom[i] != 0 && dataInt[i] >= 0) { denom[i] = dataInt[i] / denom[i]; }
	}
	
	//Multiply Ratio by prior
	apply_ratios(lengthmz, numz, blur, barr, isolength, isotopepos, isotopeval, denom, blur2);


	if (aggressiveflag ==1)
	{
		//memcpy(deltas, denom, sizeof(double)*lengthmz);
		blur_baseline(denom, lengthmz, dataMZ, mzsig, 0, filterwidth);
		//blur_baseline(denom, lengthmz, 10);
		//blur_noise(deltas, lengthmz);
		#pragma omp parallel for private(i), schedule(dynamic)
		for (i = 0; i < lengthmz; i++)
		{
			baseline[i] = baseline[i] * (denom[i]);
			//noise[i] = noise[i]*(deltas[i]);
		}
	}
	
 free(deltas);
 free(denom);
 return 0;
}

void ApplyCutoff1D(double *array, double cutoff, int lengthmz)
{
	unsigned int i;
	//#pragma omp parallel for private (i,j), schedule(dynamic)
	for (i = 0; i<lengthmz; i++)
	{
		if (array[i]<cutoff) {array[i] = 0;}
	}
}

double getfitdatspeedy(double *fitdat, double *blur, int lengthmz,int numz,int maxlength,double maxint,
	int isolength, int *isotopepos, float *isotopeval, int*starttab, int *endtab, double *mzdist, int speedyflag)
{
	unsigned int i,j,k;
	double *deltas=NULL;
	deltas=calloc(lengthmz,sizeof(double));
	if (isolength == 0){
		#pragma omp parallel for private (i,j), schedule(dynamic)
		for (i = 0; i < lengthmz; i++) //Collapse the grid into a 1D array of delta function values
		{
			double temp = 0;
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
				double topval = blur[index2D(numz, i, j)];
				for (k = 0; k < isolength; k++)
				{
					int pos = isotopepos[index3D(numz, isolength, i, j, k)];
					float val = isotopeval[index3D(numz, isolength, i, j, k)];
					deltas[pos] += topval* (double)val;
				}
			}
		}

	}
	convolve_simp(lengthmz, maxlength, starttab, endtab, mzdist, deltas, fitdat, speedyflag);
	free(deltas);
	double fitmax=0;
	//#pragma omp parallel for private (i), schedule(dynamic)
	for(i=0;i<lengthmz;i++)
	{
		if(fitdat[i]>fitmax){fitmax=fitdat[i];}
	}
	//#pragma omp parallel for private (i), schedule(dynamic)
	for(i=0;i<lengthmz;i++)
	{
    	if(fitdat[i]<0){fitdat[i]=0;}
    	else{fitdat[i]=fitdat[i]*maxint/fitmax;}
	}
	return fitmax;
}

double errfunspeedy(double *dataInt,double *fitdat, double *blur,int lengthmz,int numz,int maxlength,double maxint,
	int isolength, int *isotopepos, float *isotopeval, int*starttab, int *endtab, double *mzdist, int speedyflag)
{
    int i;
	getfitdatspeedy(fitdat, blur,lengthmz, numz,maxlength, maxint,isolength, isotopepos, isotopeval,starttab,endtab,mzdist,speedyflag);
    double error=0;
    for(i=0;i<lengthmz;i++)
    {
        error+=pow((fitdat[i]-dataInt[i]),2);
    }
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
double nativecharge(double mass,double fudge)
{
	return 0.0467*pow(mass,0.533)+fudge;
}

//Calculate Average
double Average(int length,double *xarray)
{
	double temp1=0;
	double temp2=0;
	int i;
	for(i=0;i<length;i++)
	{
		temp1+=xarray[i];
		temp2+=1;
	}
	return temp1/temp2;
}

//Calculate Standard Deviation
double StdDev(int length,double *xarray,double wmean)
{
	double temp1=0;
	double temp2=0;
	int i;
	for(i=0;i<length;i++)
	{
		temp1+=pow(xarray[i]-wmean,2);
		temp2+=1;
	}
	return sqrt(temp1/temp2);
}

//Confidence score
double Confidence(double fit, double data)
{
	return clip(1-fabs(fit-data)/data,0);
}

void KillMass(double killmass, int lengthmz, int numz, char *barr, int *nztab, double adductmass, double *dataMZ, double psfun, double mzsig)
{
	double thresh;
	if (psfun == 0){ thresh = mzsig*2.35482; }
	else{ thresh = mzsig; }
	int j, i1, i2, k;
	double testmz;
	//#pragma omp parallel for private (j,testmz,i1,i2), schedule(dynamic)
	for (j = 0; j < numz; j++)
	{
		testmz = (killmass + adductmass*nztab[j]) / (double)nztab[j];
		i1 = nearfast(dataMZ, testmz - thresh, lengthmz);
		i2 = nearfast(dataMZ, testmz + thresh, lengthmz);
		for (k = i1; k <= i2; k++){ barr[index2D(numz, k, j)] = 0; }
	}
}

void TestMassListWindowed(int lengthmz, int numz, char *barr, double *mtab, double nativezub, double nativezlb, double massub, double masslb, int *nztab, double *testmasses, int mfilelen, double mtabsig)
{
	unsigned int i, j;
	double testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(dynamic)
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

void TestMassListLimit(int lengthmz, int numz, char *barr, double *mtab, double nativezub, double nativezlb, double massub, double masslb, int *nztab, int *testmasspos, int mfilelen)
{
	unsigned int i, j, k;
	double testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(dynamic)
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

void TestMass(int lengthmz, int numz, char *barr, double *mtab, double nativezub, double nativezlb, double massub, double masslb, int *nztab)
{
	unsigned int i, j;
	double testmass, nativelimit;
#pragma omp parallel for private (i,j,testmass,nativelimit), schedule(dynamic)
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

void ManualAssign(char *manualfile, int lengthmz, int numz, double *dataMZ, char *barr, int *nztab, hid_t file_id, Config config)
{
	unsigned int i, j,k;
	int manlen = 0;
	if (config.filetype == 1)
	{
		manlen = mh5getfilelength(file_id, "/config/manuallist");
	}
	else
	{
		manlen = getfilelength(manualfile);
	}
	printf("Length of Manual List: %d \n", manlen);
	double *manualmz = malloc(sizeof(double)*manlen);
	double *manualwin = malloc(sizeof(double)*manlen);
	double *manualassign = malloc(sizeof(double)*manlen);
	if (config.filetype == 1)
	{
		mh5readfile3d(file_id, "/config/manuallist", manlen, manualmz, manualwin, manualassign);
	}
	else
	{
		readfile3(manualfile, manlen, manualmz, manualwin, manualassign);
	}
	
	for( i = 0; i < manlen; i++){
		if (manualassign[i] < 0){ manualmz[i] = -manualmz[i]*1000.0; }
		//Cheating a bit...make the manualmz very negative if deassign, so that they aren't discovered by the assign loop
		}

	double testmz;
	int closest;
	#pragma omp parallel for private (i,j,testmz,closest), schedule(dynamic)
	for(i=0;i<lengthmz;i++)
	{
		for(j=0;j<numz;j++)
		{
			//Manual Assign: The nearest value wins in the case of overlap
			testmz=dataMZ[i];
			closest=nearunsorted(manualmz,testmz,manlen);
			if (fabs(manualmz[closest] - testmz)<manualwin[closest] && (double)nztab[j] == manualassign[closest] && manualassign[closest]>0)
				{
				barr[index2D(numz, i, j)]=1;
				}
			else if (fabs(manualmz[closest] - testmz)<manualwin[closest] && (double)nztab[j] != manualassign[closest] && manualassign[closest]>0)
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

void MakeBlur(int lengthmz,int numz,int numclose,char *barr,int *closezind,int *closemind,double *mtab,double molig,double adductmass,int *nztab,double *dataMZ,int *closeind)
{
	  unsigned int i,j,k;
	  #pragma omp parallel for private (i,j,k), schedule(dynamic)
	  for(i=0;i<lengthmz;i++)
	    {
	      for(j=0;j<numz;j++)
	        {
	            if(barr[index2D(numz,i,j)]==1)
	            {
	              for(k=0;k<numclose;k++)
	                {
	                  int indz=(int)(j+closezind[k]);
					  int ind;
					  if (indz < 0 || indz >= numz) { closeind[index3D(numz, numclose, i, j, k)] = -1; }
	                  else
	                  {
	                      double point=(double)((mtab[index2D(numz, i, j)]+closemind[k]*molig+adductmass*(nztab[j]+closezind[k]))/(double)(nztab[j]+closezind[k]));
	                      if(point<dataMZ[0]||point>dataMZ[lengthmz-1]){ closeind[index3D(numz, numclose, i, j, k)] = -1; }
						  else {
							  ind = nearfast(dataMZ, point, lengthmz);
							  closeind[index3D(numz, numclose, i, j, k)] = index2D(numz, ind, indz);
						  }
	                  }
	                }
	            }
	            else
	            {
	                for (k=0;k<numclose;k++)
	                {
						closeind[index3D(numz, numclose, i, j, k)] = -1;
	                }
	            }
	        }
	    }
}



void MakePeakShape2D(int lengthmz,int maxlength,int *starttab,int *endtab,double *dataMZ,double mzsig,int psfun,int speedyflag,double *mzdist)
{
	unsigned int i,j;
	#pragma omp parallel for private (i,j), schedule(dynamic)
	for(i=0;i<lengthmz;i++)
	  {
		  int start = starttab[i];
		  int end = endtab[i];
		  for (j = start; j <= end; j++)
		  {	
			  int j2 = fixk(j, lengthmz);
			  mzdist[index2D(maxlength, i, j2 - start)] = mzpeakshape(dataMZ[i], dataMZ[j2], mzsig, psfun);
		  }
	  }
}

void MakePeakShape1D(double *dataMZ,double threshold,int lengthmz,int speedyflag,double mzsig,int psfun,double *mzdist)
{
	double binsize=dataMZ[1]-dataMZ[0];
	double newrange=threshold/binsize;
	int n;
	for(n=(int)-newrange;n<(int)newrange;n++)
	{
		mzdist[indexmod(lengthmz,0,n)]=mzpeakshape(0,n*binsize,mzsig,psfun);
	}
	printf("\nNotice: Assuming linearized data. \n\n");
}


double Reconvolve(int lengthmz, int numz, int maxlength,int *starttab, int *endtab, double *mzdist,double *blur, double *newblur,int speedyflag)
{
	double newblurmax=0;
	unsigned int i,j,k;
	if (speedyflag == 0){
		#pragma omp parallel for private (i,j,k), schedule(dynamic)
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{

				double cv = 0;
				for (k = starttab[i]; k <= endtab[i]; k++)
				{
					int k2 = fixk(k, lengthmz);
					if (blur[index2D(numz, k2, j)] != 0)
					{
						int start = starttab[k2];
						cv += blur[index2D(numz, k2, j)] * mzdist[index2D(maxlength, k2, i-start)];
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
		#pragma omp parallel for private (i,j,k), schedule(dynamic)
		for (i = 0; i<lengthmz; i++)
		{
			for (j = 0; j<numz; j++)
			{

				double cv = 0;
				
				for (k = starttab[i]; k <= endtab[i]; k++)
				{
					if (blur[index2D(numz, k, j)] != 0)
					{
						cv += blur[index2D(numz, k, j)] * mzdist[indexmod(lengthmz, k, i)];
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


void PerfectTransform(int lengthmz, int numz, double *mtab, double massmax, double massmin, int maaxle, double *massaxis, double *massaxisval, double*blur, double *massgrid)
{
	int *indexes = NULL;
	indexes = calloc(lengthmz*numz, sizeof(double));
	int i, j;
	#pragma omp parallel for private(i,j), schedule(dynamic)
	for (i = 0; i < lengthmz; i++)
	{
		for (j = 0; j < numz; j++)
		{
			double testmass = mtab[index2D(numz, i, j)];
			if (testmass<massmax&&testmass>massmin) {
				indexes[index2D(numz, i, j)] = nearfast(massaxis, testmass, maaxle);
			}
		}
	}
	
	for (i = 0; i<lengthmz; i++)
	{
		for (j = 0; j<numz; j++)
		{
			int index = indexes[index2D(numz, i, j)];
			double testmass = mtab[index2D(numz, i, j)];
			if (testmass<massmax&&testmass>massmin) {
				int index = nearfast(massaxis, testmass, maaxle);
				double newval = blur[index2D(numz, i, j)];
				if (massaxis[index] == testmass) {
					massaxisval[index] += newval;
					massgrid[index2D(numz, index, j)] += newval;
				}

				if (massaxis[index] < testmass &&index<maaxle - 2)
				{
					int index2 = index + 1;
					double interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
					massaxisval[index] += (1.0 - interpos)*newval;
					massgrid[index2D(numz, index, j)] += (1.0 - interpos)*newval;
					massaxisval[index2] += (interpos)*newval;
					massgrid[index2D(numz, index2, j)] += (interpos)*newval;
				}

				if (massaxis[index] > testmass&&index>0)
				{
					int index2 = index - 1;
					double interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
					massaxisval[index] += (1 - interpos)*newval;
					massgrid[index2D(numz, index, j)] += (1 - interpos)*newval;
					massaxisval[index2] += (interpos)*newval;
					massgrid[index2D(numz, index2, j)] += (interpos)*newval;
				}
			}
		}
	}
	free(indexes);
}


void IntegrateTransform(int lengthmz, int numz,double *mtab,double massmax, double massmin, int maaxle, double *massaxis,double *massaxisval,double*blur,double *massgrid)
{
	for(int i=0;i<lengthmz;i++)
		{
			for(int j=0;j<numz;j++)
			{
			double testmass=mtab[index2D(numz,i,j)];
			if(testmass<massmax&&testmass>massmin){
				int index=nearfast(massaxis,testmass,maaxle);
				double newval=blur[index2D(numz,i,j)];
				if (massaxis[index] == testmass){
					massaxisval[index]+=newval;
					massgrid[index2D(numz,index,j)]+=newval;
				}
				
				if (massaxis[index] < testmass &&index<maaxle-2)
				{
					int index2 = index + 1;
					double interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
					massaxisval[index] += (1.0 - interpos)*newval;
					massgrid[index2D(numz, index, j)] += (1.0 - interpos)*newval;
					massaxisval[index2] += (interpos)*newval;
					massgrid[index2D(numz, index2, j)] += (interpos)*newval;
				}
				
				if (massaxis[index] > testmass&&index>0)
				{
					int index2 = index - 1;
					double interpos = LinearInterpolatePosition(massaxis[index], massaxis[index2], testmass);
					massaxisval[index] += (1 - interpos)*newval;
					massgrid[index2D(numz, index, j)] += (1 - interpos)*newval;
					massaxisval[index2] += (interpos)*newval;
					massgrid[index2D(numz, index2, j)] += (interpos)*newval;
				}			
			}
			}
		}
}

void InterpolateTransform(int maaxle,int numz,int lengthmz, int *nztab,double *massaxis,double adductmass,double *dataMZ,double *dataInt,double *massgrid, double *massaxisval, double *blur)
{
	double startmzval = dataMZ[0];
	double endmzval = dataMZ[lengthmz - 1];
	//#pragma omp parallel for schedule(dynamic)
	for(int i=0;i<maaxle;i++)
	{
		double val=0;
		double newval=0;
		for(int j=0;j<numz;j++)
		{
			double mztest=(massaxis[i]+nztab[j]*adductmass)/(double)nztab[j];

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
						double mu=(mztest-dataMZ[index])/(dataMZ[index2]-dataMZ[index]);
						double y0=blur[index2D(numz,index-1,j)];
						double y1=blur[index2D(numz,index,j)];
						double y2=blur[index2D(numz,index2,j)];
						double y3=blur[index2D(numz,index2+1,j)];
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




void KillB(double *I, char *B, double intthresh, int lengthmz, int numz, const int isolength, int *isotopepos, float *isotopeval)
{
	unsigned int i, j, k;
	if (isolength == 0){
		for (i = 0; i < lengthmz; i++)
		{
			for (j = 0; j < numz; j++)
			{
				if (I[i] < intthresh){ B[index2D(numz, i, j)] = 0; }
			}
		}
	}
	else
	{
		float cutoff = 0.5;
		printf("Removing species where isotope fall below %f\n", cutoff * 100);
#pragma omp parallel for private (i,j,k), schedule(dynamic)
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
						if (I[pos] < intthresh){ B[index2D(numz, i, j)] = 0; }
					}
				}

			}
		}
	}
}


double isotopemid(double mass,double *isoparams)
{
double a,b,c;
a=isoparams[4];
b=isoparams[5];
c=isoparams[6];
return a+b*pow(mass,c);
}

double isotopesig(double mass,double *isoparams)
{
double a,b,c;
a=isoparams[7];
b=isoparams[8];
c=isoparams[9];
return a+b*pow(mass,c);
}

double isotopealpha(double mass,double *isoparams)
{
double a,b;
a=isoparams[0];
b=isoparams[1];
return a*exp(-mass*b);
}

double isotopebeta(double mass,double *isoparams)
{
double a,b;
a=isoparams[2];
b=isoparams[3];
return a*exp(-mass*b);
}

int setup_isotopes(double *isoparams,int *isotopepos,float *isotopeval,double *mtab,int *ztab,char *barr,double *dataMZ,int lengthmz,int numz)
{
    double minmass=100000000;
    double maxmass=1;
    int i;
    for(i=0;i<lengthmz*numz;i++)
    {
        if(barr[i]==1)
        {
			double mass = mtab[i];
            if(mass<minmass){minmass=mass;}
            if(mass>maxmass){maxmass=mass;}
        }
    }

    double minmid=isotopemid(minmass,isoparams);
    double minsig=isotopesig(minmass,isoparams);
    double maxmid=isotopemid(maxmass,isoparams);
    double maxsig=isotopesig(maxmass,isoparams);

    int isostart=(int)(minmid-4*minsig);
    int isoend=(int)(maxmid+4*maxsig);
    if(isostart<0){isostart=0;}
	if (isoend < 4){ isoend = 4; }
    int isolength=isoend-isostart;
	return isolength;
}

void make_isotopes(double *isoparams, int *isotopepos, float *isotopeval, double *mtab, int *ztab, char *barr, double *dataMZ, int lengthmz, int numz)
{
	double minmass = 100000000;
	double maxmass = 1;
	int i, j, k;
	for (i = 0; i<lengthmz*numz; i++)
	{
		if (barr[i] == 1)
		{
			double mass = mtab[i];
			if (mass<minmass){ minmass = mass; }
			if (mass>maxmass){ maxmass = mass; }
		}
	}
	double massdiff = 1.0026;

	double minmid = isotopemid(minmass, isoparams);
	double minsig = isotopesig(minmass, isoparams);
	double maxmid = isotopemid(maxmass, isoparams);
	double maxsig = isotopesig(maxmass, isoparams);

	//int isostart = (int)(minmid - 4 * minsig);
	int isostart = 0;
	int isoend = (int)(maxmid + 4 * maxsig);
	//if (isostart<0){ isostart = 0; }
	if (isoend < 4){ isoend = 4; }
	int isolength = isoend - isostart;
	double *isorange = NULL;
	int *isoindex = NULL;
	isorange = calloc(isolength, sizeof(double));
	isoindex = calloc(isolength, sizeof(int));
	for (i = 0; i<isolength; i++)
	{
		isorange[i] = (isostart + i)*massdiff;
		isoindex[i] = (isostart + i);
	}
	#pragma omp parallel for private (i,j,k), schedule(dynamic)
	for (i = 0; i<lengthmz; i++)
	{
		for (j = 0; j<numz; j++)
		{
			if (barr[index2D(numz, i, j)] == 1)
			{
				double mz = dataMZ[i];
				int z = ztab[j];
				for (k = 0; k<isolength; k++)
				{
					double newmz = mz + (isorange[k] / ((double)z));
					int pos = nearfast(dataMZ, newmz, lengthmz);
					isotopepos[index3D(numz, isolength, i, j, k)] = pos;
				}
			}
		}
	}
	#pragma omp parallel for private (i,j,k), schedule(dynamic)
	for (i = 0; i<lengthmz; i++)
	{
		for (j = 0; j<numz; j++)
		{
			if (barr[index2D(numz, i, j)] == 1)
			{
				double mass = mtab[index2D(numz, i, j)];
				double mid = isotopemid(mass, isoparams);
				double sig = isotopesig(mass, isoparams);
				double alpha = isotopealpha(mass, isoparams);
				double amp = (1.0 - alpha) / (sig*2.50662827);
				double beta = isotopebeta(mass, isoparams);
				double tot = 0;
				for (k = 0; k<isolength; k++)
				{
					double e = alpha*exp(-isoindex[k]*beta);
					double g = amp *exp(-pow(isoindex[k] - mid, 2) / (2 * pow(sig, 2)));
					double temp = e + g;
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

void isotope_dist(double mass,int isolength, int *isoindex, double *isovals, double *isoparams)
{
	double mid = isotopemid(mass, isoparams);
	double sig = isotopesig(mass, isoparams);
	double alpha = isotopealpha(mass, isoparams);
	double amp = 1.0 - alpha;
	double beta = isotopebeta(mass, isoparams);
	double tot = 0;
	int k;
	for (k = 0; k<isolength; k++)
	{
		double e = alpha*exp(-isoindex[k] * beta);
		double g = amp / (sig*2.50662827)*exp(-pow(isoindex[k] - mid, 2) / (2 * pow(sig, 2)));
		tot += e + g;
		isovals[k] = e + g;
	}
	for (k = 0; k<isolength; k++)
	{
		if (tot > 0){ isovals[k] = isovals[k] / tot; }
	}
}

void textvectorprint(double *arr, int length)
{
	int levels = 20;
	double grad = 1.0 / ((double)levels);
	int i, j;
	printf("\n");
	double max = 0;
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

void test_isotopes(double mass, double *isoparams)
{
	int i;
	double maxmid = isotopemid(mass, isoparams);
	double maxsig = isotopesig(mass, isoparams);

	int isostart =0;
	int isoend = (int)(maxmid + 4 * maxsig);
	if (isoend < 4){ isoend = 4; }
	int isolength = isoend - isostart;
	double *isorange = NULL;
	int *isoindex = NULL;
	isorange = calloc(isolength, sizeof(double));
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

void deepcopy(double *out, const double *in, int len)
{
	for (int i = 0; i < len; i++)
	{
		out[i] = in[i];
	}
}

void simp_norm(const int length, double *data)
{
	double norm = Max(data, length);
	//printf("\nMax: %f %d\n", norm, length);
	if(norm>0){
		for (int i = 0; i<length; i++)
		{
			data[i] = data[i] / norm;
		}
	}
	return;
}



void charge_scaling(double *blur, const int *nztab, const int lengthmz, const int numz)
{
	for (int i = 0; i < lengthmz; i++)
	{
		for (int j = 0; j < numz; j++)
		{
			int charge = nztab[j];
			double z = (double)charge;
			blur[index2D(numz, i, j)] /= z;
		}
	}
	return; 
}
