//
// Created by Michael Marty on 7/7/2025.
//

#include "udstruct.h"


Input SetupInputs() {
    Input inp;
    inp.dataMZ = NULL;
    inp.dataInt = NULL;
    inp.testmasses = NULL;
    inp.nztab = NULL;
    inp.mtab = NULL;
    inp.barr = NULL;
    inp.fwhmlist = NULL;
    return inp;
}

void FreeInputs(const Input inp) {
    free(inp.dataMZ);
    free(inp.dataInt);
    free(inp.nztab);
    free(inp.mtab);
    free(inp.testmasses);
    free(inp.barr);
    free(inp.fwhmlist);
}

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
    decon.starttab = NULL;
    decon.endtab = NULL;
    decon.mzdist = NULL;
    decon.rmzdist = NULL;
    decon.error = 0;
    decon.rsquared = 0;
    decon.iterations = 0;
    decon.uniscore = 0;
    decon.conv = 0;
    decon.threshold = 0;
    decon.mlen = 0;
    decon.plen = 0;
    decon.scanindex = 0;
    return decon;
}

void FreeDecon(const Decon decon) {
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
    free(decon.mzdist);
    free(decon.rmzdist);
    free(decon.starttab);
    free(decon.endtab);
}


void SetDefaultConfig(Config *config) {
    strcpy(config->infile, "default_data.txt");
    strcpy(config->outfile, "default_output");
    strcpy(config->mfile, "default_mfile->txt");
    config->numit = 50;
    config->endz = 100;
    config->startz = 1;
    config->zsig = 1;
    config->psig = 1;
    config->beta = 0;
    config->mzsig = 1;
    config->msig = 0;
    config->molig = 0;
    config->massub = 5000000;
    config->masslb = 100;
    config->psfun = 0;
    config->zpsfun = 0;
    config->psmzthresh = 0;
    config->mtabsig = 0;
    config->mflag = 0;
    config->massbins = 100;
    config->limitflag = 0;
    config->psthresh = 6;
    config->speedyflag = 0;
    config->aggressiveflag = 0;
    config->adductmass = 1.007276467f;
    config->rawflag = 1;
    config->nativezub = 100;
    config->nativezlb = -200;
    config->poolflag = 1;
    config->manualflag = 0;
    config->intthresh = 0;
    config->peakshapeinflate = 1;
    config->fixedmassaxis = 0;
    config->isotopemode = 0;
    config->filetype = 0;
    config->imflag = 0;
    config->linflag = -1;
    //IM Parameters
    config->dtsig = 0.2f;
    config->csig = 1;
    config->ccsub = 20000;
    config->ccslb = -20000;
    config->ccsbins = 100;
    config->temp = 0;
    config->press = 2;
    config->volt = 50;
    config->tcal1 = 0.3293f;
    config->tcal2 = 6.3597f;
    config->tcal3 = 0;
    config->tcal4 = 0;
    config->twaveflag = -1;
    config->hmass = 4.002602f;
    config->to = 0.97f;
    config->len = 0.18202f;
    config->edc = 1.57f;
    config->nativeccsub = 20000;
    config->nativeccslb = -20000;
    config->zout = 0;
    config->baselineflag = 1;
    config->noiseflag = 0;
    config->metamode = -2;
    config->minmz = -1;
    config->maxmz = -1;
    config->mzres = -1;
    config->mzbins = 0;
    config->bsub = 0;
    config->datareduction = 0;
    config->peakwin = 500;
    config->peakthresh = 0.1f;
    config->exchoice = 0;
    config->exchoicez = 1;
    config->exthresh = 10;
    config->exnorm = 0;
    config->exnormz = 0;
    config->peaknorm = 1;
    config->exwindow = 0;
    config->orbimode = 0;
    config->datanorm = 1;
    //Experimental
    config->filterwidth = 20;
    config->zerolog = -12;
    config->lengthmz = 0;
    config->mfilelen = 0;
    config->isolength = 0;
    // DoubleDec
    config->doubledec = 0; // Consider initializing kernel as well?
    config->file_id = 0;
    // Other
    config->silent = 0;
    config->cdmsflag = 0;
    config->variablepw = 0;
    config->minratio = 0.0f;
}

void PostImport(Config *config) {
    //Convert gaussian FWHM to sigma
    if (config->psfun == 0) { config->mzsig = config->mzsig / 2.35482f; }
    if (config->zpsfun == 0) { config->csig = config->csig / 2.35482f; }
    config->dtsig = config->dtsig / 2.35482f;
    //Check whether to turn off or on the config->limitflag. Limit flag on means it will select only the mass values closest to the mass values from mfile.
    if (config->mflag == 1 && config->mtabsig == 0) { config->limitflag = 1; }

    config->numz = config->endz - config->startz + 1;

    //If linflag is active, overwrite speedy}
    if (config->linflag != -1) {
        if (config->linflag != 2) { config->speedyflag = 1; } else { config->speedyflag = 0; }
    }
    if (config->speedyflag == 1) { printf("Speedy mode: Assuming Linearized Data\n"); }
    //
    // if (config->filetype == 0) {
    //     //Print inputs for check
    //     printf("infile = %s\n", config->infile);
    //     //printf("outfile = %s\n", config->outfile);
    //     //printf("\n");
    // }

    //Check to see if the mass axis should be fixed
    if (config->massub < 0 || config->masslb < 0) {
        config->fixedmassaxis = 1;
        config->massub = fabsf(config->massub);
        config->masslb = fabsf(config->masslb);
    }

    if (config->twaveflag == -1 && config->imflag == 1) {
        printf("\n\nNeed to define twaveflag for CCS calculation\n\n");
    }

    //bug proofing so we don't get a 1/0 problem
    //if (config->msig == 0) { config->msig = 0.00001; }
    //if (config->zsig == 0) { config->zsig = 0.00001; }
    if (config->massbins == 0) { config->massbins = 1; }

    //Convert aggressiveflag to baselineflag
    if (config->aggressiveflag == 1 || config->aggressiveflag == 2) { config->baselineflag = 1; } else {
        config->baselineflag = 0;
    }

    //Experimental correction. Not sure why this is necessary.
    if (config->psig < 0) { config->mzsig /= 3; }

    //Sets a threshold for m/z values to check. Things that are far away in m/z space don't need to be considered in the iterations.
    config->psmzthresh = config->psthresh * fabsf(config->mzsig) * config->peakshapeinflate;
}

Config LoadConfig(Config config, const char *filename) {
    // We assume argv[1] is a filename to open
    FILE *file;
    file = fopen(filename, "r");
    if (file == 0) {
        printf("Error Opening Configuration File: %s\n", filename);
        exit(1);
    }

    //Read parameters from configuration file
    //Configuration file should be formatted: name value
    else {
        char x[501] = {0};
        char y[501] = {0};
        char * endptr;
        //printf("\nRead from file:");
        while (fscanf(file, "%s %500[^\n]", x, y) != EOF) {
            //printf( "read in: %s %s \n", x,y );
            if (strstr(x, "input") != NULL) { strcpy(config.infile, y); } // printf(" input");
            if (strstr(x, "output") != NULL) { strcpy(config.outfile, y); } //printf(" %s output\n", config.outfile); }
            if (strstr(x, "mfile") != NULL) {
                strcpy(config.mfile, y);
                config.mflag = 1;
            } // printf(" mfile"); }
            if (strstr(x, "numit") != NULL) { config.numit = strtol(y, &endptr, 10); } // printf(" numit"); }
            //if (strstr(x, "numz") != NULL){ config.numz = strtol(y, &endptr, 10); printf(" numz"); }
            if (strstr(x, "startz") != NULL) { config.startz = strtol(y, &endptr, 10); } // printf(" startz"); }
            if (strstr(x, "endz") != NULL) { config.endz = strtol(y, &endptr, 10); } // printf(" endz"); }
            if (strstr(x, "zzsig") != NULL) { config.zsig = strtof(y, &endptr); } // printf(" zzsig"); }
            if (strstr(x, "psig") != NULL) { config.psig = strtof(y, &endptr); } // printf(" psig"); }
            if (strstr(x, "beta") != NULL) { config.beta = strtof(y, &endptr); } // printf(" beta"); }
            if (strstr(x, "mzsig") != NULL) { config.mzsig = strtof(y, &endptr); } // printf(" mzsig"); }
            if (strstr(x, "msig") != NULL) { config.msig = strtof(y, &endptr); } // printf(" msig"); }
            if (strstr(x, "molig") != NULL) { config.molig = strtof(y, &endptr); } // printf(" molig"); }
            if (strstr(x, "massub") != NULL) { config.massub = strtof(y, &endptr); } // printf(" massub"); }
            if (strstr(x, "masslb") != NULL) { config.masslb = strtof(y, &endptr); } // printf(" masslb"); }
            if (strstr(x, "psfun") != NULL) { config.psfun = strtol(y, &endptr, 10); } // printf(" psfun"); }
            if (strstr(x, "zpsfn") != NULL) { config.zpsfun = strtol(y, &endptr, 10); }
            if (strstr(x, "mtabsig") != NULL) { config.mtabsig = strtof(y, &endptr); } // printf(" mtabsig"); }
            if (strstr(x, "massbins") != NULL) { config.massbins = strtof(y, &endptr); } // printf(" massbins"); }
            if (strstr(x, "psthresh") != NULL) { config.psthresh = strtof(y, &endptr); } // printf(" psthresh"); }
            if (strstr(x, "speedy") != NULL) { config.speedyflag = strtol(y, &endptr, 10); } //  printf(" speedy"); }
            if (strstr(x, "aggressive") != NULL) { config.aggressiveflag = strtol(y, &endptr, 10); } //  printf(" aggressive"); }
            if (strstr(x, "adductmass") != NULL) { config.adductmass = strtof(y, &endptr); } //  printf(" adductmass"); }
            if (strstr(x, "rawflag") != NULL) { config.rawflag = strtol(y, &endptr, 10); } // printf(" rawflag"); }
            if (strstr(x, "nativezub") != NULL) { config.nativezub = strtof(y, &endptr); } // printf(" nativezub"); }
            if (strstr(x, "nativezlb") != NULL) { config.nativezlb = strtof(y, &endptr); } // printf(" nativezlb"); }
            if (strstr(x, "poolflag") != NULL) { config.poolflag = strtol(y, &endptr, 10); } // printf(" poolflag"); }
            if (strstr(x, "manualfile") != NULL) {
                config.manualflag = 1;
                strcpy(config.manualfile, y);
            } //  printf(" manualfile"); }
            if (strstr(x, "intthresh") != NULL) { config.intthresh = strtof(y, &endptr); } // printf(" intthresh"); }
            if (strstr(x, "peakshapeinflate") != NULL) { config.peakshapeinflate = strtof(y, &endptr); }
            // printf(" peakshapeinflate"); }
            if (strstr(x, "isotopemode") != NULL) { config.isotopemode = strtol(y, &endptr, 10); } // printf(" isotopemode"); }
            if (strstr(x, "orbimode") != NULL) { config.orbimode = strtol(y, &endptr, 10); } // printf(" orbimode"); }
            if (strstr(x, "imflag") != NULL) { config.imflag = strtol(y, &endptr, 10); } // printf(" imflag"); }
            if (strstr(x, "cdmsflag") != NULL) { config.cdmsflag = strtol(y, &endptr, 10); } // printf(" imflag"); }
            if (strstr(x, "linflag") != NULL) { config.linflag = strtol(y, &endptr, 10); } // printf(" linflag"); }
            //IM Parameters
            if (strstr(x, "csig") != NULL) { config.csig = strtof(y, &endptr); } // printf(" csig"); }
            if (strstr(x, "dtsig") != NULL) { config.dtsig = strtof(y, &endptr); } // printf(" dtsig"); }
            if (strstr(x, "ccsub") != NULL) { config.ccsub = strtof(y, &endptr); } // printf(" ccsub"); }
            if (strstr(x, "ccslb") != NULL) { config.ccslb = strtof(y, &endptr); } // printf(" ccslb"); }
            if (strstr(x, "ccsbins") != NULL) { config.ccsbins = strtof(y, &endptr); } // printf(" ccsbins"); }
            if (strstr(x, "temp") != NULL) { config.temp = strtof(y, &endptr); } // printf(" temp"); }
            if (strstr(x, "pressure") != NULL) { config.press = strtof(y, &endptr); } // printf(" pressure"); }
            if (strstr(x, "volt") != NULL) { config.volt = strtof(y, &endptr); } // printf(" volt"); }
            if (strstr(x, "gasmass") != NULL) { config.hmass = strtof(y, &endptr); } // printf(" gasmass"); }
            if (strstr(x, "tnaught") != NULL) { config.to = strtof(y, &endptr); } // printf(" to"); }
            if (strstr(x, "tcal1") != NULL) { config.tcal1 = strtof(y, &endptr); } // printf(" tcal1"); }
            if (strstr(x, "tcal2") != NULL) { config.tcal2 = strtof(y, &endptr); } // printf(" tcal2"); }
            if (strstr(x, "tcal3") != NULL) { config.tcal3 = strtof(y, &endptr); } // printf(" tcal3"); }
            if (strstr(x, "tcal4") != NULL) { config.tcal4 = strtof(y, &endptr); } // printf(" tcal4"); }
            if (strstr(x, "edc") != NULL) { config.edc = strtof(y, &endptr); } // printf(" edc"); }
            if (strstr(x, "zout") != NULL) { config.zout = strtol(y, &endptr, 10); } // printf(" zout"); }
            if (strstr(x, "twaveflag") != NULL) { config.twaveflag = strtol(y, &endptr, 10); } // printf(" twaveflag"); }
            if (strstr(x, "ubnativeccs") != NULL) { config.nativeccsub = strtof(y, &endptr); } // printf(" ubnativeccs"); }
            if (strstr(x, "lbnativeccs") != NULL) { config.nativeccslb = strtof(y, &endptr); } // printf(" lbnativeccs"); }
            if (strstr(x, "driftlength") != NULL) { config.len = strtof(y, &endptr); } // printf(" driftlength"); }
            if (strstr(x, "baselineflag") != NULL) { config.baselineflag = strtol(y, &endptr, 10); } // printf(" baselineflag"); }
            if (strstr(x, "noiseflag") != NULL) { config.noiseflag = strtol(y, &endptr, 10); } // printf(" noiseflag"); }
            //Experimental
            if (strstr(x, "filterwidth") != NULL) { config.filterwidth = strtol(y, &endptr, 10); } // printf(" filterwidth"); }
            if (strstr(x, "zerolog") != NULL) { config.zerolog = strtof(y, &endptr); } // printf(" zerolog"); }
            //Peak Parameters
            if (strstr(x, "peakwindow") != NULL) { config.peakwin = strtof(y, &endptr); } // printf(" peakwindow"); }
            if (strstr(x, "peakthresh") != NULL) { config.peakthresh = strtof(y, &endptr); } // printf(" peakthresh"); }
            if (strstr(x, "peaknorm") != NULL) { config.peaknorm = strtol(y, &endptr, 10); } // printf(" peaknorm"); }
            // DoubleDec Parameters
            if (strstr(x, "doubledec") != NULL) { config.doubledec = strtol(y, &endptr, 10); }
            if (strstr(x, "kernel") != NULL) { strcpy(config.kernel, y); }
        }
        //printf("\n\n");
    }
    fclose(file);

    PostImport(&config);

    return config;
}


//Gives command line help options
void PrintHelp() {
    printf(
        "\nUniDec: Universal Deconvolution of Mass and Ion Mobility Spectra\n\nWritten by Michael Marty\n\twith contributions from Erik Marklund and Andrew Baldwin\n");
    printf("Copyright University of Oxford 2016 and University of Arizona 2017\n");
    printf("\nUniDec runs off a configuration text or hdf5 file.\nExample usage: UniDec.exe conf.txt\n");
    printf("Text configuration file should be written for each line: keyword argument\n");
    printf(
        "HDF5 files should consult the Python API source code for a model of directory and metadata construction.\n");
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
    printf(
        "\t\"psfun\" \tPeak shape function\n\t\t\t\t0=Gaussian\n\t\t\t\t1=Lorentzian\n\t\t\t\t2=Split Gaussian\\Lorentzian\n");
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
    printf(
        "\t\"isotopemode\" \t0=off 1=monoisotopic 2=average: Uses isotope distributions in deconvolution (MS Only)\n");
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
    printf("\t\t\t\t4=SLIM T-Wave 2nd Order Polynomial Calibration\n");
    printf("\t\t\t\t\t\"tcal3\"=Calibration paramter 3 (P3)\n");
    printf("\t\t\t\t5=SLIM T-Wave 3rd Order Polynomial Calibration\n");
    printf("\t\t\t\t\t\"tcal4\"=Calibration parmater 4 (P4)\n");
    printf("\nEnjoy! Please report bugs to Michael Marty (mtmarty@utexas.edu) commit date 7/7/2025\n");
    //printf("\nsize of: %d",sizeof(char));

    /*
    float a = 119.528155;
    float b = 0.32967;
    float c = -0.000115585597;
    float d = 2.4120647e-8;
    float mass = 760;
    float z = 1;

    float edc = 0;
    float dt = 778;
    float ccs = 310;
    float hmass = 28;

    a = 153.38;
    b = 0.58873;
    c = -3.4613e-4;
    d = 1.2183e-7;


    mass = 50000;
    z = 15;
    dt = 350;
    ccs = 4755;

    printf("Mass: %f Charge: %f GasM: %f\n", mass, z, hmass);
    printf("a=%f b=%f c=%f d=%e\n", a, b, c, d);


    float test2 = calcCCSSLIMpoly2(mass, z, dt, a, b, c, hmass, edc);
    printf("Test 2nd Order DT: %f => CCS: %f\n", dt, test2);
    test2 = calcDtSLIMpoly2(mass, z, ccs, a, b, c, hmass, edc);
    printf("Test 2nd Order CCS: %f => DT: %f\n", ccs, test2);

    float test=calcCCSSLIMpoly3(mass, z, dt, a, b, c, d, hmass, edc);
    printf("Test 3rd Order DT: %f => CCS: %f\n", dt, test);
    test = calcDtSLIMpoly3(mass, z, ccs, a, b, c, d, hmass, edc);
    printf("Test 3rd Order CCS: %f => Dt: %f\n", ccs, test);*/
}

void PrintConfig(Config config) {
    printf("infile: %s\n", config.infile);
    printf("outfile: %s\n", config.outfile);
    printf("mfile: %s\n", config.mfile);
    printf("numit: %d\n", config.numit);
    printf("startz: %d\n", config.startz);
    printf("endz: %d\n", config.endz);
    printf("zsig: %f\n", config.zsig);
    printf("psig: %f\n", config.psig);
    printf("beta: %f\n", config.beta);
    printf("mzsig: %f\n", config.mzsig);
    printf("msig: %f\n", config.msig);
    printf("molig: %f\n", config.molig);
    printf("massub: %f\n", config.massub);
    printf("masslb: %f\n", config.masslb);
    printf("psfun: %d\n", config.psfun);
    printf("zpsfun: %d\n", config.zpsfun);
    printf("mtabsig: %f\n", config.mtabsig);
    printf("massbins: %f\n", config.massbins);
    printf("psthresh: %f\n", config.psthresh);
    printf("speedyflag: %d\n", config.speedyflag);
    printf("aggressiveflag: %d\n", config.aggressiveflag);
    printf("adductmass: %f\n", config.adductmass);
    printf("rawflag: %d\n", config.rawflag);
    printf("nativezub: %f\n", config.nativezub);
    printf("nativezlb: %f\n", config.nativezlb);
    printf("poolflag: %d\n", config.poolflag);
    printf("manualflag: %d manualfile:%s \n", config.manualflag, config.manualfile);
    printf("intthresh:%f \n", config.intthresh);
    printf("peakshapeinflate:%f \n", config.peakshapeinflate);
    //IM Parameters
    if (config.imflag == 1) {
        printf("\nIM-MS Parameters:\n");
        printf("\tcsig:%f \n\tdtsig:%f \n\tccsub:%f \n\tccslb:%f \n\tccsbins:%f \n",
               config.csig, config.dtsig, config.ccsub, config.ccslb, config.ccsbins);

        printf("\ttemp:%f \n\tpress:%f \n\tvolt:%f \n\ttnaught:%f \n",
               config.temp, config.press, config.volt, config.to);

        printf("\ttcal1:%f \n\ttcal2:%f \n\ttcal3:%f \n\ttcal4:%f \n",
               config.tcal1, config.tcal2, config.tcal3, config.tcal4);

        printf("\thmass:%f \n\tedc:%f \n", config.hmass, config.edc);

        printf("\ttwaveflag:%d \n\tubnativeccs:%f \n\tlbnativeccs:%f \n",
               config.twaveflag, config.nativeccsub, config.nativeccslb);
    }

    printf("isotopemode: %d\n", config.isotopemode);
    printf("orbimode: %d\n", config.orbimode);
    printf("imflag: %d\n", config.imflag);
    printf("cdmsflag: %d\n", config.cdmsflag);
    printf("linflag: %d\n", config.linflag);
    printf("fixedmassaxis: %d\n", config.fixedmassaxis);
    printf("filetype: %d\n", config.filetype);
    //Experimental
    printf("\nExperimental Parameters:\n");
    printf("\tfilterwidth: %d \n\tzerolog: %.1e \n\tlengthmz: %d \n\tmfilelen: %d \n\tisolength: %d \n",
           config.filterwidth, config.zerolog, config.lengthmz, config.mfilelen, config.isolength);
    printf("\tdoubledec:%d \n", config.doubledec);
    printf("\tkernel:%s \n", config.kernel);
    printf("\tsilent:%d \n", config.silent);
    printf("\tbaselineflag: %d \n", config.baselineflag);
    printf("\tnoiseflag: %d \n", config.noiseflag);
    printf("\tmetamode: %d \n", config.metamode);
    printf("\tminmz: %f \n", config.minmz);
    printf("\tmaxmz: %f \n", config.maxmz);
    printf("\tmzres: %f \n", config.mzres);
    printf("\tmzbins: %f \n", config.mzbins);
    printf("\tbsub: %f \n", config.bsub);
    printf("\tdatareduction: %f \n", config.datareduction);
    printf("\tpeakwin: %f \n", config.peakwin);
    printf("\tpeakthresh: %.2f \n", config.peakthresh);
    printf("\texchoice: %d \n", config.exchoice);
    printf("\texchoicez: %d \n", config.exchoicez);
    printf("\texthresh: %f \n", config.exthresh);
    printf("\texnorm: %d \n", config.exnorm);
    printf("\texnormz: %d \n", config.exnormz);
    printf("\texwindow: %f \n", config.exwindow);
    printf("\tdatanorm: %d \n", config.datanorm);
    printf("\tsilent: %d \n", config.silent);
    printf("\n");
}