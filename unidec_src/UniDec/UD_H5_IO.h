/*
* UD_H5_IO.h
*
*  Created on : 19 December 2016
* Author : Michael.Marty
*/

//
//
// Copyright 2016 University of Arizona
//
//

#ifndef UD_H5_IO
#define UD_H5_IO

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int check_group(hid_t file_id, char *dataname)
{
	int status = H5Lexists(file_id, dataname, H5P_DEFAULT);
	if (status == 0) { printf("Dataset %s does not exist. Error.\n", dataname); exit(12); }
	return status;
}

int check_group_noexit(hid_t file_id, char* dataname)
{
	int status = H5Lexists(file_id, dataname, H5P_DEFAULT);
	if (status == 0) { printf("Dataset %s does not exist.\n", dataname);}
	return status;
}

void delete_group(hid_t file_id, char* dataname)
{
	if (H5LTpath_valid(file_id, dataname, 1)) {
		H5Ldelete(file_id, dataname, H5P_DEFAULT);
	}
}


int mh5getfilelength(hid_t file_id, char *dataname)
{
	check_group(file_id, dataname);
	hsize_t     dims[2];
	H5LTget_dataset_info(file_id, dataname, dims, NULL, NULL);
	// printf("%d %d\n", dims[0], dims[1]);
	return dims[0];
}

/*
void mh5readfile2d_grid(hid_t file_id, char* dataname, int length1, int length2, float* data1)
{
	//Unfinished
	const hsize_t length[2] = { length1 ,length2 };
	if (H5LTpath_valid(file_id, dataname, 1)) {
		H5Ldelete(file_id, dataname, H5P_DEFAULT);
	}
	H5LTmake_dataset_float(file_id, dataname, 2, length, data1);
}*/

void mh5readfile2d(hid_t file_id, char *dataname, int lengthmz, float *dataMZ, float *dataInt)
{
	check_group(file_id, dataname);
	int i, j;
	float *data;
	data = calloc(lengthmz * 2, sizeof(float));
	H5LTread_dataset_float(file_id, dataname, data);
	for (i = 0; i<lengthmz; i++)
	{
		dataMZ[i] = data[i * 2];
		dataInt[i] = data[i * 2 + 1];
	}
	free(data);
}

void mh5readfile3d(hid_t file_id, char *dataname, int lengthmz, float *dataMZ, float *dataInt, float *data3)
{
	check_group(file_id, dataname);
	int i, j;
	float *data;
	data = calloc(lengthmz * 3, sizeof(float));
	H5LTread_dataset_float(file_id, dataname, data);
	for (i = 0; i<lengthmz; i++)
	{
		dataMZ[i] = data[i * 3];
		dataInt[i] = data[i * 3 + 1];
		data3[i] = data[i * 3 + 2];
	}
	free(data);
}

void mh5readfile1d(hid_t file_id, char *dataname, float *data)
{
	check_group(file_id, dataname);
	H5LTread_dataset_float(file_id, dataname, data);
}

void mh5writefile1d(hid_t file_id, char *dataname,int length, float *data1)
{
	const hsize_t length2[1] = { length};
	if (H5LTpath_valid(file_id, dataname, 1)){
		H5Ldelete(file_id, dataname, H5P_DEFAULT);
	}
	H5LTmake_dataset_float(file_id, dataname, 1, length2, data1);
}

void mh5writefile2d(hid_t file_id, char *dataname, int length, float *data1, float *data2)
{
	const hsize_t length2[2] = { length ,2};
	if (H5LTpath_valid(file_id, dataname, 1)){
		H5Ldelete(file_id, dataname, H5P_DEFAULT);
	}
	float *data;
	int i;
	data = calloc(length * 2, sizeof(float));
	for (i = 0; i<length; i++)
	{
		data[2*i] = data1[i];
		data[2*i+1] = data2[i];
	}
	H5LTmake_dataset_float(file_id, dataname, 2, length2, data);
	free(data);
}

void mh5writefile2d_grid(hid_t file_id, char *dataname, int length1, int length2, float *data1)
{
	const hsize_t length[2] = { length1 ,length2 };
	if (H5LTpath_valid(file_id, dataname, 1)) {
		H5Ldelete(file_id, dataname, H5P_DEFAULT);
	}
	H5LTmake_dataset_float(file_id, dataname, 2, length, data1);
}

int int_attr(hid_t file_id, char *gname, char * attr_name, int value)
{
	hid_t dset_id = H5Gopen1(file_id, gname);
	herr_t x = H5LTfind_attribute(dset_id, attr_name);
	herr_t y;
	int attr[1] = { 0 };
	if (x == 1){
		y=H5LTget_attribute_int(file_id, gname, attr_name, attr);
		if (y < 0) { printf("Error with: %s\n", attr_name); }
		else{ value = attr[0]; }		
	}
	H5Gclose(dset_id);
	return value;
}

float float_attr(hid_t file_id, char *gname, char * attr_name, float value)
{
	hid_t dset_id = H5Gopen1(file_id, gname);
	herr_t x = H5LTfind_attribute(dset_id, attr_name);
	herr_t y;
	float attr[1] = { 0 };
	if (x == 1){
		y=H5LTget_attribute_float(file_id, gname, attr_name, attr);
		if (y < 0) { printf("Error with: %s\n", attr_name); }
		else { value = attr[0]; }
	}
	H5Gclose(dset_id);
	return value;
}

void write_attr_int(hid_t file_id, char *gname, char *attr_name, int value)
{
	int buff[1] = { value };
	H5LTset_attribute_int(file_id, gname, attr_name, buff, 1);
}

void write_attr_float(hid_t file_id, char *gname, char *attr_name, float value)
{
	float buff[1] = { value };
	H5LTset_attribute_float(file_id, gname, attr_name, buff, 1);
}

Config mh5LoadConfig(Config config, hid_t file_id)
{
	
	config.numit=int_attr(file_id, "/config", "numit", config.numit);
	config.endz=int_attr(file_id, "/config", "endz", config.endz);
	config.startz=int_attr(file_id, "/config", "startz", config.startz);
	config.psfun=int_attr(file_id, "/config", "psfun", config.psfun);
	config.linflag=int_attr(file_id, "/config", "linflag", config.linflag);
	config.aggressiveflag=int_attr(file_id, "/config", "aggressive", config.aggressiveflag);
	config.rawflag=int_attr(file_id, "/config", "rawflag", config.rawflag);
	config.poolflag=int_attr(file_id, "/config", "poolflag", config.poolflag);
	config.isotopemode=int_attr(file_id, "/config", "isotopemode", config.isotopemode);
	config.orbimode = int_attr(file_id, "/config", "orbimode", config.orbimode);
	config.manualflag=int_attr(file_id, "/config", "manualfileflag", config.manualflag);
	config.mflag = int_attr(file_id, "/config", "mfileflag", config.mflag);

	config.zsig=float_attr(file_id, "/config", "zzsig", config.zsig);
	config.psig = float_attr(file_id, "/config", "psig", config.psig);
	config.beta = float_attr(file_id, "/config", "beta", config.beta);
	config.mzsig=float_attr(file_id, "/config", "mzsig", config.mzsig);
	config.msig=float_attr(file_id, "/config", "msig", config.msig);
	config.molig=float_attr(file_id, "/config", "molig", config.molig);
	config.massub=float_attr(file_id, "/config", "massub", config.massub);
	config.masslb=float_attr(file_id, "/config", "masslb", config.masslb);
	config.mtabsig=float_attr(file_id, "/config", "mtabsig", config.mtabsig);
	config.psthresh=float_attr(file_id, "/config", "psthresh", config.psthresh);
	config.massbins=float_attr(file_id, "/config", "massbins", config.massbins);
	config.adductmass=float_attr(file_id, "/config", "adductmass", config.adductmass);
	config.nativezub=float_attr(file_id, "/config", "nativezub", config.nativezub);
	config.nativezlb=float_attr(file_id, "/config", "nativezlb", config.nativezlb);
	config.intthresh=float_attr(file_id, "/config", "intthresh", config.intthresh);
	config.peakshapeinflate=float_attr(file_id, "/config", "peakshapeinflate", config.peakshapeinflate);
	//IM Parameters
	config.csig =float_attr(file_id, "/config", "csig", config.csig);
	config.dtsig = float_attr(file_id, "/config", "dtsig", config.dtsig);
	config.ccsub = float_attr(file_id, "/config", "ccsub", config.ccsub);
	config.ccslb = float_attr(file_id, "/config", "ccslb", config.ccslb);
	config.ccsbins = float_attr(file_id, "/config", "ccsbins", config.ccsbins);
	config.temp = float_attr(file_id, "/config", "temp", config.temp);
	config.press = float_attr(file_id, "/config", "pressure", config.press);
	config.volt = float_attr(file_id, "/config", "volt", config.volt);
	config.hmass = float_attr(file_id, "/config", "gasmass", config.hmass);
	config.to = float_attr(file_id, "/config", "tnaught", config.to);
	config.tcal1 = float_attr(file_id, "/config", "tcal1", config.tcal1);
	config.tcal2 = float_attr(file_id, "/config", "tcal2", config.tcal2);
	config.edc = float_attr(file_id, "/config", "edc", config.edc);
	config.nativeccsub = float_attr(file_id, "/config", "ubnativeccs", config.nativeccsub);
	config.nativeccslb = float_attr(file_id, "/config", "lbnativeccs", config.nativeccslb);
	config.len = float_attr(file_id, "/config", "driftlength", config.len);

	config.zout = int_attr(file_id, "/config", "zout", config.zout);
	config.twaveflag = int_attr(file_id, "/config", "twaveflag", config.twaveflag);
	config.baselineflag = int_attr(file_id, "/config", "baselineflag", config.baselineflag);
	config.noiseflag = int_attr(file_id, "/config", "noiseflag", config.noiseflag);

	config.metamode = int_attr(file_id, "/config", "metamode", config.metamode);
	config.minmz = float_attr(file_id, "/config", "minmz", config.minmz);
	config.maxmz = float_attr(file_id, "/config", "maxmz", config.maxmz);
	config.mzbins = int_attr(file_id, "/config", "mzbins", config.mzbins);
	config.bsub = float_attr(file_id, "/config", "subbuff", config.bsub);
	config.datareduction = float_attr(file_id, "/config", "reductionpercent", config.datareduction);

	config.peakwin = float_attr(file_id, "/config", "peakwindow", config.peakwin);
	config.peakthresh = float_attr(file_id, "/config", "peakthresh", config.peakthresh);

	config.exwindow = float_attr(file_id, "/config", "exwindow", config.exwindow);
	config.exchoice = int_attr(file_id, "/config", "exchoice", config.exchoice);
	config.exchoicez = int_attr(file_id, "/config", "exchoicez", config.exchoicez);
	config.exthresh = float_attr(file_id, "/config", "exthresh", config.exthresh);
	config.exnorm = int_attr(file_id, "/config", "exnorm", config.exnorm);
	config.exnormz = int_attr(file_id, "/config", "exnormz", config.exnormz);
	config.peaknorm = int_attr(file_id, "/config", "peaknorm", config.peaknorm);
	config.datanorm = int_attr(file_id, "/config", "datanorm", config.datanorm);

	//Experimental
	config.filterwidth = int_attr(file_id, "/config", "filterwidth", config.filterwidth);
	config.zerolog = int_attr(file_id, "/config", "zerolog", config.zerolog);

	config=PostImport(config);

	return config;

}

void strjoin(const char *s1, const char *s2, char *newstring)
{
	strcpy(newstring, s1);
	strcat(newstring, s2);
}

void makegroup(hid_t file_id, char *dataset)
{
	hid_t group_id;
	htri_t status;
	status = H5Lexists(file_id, dataset, H5P_DEFAULT);
	if (status == 0) {
		group_id = H5Gcreate2(file_id, dataset, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		H5Gclose(group_id);
		//printf("Peaks Group not found. Creating...\n");
	}
}

void set_needs_grids(hid_t file_id)
{
	write_attr_int(file_id, "/config", "gridsflag", 0);
}

void set_got_grids(hid_t file_id)
{
	write_attr_int(file_id, "/config", "gridsflag", 1);
}

int question_grids(hid_t file_id)
{
	//int status = check_group(file_id, "/ms_dataset/mass_axis");
	int status = H5Lexists(file_id, "/ms_dataset/mass_axis", H5P_DEFAULT);
	if (status == 0) { write_attr_int(file_id, "/config", "gridsflag", 0); }//printf("Grids not found %d\n", status);}
	return int_attr(file_id, "/config", "gridsflag", 0);
}

int question_rawgrid(hid_t file_id)
{
	//int status = check_group(file_id, "/ms_dataset/mass_axis");
	int status = H5Lexists(file_id, "/ms_dataset/raw_axis", H5P_DEFAULT);
	if (status == 0) { write_attr_int(file_id, "/config", "gridsflag", 0); }//printf("Grids not found %d\n", status);}
	return int_attr(file_id, "/config", "gridsflag", 0);
}

#endif