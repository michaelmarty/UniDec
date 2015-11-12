#include <stdio.h>
#include <stdlib.h>
#include <math.h>



// Written by MTM
//
// Helper functions for the python GUI
//
// Useful for speeding up some slow python parts
//


double mzpeakshape(double x,double y,double sig,int psfun)
{
    double result;
    if(psfun==0)
    {
		double sig2 = sig / 2.35482;
        result=exp(-(pow(x-y,2))/(2*sig2*sig2));
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
        //printf("Invalid Peak Function");
        exit(14);
    }
    return result;
}

int nearfast(double *dataMZ,double point, int numdat)
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

void convolve(double *xvals,double * input,double *output, int psfun,double sig,int len){
    int i,k;
    double window;
    if(psfun==0)
    {
    window=sig*5;
    }
    else{window=sig*15;}
    #pragma omp parallel for private(i,k), schedule(dynamic)
    for(i=0;i<len;i++)
    {
        double temp=0;
        double start=xvals[i]-window;
        double end=xvals[i]+window;
        int startpos=nearfast(xvals,start,len);
        int endpos=nearfast(xvals,end,len);
        for(k=startpos;k<endpos;k++)
        {
            //double diff=abs(xvals[k]-xvals[i]);
            //if(diff<window){
            temp+=input[k]*mzpeakshape(xvals[k],xvals[i],sig,psfun);
            //}
        }
        output[i]=temp;
    }
}


void distance(double *points, double * output, int len){
		int i, k;
		#pragma omp parallel for private(i,k), schedule(dynamic)
		for (i = 0; i < len; i++)
		{
			for (k = 0; k < len; k++)
			{
				double x = pow(points[i * 3] - points[k * 3], 2);
				double y = pow(points[i * 3 + 1] - points[k * 3 + 1], 2);
				double z = pow(points[i * 3 + 2] - points[k * 3 + 2], 2);
				//double result=sqrt(x+y+z);
				output[i*len + k] = sqrt(x + y + z);
				//output[k*len+i]=result;

			}
		}
	}


void inversemask(double *mzgrid,double *massgrid,double *rands,double *totals,double peakwidth,int psfun,int bnum,int gdim1,int gdim2)
{
int i,j,k;
#pragma omp parallel for private(i,j,k),schedule(dynamic)
for(i=0;i<bnum;i++){
    double blackout=rands[i];
    for(j=0;j<gdim1;j++){
        for(k=0;k<gdim2;k++){
            double mz=mzgrid[j*gdim2+k];
            if(fabs(mz-blackout)<peakwidth*5.0)
            {
                double val=mzpeakshape(blackout,mz,peakwidth,psfun);
                totals[i*gdim1+j]+=val*massgrid[j*gdim2+k];
            }
        }
    }
}
}

