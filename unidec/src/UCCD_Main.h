/*
 * UCCD_Main.h
 *
 * Charge-detection deconvolution with an additional chromatography axis.
 */

#ifndef UCCD_MAIN_H
#define UCCD_MAIN_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "udcore.h"

/*
 * UCCD data are dense and ordered as chromatography, m/z, charge.  The input
 * text file has four columns in that same order:
 *
 *     chromatography_coordinate  m/z  charge  intensity
 */

void blur_it_UCCD(float *output, const float *input, const int *upinds,
                  const int *loinds, int length, float floor);

void setup_blur_z_UCCD(int *zupind, int *zloind, const float *mzdat,
                       const float *zdat, int scan_length, float adductmass,
                       const float mzranges[4], const int size[3]);

void setup_blur_m_UCCD(int *mupind, int *mloind, const float *mzdat,
                       const float *zdat, int scan_length, float adductmass,
                       const float mzranges[4], const int size[3], float molig);

void make_kernel3D_UCCD(float *peak, const int size[3], const float *chromext,
                        const float *mzext, const float *zext, float chromsig,
                        float mzsig, float zsig, int psfun, int zpsfun);

void precompute_fft3D_UCCD(const float *list, const int size[3],
                           fftwf_complex *output);

void fftconvolve3D_precomputed_UCCD(float *corr, const float *list,
                                    const fftwf_complex *kernel_fft,
                                    const int size[3], fftwf_plan forward_plan,
                                    fftwf_plan backward_plan,
                                    fftwf_complex *work,
                                    fftwf_complex *transform);

void complex_conjugate_UCCD(const fftwf_complex *input, fftwf_complex *output,
                            int length);

/*
 * config.csig is the chromatographic peak width, expressed in chromatography-
 * coordinate units.  It uses the peak shape selected by config.psfun, just as
 * config.mzsig does for the m/z axis.  Zero disables chromatography broadening.
 * Charge-axis smoothing remains controlled by config.zsig.
 */
int run_unidec_UCCD(int argc, char *argv[], Config config);

#endif
