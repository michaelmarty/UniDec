/*
 * UCCD_Main.h
 *
 * Charge-detection deconvolution with an additional chromatography axis.
 */

#ifndef UCCD_MAIN_H
#define UCCD_MAIN_H

#include <math.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "udcore.h"

/*
 * UCCD data are stored in a sparse little-endian binary format. The header is:
 *
 *     char[8] magic ("UCCDBIN1")
 *     uint32 chromatography_count, mz_count, charge_count
 *     uint64 nonzero_count
 *     float32 chromatography_axis[], mz_axis[], charge_axis[]
 *     { uint32 flat_index; float32 intensity; } nonzero_values[]
 *
 * flat_index addresses dense [chromatography][m/z][charge] storage, with
 * charge contiguous.
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

/*
 * config.dtsig is the chromatographic peak width, expressed in chromatography-
 * coordinate units. config.mzsig and config.csig are the m/z and charge peak
 * widths, respectively. Zero disables broadening on the corresponding axis.
 * Iterative charge-state smoothing remains controlled by config.zsig.
 */
int run_unidec_UCCD(int argc, char *argv[], Config config);

#endif
