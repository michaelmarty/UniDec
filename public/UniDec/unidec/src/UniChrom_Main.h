/*
 * Coupled chromatographic UniDec deconvolution for MetaUniDec HDF5 files.
 */

#ifndef UNICHROM_MAIN_H
#define UNICHROM_MAIN_H

#include "MetaUniDec_Main.h"
#include "UCCD_Main.h"

/*
 * config.dtsig is the chromatography peak width in scans.
 * UClineardecon=1 uses a common linear m/z grid and FFT convolution.
 * UClineardecon=0 keeps each processed m/z axis and uses direct convolution.
 */
int run_chromatogram(int argc, char *argv[], Config config);

#endif
