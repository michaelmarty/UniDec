/*
 * Coupled chromatographic UniDec deconvolution for MetaUniDec HDF5 files.
 */

#ifndef UNICHROM_MAIN_H
#define UNICHROM_MAIN_H

#include "MetaUniDec_Main.h"
#include "UCCD_Main.h"

/*
 * config.dtsig is the chromatography peak width in scans.
 * The HDF5 spectra are coupled as [chromatography][m/z][charge] and processed
 * with three-dimensional FFTs.
 */
int run_chromatogram(int argc, char *argv[], Config config);

#endif
