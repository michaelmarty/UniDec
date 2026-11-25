/*
 * Coupled chromatographic UniDec deconvolution for MetaUniDec HDF5 files.
 */

#ifndef UNICHROM_MAIN_H
#define UNICHROM_MAIN_H

#include "MetaUniDec_Main.h"

/*
 * config.dtsig is the chromatography peak width in scans (UCtype=0) or
 * retention-time units (UCtype=1). Time mode always uses direct convolution.
 * UClineardecon=1 uses a common linear m/z grid and FFT convolution.
 * UClineardecon=0 keeps each processed m/z axis and uses direct convolution.
 */
int run_chromatogram(int argc, char *argv[], Config config);

int setup_blur_UniChrom(Config config, float *mz_axis, int mz_count,
                        const char *allowed, int batch_count, IntraDecon *intra);

void transform_mass_grid_UniChrom(Config config, const float *cube,
                                  const float *mz_axis, int scan_count,
                                  int mz_count, int charge_count,
                                  float mass_min, float mass_max,
                                  float *mass_axis, int mass_count,
                                  float *mass_grid);

#endif
