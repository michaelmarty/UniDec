/*
 * Coupled chromatographic extension of UniDec.
 *
 * The observed MetaUniDec spectra are first interpolated onto one m/z axis by
 * make_grid().  The latent [chromatography][m/z][charge] cube is then updated
 * by Richardson-Lucy iterations.  Its forward model convolves the cube in
 * chromatography and m/z, sums charge, and compares that projection with the
 * observed [chromatography][m/z] grid.
 */

#include "UniChrom_Main.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define UNICHROM_OMP_MIN_LENGTH 32768

typedef struct {
    int real_length;
    int spectrum_length;
    int normalization;
    float *real_work;
    fftwf_complex *spectrum_work;
    fftwf_plan forward_plan;
    fftwf_plan backward_plan;
} UniChromFFT;


static int initialize_fft_UniChrom(UniChromFFT *context, const int size[3])
{
    memset(context, 0, sizeof(*context));
    const int64_t real_length = (int64_t)size[2] * size[0] * size[1];
    const int64_t spectrum_length = (int64_t)size[2] * size[0] * (size[1] / 2 + 1);
    if (real_length <= 0 || real_length > INT_MAX ||
        spectrum_length <= 0 || spectrum_length > INT_MAX) {
        fprintf(stderr, "UniChrom FFT cube is too large\n");
        return 0;
    }
    context->real_length = (int)real_length;
    context->spectrum_length = (int)spectrum_length;
    context->normalization = context->real_length;
    context->real_work = fftwf_malloc((size_t)context->real_length * sizeof(float));
    context->spectrum_work = fftwf_malloc(
        (size_t)context->spectrum_length * sizeof(fftwf_complex));
    if (context->real_work == NULL || context->spectrum_work == NULL) {
        fprintf(stderr, "Unable to allocate UniChrom FFT workspace\n");
        return 0;
    }
    context->forward_plan = fftwf_plan_dft_r2c_3d(
        size[2], size[0], size[1], context->real_work,
        context->spectrum_work, FFTW_ESTIMATE);
    context->backward_plan = fftwf_plan_dft_c2r_3d(
        size[2], size[0], size[1], context->spectrum_work,
        context->real_work, FFTW_ESTIMATE);
    if (context->forward_plan == NULL || context->backward_plan == NULL) {
        fprintf(stderr, "Unable to create UniChrom 3-D FFT plans\n");
        return 0;
    }
    return 1;
}


static void destroy_fft_UniChrom(UniChromFFT *context)
{
    if (context->forward_plan != NULL) { fftwf_destroy_plan(context->forward_plan); }
    if (context->backward_plan != NULL) { fftwf_destroy_plan(context->backward_plan); }
    fftwf_free(context->real_work);
    fftwf_free(context->spectrum_work);
    memset(context, 0, sizeof(*context));
}


static void execute_convolution_UniChrom(const float *input,
                                         const fftwf_complex *kernel_fft,
                                         const int conjugate_kernel,
                                         UniChromFFT *context)
{
    memcpy(context->real_work, input,
           (size_t)context->real_length * sizeof(float));
    fftwf_execute(context->forward_plan);
    #pragma omp parallel for schedule(static) if(context->spectrum_length >= UNICHROM_OMP_MIN_LENGTH)
    for (int i = 0; i < context->spectrum_length; i++) {
        const float a = context->spectrum_work[i][0];
        const float b = context->spectrum_work[i][1];
        const float c = kernel_fft[i][0];
        const float d = conjugate_kernel ? -kernel_fft[i][1] : kernel_fft[i][1];
        context->spectrum_work[i][0] = a * c - b * d;
        context->spectrum_work[i][1] = b * c + a * d;
    }
    fftwf_execute(context->backward_plan);
}


static int make_kernel_fft_UniChrom(fftwf_complex *kernel_fft,
                                    UniChromFFT *context, const int size[3],
                                    const float *mz_axis, const float *charge_axis,
                                    const Config config)
{
    /* Charge is latent rather than measured, so the instrumental kernel is a
     * delta on charge.  Charge/mass regularization is applied separately. */
    make_kernel3D_UCCD(context->real_work, size, mz_axis, charge_axis,
                       config.dtsig, config.mzsig, 0,
                       config.psfun, config.zpsfun);
    fftwf_execute(context->forward_plan);
    memcpy(kernel_fft, context->spectrum_work,
           (size_t)context->spectrum_length * sizeof(fftwf_complex));
    if (kernel_fft[0][0] == 0 || !isfinite(kernel_fft[0][0])) {
        fprintf(stderr, "Invalid UniChrom peak-shape kernel\n");
        return 0;
    }
    return 1;
}


static void normalize_kernel_fft_UniChrom(fftwf_complex *kernel_fft,
                                          const int length)
{
    const float inverse_sum = 1.0f / kernel_fft[0][0];
    #pragma omp parallel for schedule(static) if(length >= UNICHROM_OMP_MIN_LENGTH)
    for (int i = 0; i < length; i++) {
        kernel_fft[i][0] *= inverse_sum;
        kernel_fft[i][1] *= inverse_sum;
    }
}


static void make_allowed_grid_UniChrom(const Config config,
                                       const float *mz_axis, const float *mz_sum,
                                       const int mz_count, char *allowed)
{
    for (int mz = 0; mz < mz_count; mz++) {
        for (int charge_index = 0; charge_index < config.numz; charge_index++) {
            const int charge = config.startz + charge_index;
            const float mass = calcmass(mz_axis[mz], charge, config.adductmass);
            const float native_limit = nativecharge(mass, 0);
            allowed[index2D(config.numz, mz, charge_index)] =
                mz_sum[mz] > config.intthresh &&
                mass > config.masslb && mass < config.massub &&
                (float)charge < native_limit + config.nativezub &&
                (float)charge > native_limit + config.nativezlb;
        }
    }
}


static void clear_disallowed_UniChrom(float *cube, const char *allowed,
                                      const int scan_count, const int scan_length)
{
    #pragma omp parallel for schedule(static) if((int64_t)scan_count * scan_length >= UNICHROM_OMP_MIN_LENGTH)
    for (int scan = 0; scan < scan_count; scan++) {
        float *scan_data = cube + (size_t)scan * scan_length;
        for (int i = 0; i < scan_length; i++) {
            if (!allowed[i] || !isfinite(scan_data[i]) || scan_data[i] < 0) {
                scan_data[i] = 0;
            }
        }
    }
}


static void scale_merged_sum_UniChrom(float *sum, const int sum_length,
                                      const float *grid, const int grid_length,
                                      const int datanorm)
{
    if (datanorm == 1) {
        norm1d(sum, sum_length);
        return;
    }
    const float grid_max = Max(grid, grid_length);
    const float sum_max = Max(sum, sum_length);
    if (sum_max > 0) { Normalize(sum_length, sum, sum_max / grid_max); }
}


static void find_mass_axis_UniChrom(const Config config, const float *cube,
                                    const char *allowed, const float *mz_axis,
                                    const int size[3], float *mass_min,
                                    float *mass_max)
{
    if (config.fixedmassaxis) {
        *mass_min = config.masslb;
        *mass_max = config.massub;
        return;
    }
    const float cube_max = Max(cube, size[2] * size[0] * size[1]);
    const float cutoff = cube_max * 0.000001f;
    *mass_min = config.massub;
    *mass_max = config.masslb;
    for (int mz = 0; mz < size[0]; mz++) {
        for (int charge_index = 0; charge_index < size[1]; charge_index++) {
            const int grid_index = index2D(size[1], mz, charge_index);
            if (!allowed[grid_index]) { continue; }
            float maximum = 0;
            for (int scan = 0; scan < size[2]; scan++) {
                const float value = cube[(size_t)scan * size[0] * size[1] + grid_index];
                if (value > maximum) { maximum = value; }
            }
            if (maximum > cutoff) {
                const int charge = config.startz + charge_index;
                const float mass = calcmass(mz_axis[mz], charge, config.adductmass);
                const float low = mass - config.psmzthresh * (float)charge;
                const float high = mass + config.psmzthresh * (float)charge + config.massbins;
                if (low < *mass_min) { *mass_min = low; }
                if (high > *mass_max) { *mass_max = high; }
            }
        }
    }
    if (*mass_max <= *mass_min) {
        *mass_min = config.masslb;
        *mass_max = config.massub;
    }
    *mass_min = floorf(*mass_min / config.massbins) * config.massbins;
    *mass_max = ceilf(*mass_max / config.massbins) * config.massbins;
}


static void transform_mass_grid_UniChrom(const Config config, const float *cube,
                                         const char *allowed, const float *mz_axis,
                                         const int size[3], const float mass_min,
                                         const int mass_count, float *mass_grid)
{
    const int scan_length = size[0] * size[1];
    #pragma omp parallel for schedule(static) if(size[2] > 1)
    for (int scan = 0; scan < size[2]; scan++) {
        const float *scan_cube = cube + (size_t)scan * scan_length;
        float *scan_mass = mass_grid + (size_t)scan * mass_count;
        for (int mz = 0; mz < size[0]; mz++) {
            for (int charge_index = 0; charge_index < size[1]; charge_index++) {
                const int grid_index = index2D(size[1], mz, charge_index);
                const float value = scan_cube[grid_index];
                if (!allowed[grid_index] || value <= 0) { continue; }
                const int charge = config.startz + charge_index;
                const float mass = calcmass(mz_axis[mz], charge, config.adductmass);
                const float position = (mass - mass_min) / config.massbins;
                const int lower = (int)floorf(position);
                const float fraction = position - (float)lower;
                if (lower >= 0 && lower < mass_count) {
                    scan_mass[lower] += value * (1.0f - fraction);
                }
                if (lower + 1 >= 0 && lower + 1 < mass_count) {
                    scan_mass[lower + 1] += value * fraction;
                }
            }
        }
    }
}


static void write_outputs_UniChrom(const Config config, const float *cube,
                                   const char *allowed, const float *mz_axis,
                                   const int size[3])
{
    const int scan_count = size[2];
    const int mz_count = size[0];
    const int charge_count = size[1];
    const int scan_length = mz_count * charge_count;
    float *mz_grid = calloc((size_t)scan_count * mz_count, sizeof(float));
    float *mz_sum = calloc((size_t)mz_count, sizeof(float));
    if (mz_grid == NULL || mz_sum == NULL) {
        fprintf(stderr, "Unable to allocate UniChrom m/z outputs\n");
        exit(11);
    }
    #pragma omp parallel for schedule(static) if(scan_count > 1)
    for (int scan = 0; scan < scan_count; scan++) {
        const float *scan_cube = cube + (size_t)scan * scan_length;
        for (int mz = 0; mz < mz_count; mz++) {
            float value = 0;
            for (int charge = 0; charge < charge_count; charge++) {
                const int index = index2D(charge_count, mz, charge);
                if (allowed[index]) { value += scan_cube[index]; }
            }
            mz_grid[(size_t)scan * mz_count + mz] = value;
        }
    }
    for (int scan = 0; scan < scan_count; scan++) {
        for (int mz = 0; mz < mz_count; mz++) {
            mz_sum[mz] += mz_grid[(size_t)scan * mz_count + mz];
        }
    }
    scale_merged_sum_UniChrom(mz_sum, mz_count, mz_grid,
                              scan_count * mz_count, config.datanorm);

    float mass_min = 0;
    float mass_max = 0;
    find_mass_axis_UniChrom(config, cube, allowed, mz_axis, size,
                            &mass_min, &mass_max);
    const int mass_count = 1 + (int)((mass_max - mass_min) / config.massbins);
    float *mass_axis = calloc((size_t)mass_count, sizeof(float));
    float *mass_grid = calloc((size_t)scan_count * mass_count, sizeof(float));
    float *mass_sum = calloc((size_t)mass_count, sizeof(float));
    if (mass_axis == NULL || mass_grid == NULL || mass_sum == NULL) {
        fprintf(stderr, "Unable to allocate UniChrom mass outputs (%d bins)\n", mass_count);
        exit(11);
    }
    for (int mass = 0; mass < mass_count; mass++) {
        mass_axis[mass] = mass_min + (float)mass * config.massbins;
    }
    transform_mass_grid_UniChrom(config, cube, allowed, mz_axis, size,
                                 mass_min, mass_count, mass_grid);
    for (int scan = 0; scan < scan_count; scan++) {
        for (int mass = 0; mass < mass_count; mass++) {
            mass_sum[mass] += mass_grid[(size_t)scan * mass_count + mass];
        }
    }
    scale_merged_sum_UniChrom(mass_sum, mass_count, mass_grid,
                              scan_count * mass_count, config.datanorm);

    mh5writefile1d(config.file_id, "/ms_dataset/mz_grid",
                   scan_count * mz_count, mz_grid);
    mh5writefile1d(config.file_id, "/ms_dataset/mz_axis", mz_count, mz_axis);
    mh5writefile1d(config.file_id, "/ms_dataset/mz_sum", mz_count, mz_sum);
    mh5writefile1d(config.file_id, "/ms_dataset/mass_grid",
                   scan_count * mass_count, mass_grid);
    mh5writefile1d(config.file_id, "/ms_dataset/mass_axis", mass_count, mass_axis);
    mh5writefile1d(config.file_id, "/ms_dataset/mass_sum", mass_count, mass_sum);

    for (int scan = 0; scan < scan_count; scan++) {
        char group[1024];
        char path[1024];
        snprintf(group, sizeof(group), "/ms_dataset/%d", scan);
        snprintf(path, sizeof(path), "%s/mass_data", group);
        mh5writefile2d(config.file_id, path, mass_count, mass_axis,
                       mass_grid + (size_t)scan * mass_count);
        write_attr_int(config.file_id, group, "length_mz", mz_count);
        write_attr_int(config.file_id, group, "length_mass", mass_count);
        write_attr_float(config.file_id, group, "mzsig", config.mzsig);
        write_attr_float(config.file_id, group, "dtsig", config.dtsig);
    }
    set_got_grids(config.file_id);
    printf("UniChrom outputs: %d x %d m/z grid; %d x %d mass grid\n",
           scan_count, mz_count, scan_count, mass_count);
    free(mz_grid); free(mz_sum);
    free(mass_axis); free(mass_grid); free(mass_sum);
}


int run_chromatogram(int argc, char *argv[], Config config)
{
    const clock_t starttime = clock();
    int result = 1;
    UniChromFFT fft_context;
    memset(&fft_context, 0, sizeof(fft_context));
    config.file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
    if (config.file_id < 0) {
        fprintf(stderr, "Unable to open UniChrom HDF5 file: %s\n", argv[1]);
        return 2;
    }
    const int scan_count = int_attr(config.file_id, "/ms_dataset", "num", 0);
    if (scan_count < 1 || config.numz < 1 || config.dtsig <= 0) {
        fprintf(stderr, "Invalid UniChrom scan, charge, or dtsig configuration\n");
        H5Fclose(config.file_id);
        return 2;
    }
    if (config.baselineflag != 0) {
        printf("UniChrom note: coupled baseline fitting is not yet implemented; "
               "use subbuff for preprocessing if needed.\n");
    }
    if (config.manualflag != 0 || config.mflag != 0) {
        printf("UniChrom note: manual assignments and mass-list constraints are not yet applied.\n");
    }

    /* Refresh every processed spectrum, then use the established MetaUniDec
     * merger to create one explicitly linear m/z axis. */
    printf("UniChrom: processing and linearizing %d spectra\n", scan_count);
    for (int scan = 0; scan < scan_count; scan++) {
        config.metamode = scan;
        config.silent = 1;
        process_data(argc, argv, config);
    }
    make_grid(argc, argv, config, "/processed_data", "/mz_grid", "/mz_axis", "/mz_sum");

    const int mz_count = mh5getfilelength(config.file_id, "/ms_dataset/mz_axis");
    const int observed_length = mh5getfilelength(config.file_id, "/ms_dataset/mz_grid");
    if (mz_count < 2 || observed_length != scan_count * mz_count) {
        fprintf(stderr, "Invalid linearized UniChrom m/z grid dimensions\n");
        goto cleanup_file_UniChrom;
    }
    const int size[3] = {mz_count, config.numz, scan_count};
    const int scan_length = mz_count * config.numz;
    const int cube_length = scan_length * scan_count;
    float *mz_axis = calloc((size_t)mz_count, sizeof(float));
    float *mz_sum = calloc((size_t)mz_count, sizeof(float));
    float *observed = fftwf_malloc((size_t)observed_length * sizeof(float));
    float *charge_axis = calloc((size_t)config.numz, sizeof(float));
    char *allowed = calloc((size_t)scan_length, sizeof(char));
    if (mz_axis == NULL || mz_sum == NULL || observed == NULL ||
        charge_axis == NULL || allowed == NULL) {
        fprintf(stderr, "Unable to allocate UniChrom input arrays\n");
        goto cleanup_inputs_UniChrom;
    }
    mh5readfile1d(config.file_id, "/ms_dataset/mz_axis", mz_axis);
    mh5readfile1d(config.file_id, "/ms_dataset/mz_sum", mz_sum);
    mh5readfile1d(config.file_id, "/ms_dataset/mz_grid", observed);
    for (int charge = 0; charge < config.numz; charge++) {
        charge_axis[charge] = (float)(config.startz + charge);
    }
    make_allowed_grid_UniChrom(config, mz_axis, mz_sum, mz_count, allowed);

    int allowed_count = 0;
    for (int i = 0; i < scan_length; i++) { allowed_count += allowed[i] != 0; }
    if (allowed_count == 0) {
        fprintf(stderr, "UniChrom mass/charge constraints exclude the full grid\n");
        goto cleanup_inputs_UniChrom;
    }
    printf("UniChrom dimensions: %d chromatography x %d m/z x %d charge (%d cells)\n",
           scan_count, mz_count, config.numz, cube_length);
    printf("UniChrom peak widths: chromatography %g scans, m/z %g (internal sigma units)\n",
           config.dtsig, config.mzsig);

    float *blur = fftwf_malloc((size_t)cube_length * sizeof(float));
    float *scratch = fftwf_malloc((size_t)cube_length * sizeof(float));
    float *oldblur = fftwf_malloc((size_t)cube_length * sizeof(float));
    float *predicted = calloc((size_t)observed_length, sizeof(float));
    int *zupind = calloc((size_t)scan_length, sizeof(int));
    int *zloind = calloc((size_t)scan_length, sizeof(int));
    int *mupind = calloc((size_t)scan_length, sizeof(int));
    int *mloind = calloc((size_t)scan_length, sizeof(int));
    int *nztab = calloc((size_t)config.numz, sizeof(int));
    float *smoothing_sums = calloc((size_t)config.numz, sizeof(float));
    if (blur == NULL || scratch == NULL || oldblur == NULL || predicted == NULL ||
        zupind == NULL || zloind == NULL || mupind == NULL || mloind == NULL ||
        nztab == NULL || smoothing_sums == NULL) {
        fprintf(stderr, "Unable to allocate UniChrom iteration arrays\n");
        goto cleanup_processing_UniChrom;
    }
    memset(blur, 0, (size_t)cube_length * sizeof(float));
    const float initial_divisor = (float)(config.numz + 2);
    for (int scan = 0; scan < scan_count; scan++) {
        for (int mz = 0; mz < mz_count; mz++) {
            const float value = observed[(size_t)scan * mz_count + mz] / initial_divisor;
            for (int charge = 0; charge < config.numz; charge++) {
                const int scan_index = index2D(config.numz, mz, charge);
                if (allowed[scan_index]) {
                    blur[(size_t)scan * scan_length + scan_index] = value;
                }
            }
        }
    }
    memcpy(oldblur, blur, (size_t)cube_length * sizeof(float));

    const float mzranges[4] = {mz_axis[0], mz_axis[mz_count - 1],
                               charge_axis[0], charge_axis[config.numz - 1]};
    float *mz_coordinates = calloc((size_t)scan_length, sizeof(float));
    float *z_coordinates = calloc((size_t)scan_length, sizeof(float));
    if (mz_coordinates == NULL || z_coordinates == NULL) {
        fprintf(stderr, "Unable to allocate UniChrom coordinate arrays\n");
        free(mz_coordinates); free(z_coordinates);
        goto cleanup_processing_UniChrom;
    }
    for (int mz = 0; mz < mz_count; mz++) {
        for (int charge = 0; charge < config.numz; charge++) {
            const int index = index2D(config.numz, mz, charge);
            mz_coordinates[index] = mz_axis[mz];
            z_coordinates[index] = charge_axis[charge];
            nztab[charge] = config.startz + charge;
        }
    }
    if (config.zsig != 0) {
        setup_blur_z_UCCD(zupind, zloind, mz_coordinates, z_coordinates,
                          scan_length, config.adductmass, mzranges, size);
    }
    if (config.msig != 0) {
        setup_blur_m_UCCD(mupind, mloind, mz_coordinates, z_coordinates,
                          scan_length, config.adductmass, mzranges, size, config.molig);
    }
    free(mz_coordinates); free(z_coordinates);

    if (!initialize_fft_UniChrom(&fft_context, size)) {
        goto cleanup_processing_UniChrom;
    }
    fftwf_complex *kernel_fft = fftwf_malloc(
        (size_t)fft_context.spectrum_length * sizeof(fftwf_complex));
    if (kernel_fft == NULL ||
        !make_kernel_fft_UniChrom(kernel_fft, &fft_context, size,
                                  mz_axis, charge_axis, config)) {
        fftwf_free(kernel_fft);
        goto cleanup_processing_UniChrom;
    }

    const float data_max = Max(observed, observed_length);
    const float beta_factor = data_max > 1 ? data_max : 1;
    double convergence = 0;
    int convergence_seen = 0;
    printf("UniChrom: iterating with coupled 3-D FFTs.");
    for (int iteration = 0; iteration < abs(config.numit); iteration++) {
        if (config.beta > 0 && iteration > 0) {
            #pragma omp parallel for schedule(static) if(scan_count > 1)
            for (int scan = 0; scan < scan_count; scan++) {
                softargmax(blur + (size_t)scan * scan_length,
                           mz_count, config.numz, config.beta / beta_factor);
            }
        }
        if (config.psig >= 1 && iteration > 0) {
            for (int scan = 0; scan < scan_count; scan++) {
                point_smoothing(blur + (size_t)scan * scan_length,
                                scratch + (size_t)scan * scan_length,
                                smoothing_sums, allowed, mz_count, config.numz,
                                abs((int)config.psig));
            }
        }
        if (iteration > config.suppression_startit &&
            (config.suppression_satellite > 0 || config.suppression_harmonic > 0 ||
             config.suppression_topn > 0 || config.suppression_topx > 0)) {
            apply_suppressions(blur, scratch, scan_count * mz_count, config.numz,
                               config.suppression_satellite, config.suppression_harmonic,
                               nztab, config.suppression_topn, config.suppression_topx,
                               config.suppression_percent);
        }
        if (config.zsig != 0 || config.msig != 0) {
            #pragma omp parallel for schedule(static) if(scan_count > 1)
            for (int scan = 0; scan < scan_count; scan++) {
                const size_t offset = (size_t)scan * scan_length;
                if (config.zsig != 0) {
                    blur_it_UCCD(scratch + offset, blur + offset, zupind, zloind,
                                 scan_length, config.zsig * data_max);
                    memcpy(blur + offset, scratch + offset,
                           (size_t)scan_length * sizeof(float));
                }
                if (config.msig != 0) {
                    blur_it_UCCD(scratch + offset, blur + offset, mupind, mloind,
                                 scan_length, config.msig * data_max);
                    memcpy(blur + offset, scratch + offset,
                           (size_t)scan_length * sizeof(float));
                }
            }
        }
        clear_disallowed_UniChrom(blur, allowed, scan_count, scan_length);

        execute_convolution_UniChrom(blur, kernel_fft, 0, &fft_context);
        const float fft_scale = 1.0f / (float)fft_context.normalization;
        #pragma omp parallel for schedule(static) if(observed_length >= UNICHROM_OMP_MIN_LENGTH)
        for (int index = 0; index < observed_length; index++) {
            const int scan = index / mz_count;
            const int mz = index - scan * mz_count;
            const float *row = fft_context.real_work +
                               (size_t)scan * scan_length + (size_t)mz * config.numz;
            float value = 0;
            for (int charge = 0; charge < config.numz; charge++) {
                value += row[charge] * fft_scale;
            }
            predicted[index] = value;
        }
        for (int scan = 0; scan < scan_count; scan++) {
            for (int mz = 0; mz < mz_count; mz++) {
                const int observed_index = scan * mz_count + mz;
                const float denominator = predicted[observed_index];
                const float ratio = denominator > 0 ? observed[observed_index] / denominator : 0;
                float *row = fft_context.real_work +
                             (size_t)scan * scan_length + (size_t)mz * config.numz;
                for (int charge = 0; charge < config.numz; charge++) {
                    row[charge] = ratio;
                }
            }
        }
        fftwf_execute(fft_context.forward_plan);
        #pragma omp parallel for schedule(static) if(fft_context.spectrum_length >= UNICHROM_OMP_MIN_LENGTH)
        for (int i = 0; i < fft_context.spectrum_length; i++) {
            const float a = fft_context.spectrum_work[i][0];
            const float b = fft_context.spectrum_work[i][1];
            const float c = kernel_fft[i][0];
            const float d = -kernel_fft[i][1];
            fft_context.spectrum_work[i][0] = a * c - b * d;
            fft_context.spectrum_work[i][1] = b * c + a * d;
        }
        fftwf_execute(fft_context.backward_plan);
        #pragma omp parallel for schedule(static) if(cube_length >= UNICHROM_OMP_MIN_LENGTH)
        for (int i = 0; i < cube_length; i++) {
            blur[i] *= fft_context.real_work[i] * fft_scale;
        }
        clear_disallowed_UniChrom(blur, allowed, scan_count, scan_length);

        if (abs(config.numit) < 10 || iteration == 1 || iteration % 10 == 0 ||
            iteration >= 9 * abs(config.numit) / 10) {
            double difference = 0;
            double total = 0;
            #pragma omp parallel for reduction(+:difference,total) schedule(static) if(cube_length >= UNICHROM_OMP_MIN_LENGTH)
            for (int i = 0; i < cube_length; i++) {
                const double delta = (double)blur[i] - oldblur[i];
                difference += delta * delta;
                total += blur[i];
            }
            convergence = total > 0 ? difference / total : INFINITY;
            if (convergence < 0.000001) {
                if (convergence_seen && config.numit > 0) {
                    printf(" converged in %d iterations", iteration + 1);
                    break;
                }
                convergence_seen = 1;
            } else {
                convergence_seen = 0;
            }
            memcpy(oldblur, blur, (size_t)cube_length * sizeof(float));
        }
    }
    printf(" done (metric %.6g)\n", convergence);

    const float *output_cube = blur;
    if (config.rawflag == 0 || config.rawflag == 2) {
        normalize_kernel_fft_UniChrom(kernel_fft, fft_context.spectrum_length);
        execute_convolution_UniChrom(blur, kernel_fft, 0, &fft_context);
        const float scale = 1.0f / (float)fft_context.normalization;
        #pragma omp parallel for schedule(static) if(cube_length >= UNICHROM_OMP_MIN_LENGTH)
        for (int i = 0; i < cube_length; i++) {
            scratch[i] = fft_context.real_work[i] * scale;
        }
        clear_disallowed_UniChrom(scratch, allowed, scan_count, scan_length);
        output_cube = scratch;
    }
    if (config.datanorm == 1) {
        const float output_max = Max(output_cube, cube_length);
        if (output_max > 0 && data_max > 0) {
            Normalize(cube_length, (float *)output_cube, output_max / data_max);
        }
    }
    ApplyCutoff((float *)output_cube, 0, cube_length);
    write_outputs_UniChrom(config, output_cube, allowed, mz_axis, size);
    result = 0;

    fftwf_free(kernel_fft);
cleanup_processing_UniChrom:
    destroy_fft_UniChrom(&fft_context);
    fftwf_free(blur); fftwf_free(scratch); fftwf_free(oldblur);
    free(predicted); free(zupind); free(zloind); free(mupind); free(mloind);
    free(nztab); free(smoothing_sums);
cleanup_inputs_UniChrom:
    free(mz_axis); free(mz_sum); fftwf_free(observed);
    free(charge_axis); free(allowed);
cleanup_file_UniChrom:
    H5Fclose(config.file_id);
    printf("UniChrom finished in %.3f s\n",
           (float)(clock() - starttime) / CLOCKS_PER_SEC);
    return result;
}
