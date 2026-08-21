/*
 * Charge-detection deconvolution with a chromatography dimension.
 *
 * Storage order is [chromatography][m/z][charge], with charge contiguous.
 */

#include "UCCD_Main.h"


static int readfile4_UCCD(const char *infile, const int length, float *chromdat,
                          float *mzdat, float *zdat, float *intensity)
{
    FILE *file_ptr = fopen(infile, "r");
    if (file_ptr == NULL) {
        fprintf(stderr, "Error opening %s\n", infile);
        return 0;
    }

    for (int i = 0; i < length; i++) {
        if (fscanf(file_ptr, "%f %f %f %f", &chromdat[i], &mzdat[i],
                   &zdat[i], &intensity[i]) != 4) {
            fprintf(stderr, "Error reading row %d from %s; expected four numeric columns\n",
                    i + 1, infile);
            fclose(file_ptr);
            return 0;
        }
    }
    fclose(file_ptr);
    return 1;
}


static float periodic_axis_peak_UCCD(const float *axis, const int length,
                                     const int index, const float sig,
                                     const int psfun)
{
    if (sig == 0) {
        return index == 0 ? 1.0f : 0.0f;
    }
    if (length == 1) {
        return 1.0f;
    }

    const float low = axis[0];
    const float high = 2.0f * axis[length - 1] - axis[length - 2];
    return mzpeakshape(low, axis[index], sig, psfun)
           + mzpeakshape(high, axis[index], sig, psfun);
}


void blur_it_UCCD(float *output, const float *input, const int *upinds,
                  const int *loinds, const int length, const float floor)
{
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < length; i++) {
        float i1 = input[i];
        float i2 = input[loinds[i]];
        float i3 = input[upinds[i]];
        if (floor > 0) {
            i1 = logf(i1 + floor);
            i2 = logf(i2 + floor);
            i3 = logf(i3 + floor);
            if (isnan(i1) || isinf(i1)) { i1 = 0; }
            if (isnan(i2) || isinf(i2)) { i2 = 0; }
            if (isnan(i3) || isinf(i3)) { i3 = 0; }

            float newval = expf((i1 + i2 + i3) / 3.0f) - floor;
            output[i] = newval > 0 ? newval : 0;
        } else {
            const float ratio = fabsf(floor);
            output[i] = (i1 + i2 * ratio + i3 * ratio) / 3.0f;
        }
    }
}


void setup_blur_z_UCCD(int *zupind, int *zloind, const float *mzdat,
                       const float *zdat, const int scan_length,
                       const float adductmass, const float mzranges[4],
                       const int size[3])
{
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < scan_length; i++) {
        const float mz = mzdat[i];
        const float z = zdat[i];
        const int zind = (int)(z - mzranges[2]);
        const int zint = (int)z;
        const float mass = calcmass(mz, zint, adductmass);
        const float uppermz = calcmz(mass, zint + 1, adductmass);
        const float lowermz = zint != 1 ? calcmz(mass, zint - 1, adductmass) : 0;

        const int mzupind = (int)roundf((uppermz - mzranges[0]) /
                            (mzranges[1] - mzranges[0]) * (float)(size[0] - 1));
        if (z == mzranges[3] || mzupind < 0 || mzupind >= size[0]) {
            zupind[i] = i;
        } else {
            int newind = index2D(size[1], mzupind, zind + 1);
            if (newind < 0) { newind = 0; }
            if (newind >= scan_length) { newind = scan_length - 1; }
            zupind[i] = newind;
        }

        const int mzloind = (int)roundf((lowermz - mzranges[0]) /
                            (mzranges[1] - mzranges[0]) * (float)(size[0] - 1));
        if (zind == 0 || mzloind < 0 || mzloind >= size[0]) {
            zloind[i] = i;
        } else {
            int newind = index2D(size[1], mzloind, zind - 1);
            if (newind < 0) { newind = 0; }
            if (newind >= scan_length) { newind = scan_length - 1; }
            zloind[i] = newind;
        }
    }
}


void setup_blur_m_UCCD(int *mupind, int *mloind, const float *mzdat,
                       const float *zdat, const int scan_length,
                       const float adductmass, const float mzranges[4],
                       const int size[3], const float molig)
{
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < scan_length; i++) {
        const float mz = mzdat[i];
        const int zint = (int)zdat[i];
        const int zind = (int)(zdat[i] - mzranges[2]);
        const float mass = calcmass(mz, zint, adductmass);
        const float uppermz = calcmz(mass + molig, zint, adductmass);
        const float lowermz = calcmz(mass - molig, zint, adductmass);

        const int mzupind = (int)roundf((uppermz - mzranges[0]) /
                            (mzranges[1] - mzranges[0]) * (float)(size[0] - 1));
        if (mzupind < 0 || mzupind >= size[0]) {
            mupind[i] = i;
        } else {
            mupind[i] = index2D(size[1], mzupind, zind);
        }

        const int mzloind = (int)roundf((lowermz - mzranges[0]) /
                            (mzranges[1] - mzranges[0]) * (float)(size[0] - 1));
        if (mzloind < 0 || mzloind >= size[0]) {
            mloind[i] = i;
        } else {
            mloind[i] = index2D(size[1], mzloind, zind);
        }
    }
}


void make_kernel3D_UCCD(float *peak, const int size[3], const float *chromext,
                        const float *mzext, const float *zext,
                        const float chromsig, const float mzsig,
                        const float zsig, const int psfun, const int zpsfun)
{
    const int scan_length = size[0] * size[1];
    float *scan_kernel = calloc((size_t)scan_length, sizeof(float));
    if (scan_kernel == NULL) {
        fprintf(stderr, "Error allocating UCCD scan kernel\n");
        exit(1);
    }

    MakeKernel2D(scan_kernel, size, mzext, zext, mzsig, zsig, psfun, zpsfun);
    #pragma omp parallel for schedule(dynamic)
    for (int scan = 0; scan < size[2]; scan++) {
        const float chrom_value = periodic_axis_peak_UCCD(chromext, size[2], scan,
                                                          chromsig, psfun);
        const int offset = scan * scan_length;
        for (int i = 0; i < scan_length; i++) {
            peak[offset + i] = chrom_value * scan_kernel[i];
        }
    }
    free(scan_kernel);
}


void precompute_fft3D_UCCD(const float *list, const int size[3],
                           fftwf_complex *output)
{
    const int length = size[0] * size[1] * size[2];
    fftwf_complex *input = fftwf_malloc(sizeof(fftwf_complex) * (size_t)length);
    if (input == NULL) {
        fprintf(stderr, "Error allocating UCCD FFT input\n");
        exit(1);
    }
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < length; i++) {
        input[i][0] = list[i];
        input[i][1] = 0;
    }
    fftwf_plan plan = fftwf_plan_dft_3d(size[2], size[0], size[1], input,
                                        output, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_execute(plan);
    fftwf_destroy_plan(plan);
    fftwf_free(input);
}


void fftconvolve3D_precomputed_UCCD(float *corr, const float *list,
                                    const fftwf_complex *kernel_fft,
                                    const int size[3], fftwf_plan forward_plan,
                                    fftwf_plan backward_plan,
                                    fftwf_complex *work,
                                    fftwf_complex *transform)
{
    const int length = size[0] * size[1] * size[2];
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < length; i++) {
        work[i][0] = list[i];
        work[i][1] = 0;
    }
    fftwf_execute(forward_plan);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < length; i++) {
        const float a = transform[i][0];
        const float b = transform[i][1];
        const float c = kernel_fft[i][0];
        const float d = kernel_fft[i][1];
        work[i][0] = a * c - b * d;
        work[i][1] = b * c + a * d;
    }
    fftwf_execute(backward_plan);

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < length; i++) {
        corr[i] = transform[i][0] / (float)length;
    }
}


void complex_conjugate_UCCD(const fftwf_complex *input, fftwf_complex *output,
                            const int length)
{
    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < length; i++) {
        output[i][0] = input[i][0];
        output[i][1] = -input[i][1];
    }
}


int run_unidec_UCCD(int argc, char *argv[], Config config)
{
    (void)argc;
    (void)argv;
    const time_t starttime = time(NULL);
    printf("Opening UCCD file: %s\n", config.infile);

    const int lines = getfilelength(config.infile);
    if (lines <= 0) {
        fprintf(stderr, "UCCD input contains no data\n");
        return 1;
    }

    float *chromdat = calloc((size_t)lines, sizeof(float));
    float *mzdat = calloc((size_t)lines, sizeof(float));
    float *zdat = calloc((size_t)lines, sizeof(float));
    float *dataInt = calloc((size_t)lines, sizeof(float));
    if (chromdat == NULL || mzdat == NULL || zdat == NULL || dataInt == NULL) {
        fprintf(stderr, "Error allocating UCCD data arrays\n");
        return 1;
    }
    if (!readfile4_UCCD(config.infile, lines, chromdat, mzdat, zdat, dataInt)) {
        return 2;
    }

    int size[3] = {GetSize0(mzdat, lines), GetSize1(zdat, lines), 0};
    const int scan_length = size[0] * size[1];
    if (scan_length <= 0 || lines % scan_length != 0) {
        fprintf(stderr, "UCCD input is not a complete dense m/z-by-charge grid per scan\n");
        return 2;
    }
    size[2] = lines / scan_length;
    for (int scan = 0; scan < size[2]; scan++) {
        const int offset = scan * scan_length;
        if (chromdat[offset] != chromdat[offset + scan_length - 1]) {
            fprintf(stderr, "Chromatography coordinate changes inside scan %d\n", scan);
            return 2;
        }
        for (int i = 0; i < scan_length; i++) {
            if (mzdat[offset + i] != mzdat[i] || zdat[offset + i] != zdat[i]) {
                fprintf(stderr, "m/z-charge grid differs in chromatography scan %d\n", scan);
                return 2;
            }
        }
    }
    printf("Dimensions: %d chromatography by %d m/z by %d charge: %d total\n",
           size[2], size[0], size[1], lines);

    float *chromext = calloc((size_t)size[2], sizeof(float));
    float *mzext = calloc((size_t)size[0], sizeof(float));
    float *zext = calloc((size_t)size[1], sizeof(float));
    if (chromext == NULL || mzext == NULL || zext == NULL) {
        fprintf(stderr, "Error allocating UCCD axis arrays\n");
        return 1;
    }
    for (int scan = 0; scan < size[2]; scan++) {
        chromext[scan] = chromdat[scan * scan_length];
    }
    PullXY(mzext, zext, mzdat, zdat, size);
    const float mzranges[4] = {mzext[0], mzext[size[0] - 1],
                               zext[0], zext[size[1] - 1]};
    printf("Chromatography range: %f to %f\n", chromext[0], chromext[size[2] - 1]);
    printf("MZ range: %f to %f; charge range: %f to %f\n",
           mzranges[0], mzranges[1], mzranges[2], mzranges[3]);

    const float dmax = Max(dataInt, lines);
    const float betafactor = dmax > 1 ? dmax : 1;
    int *zupind = calloc((size_t)scan_length, sizeof(int));
    int *zloind = calloc((size_t)scan_length, sizeof(int));
    int *mupind = calloc((size_t)scan_length, sizeof(int));
    int *mloind = calloc((size_t)scan_length, sizeof(int));
    char *barr = calloc((size_t)scan_length, sizeof(char));
    if (zupind == NULL || zloind == NULL || mupind == NULL || mloind == NULL ||
        barr == NULL) {
        fprintf(stderr, "Error allocating UCCD scan-processing arrays\n");
        return 1;
    }
    for (int i = 0; i < scan_length; i++) { barr[i] = 1; }
    if (config.zsig != 0) {
        setup_blur_z_UCCD(zupind, zloind, mzdat, zdat, scan_length,
                          config.adductmass, mzranges, size);
    }
    if (config.msig != 0) {
        setup_blur_m_UCCD(mupind, mloind, mzdat, zdat, scan_length,
                          config.adductmass, mzranges, size, config.molig);
    }

    float *peakshape = calloc((size_t)lines, sizeof(float));
    float *mkernel = calloc((size_t)lines, sizeof(float));
    float *blur = calloc((size_t)lines, sizeof(float));
    float *newblur = calloc((size_t)lines, sizeof(float));
    float *newblur2 = calloc((size_t)lines, sizeof(float));
    float *oldblur = calloc((size_t)lines, sizeof(float));
    fftwf_complex *work = fftwf_malloc(sizeof(fftwf_complex) * (size_t)lines);
    fftwf_complex *transform = fftwf_malloc(sizeof(fftwf_complex) * (size_t)lines);
    fftwf_complex *peakshape_fft = fftwf_malloc(sizeof(fftwf_complex) * (size_t)lines);
    fftwf_complex *inverse_peakshape_fft = fftwf_malloc(sizeof(fftwf_complex) * (size_t)lines);
    fftwf_complex *mkernel_fft = fftwf_malloc(sizeof(fftwf_complex) * (size_t)lines);
    if (peakshape == NULL || mkernel == NULL || blur == NULL || newblur == NULL ||
        newblur2 == NULL || oldblur == NULL || work == NULL || transform == NULL ||
        peakshape_fft == NULL || inverse_peakshape_fft == NULL || mkernel_fft == NULL) {
        fprintf(stderr, "Error allocating UCCD processing arrays\n");
        return 1;
    }

    fftwf_plan forward_plan = fftwf_plan_dft_3d(size[2], size[0], size[1], work,
                                                transform, FFTW_FORWARD, FFTW_ESTIMATE);
    fftwf_plan backward_plan = fftwf_plan_dft_3d(size[2], size[0], size[1], work,
                                                 transform, FFTW_BACKWARD, FFTW_ESTIMATE);
    make_kernel3D_UCCD(peakshape, size, chromext, mzext, zext,
                       config.csig, config.mzsig, 0,
                       config.psfun, config.zpsfun);
    precompute_fft3D_UCCD(peakshape, size, peakshape_fft);
    complex_conjugate_UCCD(peakshape_fft, inverse_peakshape_fft, lines);

    const size_t matsize = (size_t)lines * sizeof(float);
    memcpy(blur, dataInt, matsize);
    memcpy(oldblur, blur, matsize);
    printf("Iterating.");
    float conv = 0;
    int off = 0;
    for (int iteration = 0; iteration < config.numit; iteration++) {
        /* These operations intentionally do not cross chromatography scans. */
        for (int scan = 0; scan < size[2]; scan++) {
            const int offset = scan * scan_length;
            if (config.beta > 0) {
                softargmax(blur + offset, size[0], size[1], config.beta / betafactor);
            }
            if (config.psig > 0) {
                point_smoothing(blur + offset, barr, size[0], size[1],
                                abs((int)config.psig));
            }
            if (config.zsig != 0) {
                blur_it_UCCD(newblur + offset, blur + offset, zupind, zloind,
                             scan_length, config.zsig * dmax);
                memcpy(blur + offset, newblur + offset,
                       (size_t)scan_length * sizeof(float));
            }
            if (config.msig != 0) {
                blur_it_UCCD(newblur + offset, blur + offset, mupind, mloind,
                             scan_length, config.msig * dmax);
                memcpy(blur + offset, newblur + offset,
                       (size_t)scan_length * sizeof(float));
            }
        }

        fftconvolve3D_precomputed_UCCD(newblur, blur, peakshape_fft, size,
                                      forward_plan, backward_plan, work, transform);
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < lines; i++) {
            newblur2[i] = newblur[i] != 0 ? dataInt[i] / newblur[i] : 0;
        }
        fftconvolve3D_precomputed_UCCD(newblur, newblur2, inverse_peakshape_fft,
                                      size, forward_plan, backward_plan, work, transform);
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < lines; i++) { blur[i] *= newblur[i]; }

        if (config.numit < 10 || iteration % 10 == 0 || iteration % 10 == 1 ||
            iteration > 0.9 * config.numit) {
            float diff = 0;
            float total = 0;
            for (int i = 0; i < lines; i++) {
                diff += powf(blur[i] - oldblur[i], 2);
                total += blur[i];
            }
            if (total != 0) {
                conv = diff / total;
            } else if (conv == 12345678) {
                printf("UCCD grid is zero. Iteration: %d\n", iteration);
                break;
            } else {
                conv = 12345678;
            }
            if (conv < 0.000001) {
                if (off == 1 && config.numit > 0) {
                    printf("Converged in %d iterations.\n", iteration);
                    break;
                }
                off = 1;
            }
            memcpy(oldblur, blur, matsize);
        }
    }
    printf("Completed iterations\n");

    memcpy(newblur, blur, matsize);
    fftconvolve3D_precomputed_UCCD(newblur2, newblur, peakshape_fft, size,
                                  forward_plan, backward_plan, work, transform);
    if (config.datanorm == 1) {
        const float fitmax = Max(newblur2, lines);
        if (dmax != 0 && fitmax != 0) { Normalize(lines, newblur2, fitmax / dmax); }
    }
    ApplyCutoff(newblur2, 0, lines);
    write1D(config.outfile, "fitdat", newblur2, lines);

    if (config.rawflag == 0) {
        make_kernel3D_UCCD(mkernel, size, chromext, mzext, zext, 0,
                           config.mzsig, 0, config.psfun, config.zpsfun);
        precompute_fft3D_UCCD(mkernel, size, mkernel_fft);
        memcpy(newblur, blur, matsize);
        fftconvolve3D_precomputed_UCCD(blur, newblur, mkernel_fft, size,
                                      forward_plan, backward_plan, work, transform);
        printf("Reconvolved with m/z dimension\n");
    }
    if (config.datanorm == 1) {
        const float blurmax = Max(blur, lines);
        if (dmax != 0 && blurmax != 0) { Normalize(lines, blur, blurmax / dmax); }
    }
    ApplyCutoff(blur, 0, lines);

    char outstring[510];
    snprintf(outstring, sizeof(outstring), "%s_decon.txt", config.outfile);
    FILE *out_ptr = fopen(outstring, "w");
    if (out_ptr == NULL) {
        fprintf(stderr, "Error opening %s\n", outstring);
        return 1;
    }
    for (int i = 0; i < lines; i++) { fprintf(out_ptr, "%f\n", blur[i]); }
    fclose(out_ptr);
    printf("Wrote deconvolution output to: %s\n", outstring);

    free(chromdat); free(mzdat); free(zdat); free(dataInt);
    free(chromext); free(mzext); free(zext);
    free(zupind); free(zloind); free(mupind); free(mloind); free(barr);
    free(peakshape); free(mkernel); free(blur); free(newblur); free(newblur2); free(oldblur);
    fftwf_destroy_plan(forward_plan); fftwf_destroy_plan(backward_plan);
    fftwf_free(work); fftwf_free(transform); fftwf_free(peakshape_fft);
    fftwf_free(inverse_peakshape_fft); fftwf_free(mkernel_fft);

    printf("Done in %ds!\n", (int)difftime(time(NULL), starttime));
    return 0;
}
