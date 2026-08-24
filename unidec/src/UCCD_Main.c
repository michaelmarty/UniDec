/*
 * Charge-detection deconvolution with a chromatography dimension.
 *
 * Storage order is [chromatography][m/z][charge], with charge contiguous.
 */

#include "UCCD_Main.h"

#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef UNIDEC_USE_MKL
#include <mkl_service.h>
#endif

#define UCCD_MAGIC "UCCDBIN1"
#define UCCD_INPUT_SUFFIX "_uccd_input.bin"
#define UCCD_DECON_SUFFIX "_uccd_decon.bin"
/* Fit output is temporarily disabled. */
/* #define UCCD_FIT_SUFFIX "_uccd_fit.bin" */

#define UCCD_OMP_MIN_LENGTH 32768
#define UCCD_FFT_PLAN_FLAGS FFTW_MEASURE

typedef struct {
    uint32_t index;
    float intensity;
} UCCDRecord;

typedef struct {
    int real_length;
    int spectrum_length;
    int kernel_length;
    int normalization;
    int batched_2d;
    float *real_work;
    fftwf_complex *spectrum_work;
    fftwf_plan forward_plan;
    fftwf_plan backward_plan;
} UCCDFFTContext;

_Static_assert(sizeof(UCCDRecord) == 8, "Unexpected UCCD sparse record padding");


static int make_filename_UCCD(char *filename, const size_t length,
                              const char *outfile, const char *suffix)
{
    const int written = snprintf(filename, length, "%s%s", outfile, suffix);
    if (written < 0 || (size_t)written >= length) {
        fprintf(stderr, "UCCD file name is too long\n");
        return 0;
    }
    return 1;
}


static int read_sparse_UCCD(const char *filename, int size[3],
                            float **chromext, float **mzext, float **zext,
                            float **intensity)
{
    FILE *file_ptr = fopen(filename, "rb");
    if (file_ptr == NULL) {
        fprintf(stderr, "Error opening %s\n", filename);
        return 0;
    }

    char magic[8];
    uint32_t dimensions[3];
    uint64_t nonzero_count;
    if (fread(magic, 1, sizeof(magic), file_ptr) != sizeof(magic) ||
        memcmp(magic, UCCD_MAGIC, sizeof(magic)) != 0 ||
        fread(dimensions, sizeof(uint32_t), 3, file_ptr) != 3 ||
        fread(&nonzero_count, sizeof(uint64_t), 1, file_ptr) != 1) {
        fprintf(stderr, "Invalid UCCD binary header in %s\n", filename);
        fclose(file_ptr);
        return 0;
    }

    if (dimensions[0] == 0 || dimensions[1] == 0 || dimensions[2] == 0 ||
        dimensions[0] > INT_MAX || dimensions[1] > INT_MAX || dimensions[2] > INT_MAX) {
        fprintf(stderr, "Invalid UCCD dimensions in %s\n", filename);
        fclose(file_ptr);
        return 0;
    }
    size[2] = (int)dimensions[0];
    size[0] = (int)dimensions[1];
    size[1] = (int)dimensions[2];
    const uint64_t total = (uint64_t)size[2] * (uint64_t)size[0] * (uint64_t)size[1];
    if (total > INT_MAX || total > UINT32_MAX || nonzero_count > total) {
        fprintf(stderr, "UCCD cube is too large or has an invalid sparse count\n");
        fclose(file_ptr);
        return 0;
    }

    *chromext = calloc((size_t)size[2], sizeof(float));
    *mzext = calloc((size_t)size[0], sizeof(float));
    *zext = calloc((size_t)size[1], sizeof(float));
    *intensity = fftwf_malloc((size_t)total * sizeof(float));
    if (*chromext == NULL || *mzext == NULL || *zext == NULL || *intensity == NULL) {
        fprintf(stderr, "Error allocating UCCD input arrays\n");
        fclose(file_ptr);
        return 0;
    }
    memset(*intensity, 0, (size_t)total * sizeof(float));

    if (fread(*chromext, sizeof(float), (size_t)size[2], file_ptr) != (size_t)size[2] ||
        fread(*mzext, sizeof(float), (size_t)size[0], file_ptr) != (size_t)size[0] ||
        fread(*zext, sizeof(float), (size_t)size[1], file_ptr) != (size_t)size[1]) {
        fprintf(stderr, "Incomplete UCCD axes in %s\n", filename);
        fclose(file_ptr);
        return 0;
    }

    const size_t buffer_length = 65536;
    UCCDRecord *records = malloc(buffer_length * sizeof(UCCDRecord));
    if (records == NULL) {
        fprintf(stderr, "Error allocating UCCD sparse read buffer\n");
        fclose(file_ptr);
        return 0;
    }
    uint64_t records_left = nonzero_count;
    while (records_left > 0) {
        const size_t count = records_left < buffer_length ? (size_t)records_left : buffer_length;
        if (fread(records, sizeof(UCCDRecord), count, file_ptr) != count) {
            fprintf(stderr, "Incomplete UCCD sparse data in %s\n", filename);
            free(records);
            fclose(file_ptr);
            return 0;
        }
        for (size_t i = 0; i < count; i++) {
            if (records[i].index >= total) {
                fprintf(stderr, "Invalid UCCD sparse index in %s\n", filename);
                free(records);
                fclose(file_ptr);
                return 0;
            }
            (*intensity)[records[i].index] = records[i].intensity;
        }
        records_left -= count;
    }
    free(records);
    fclose(file_ptr);
    return 1;
}


static int write_sparse_UCCD(const char *filename, const int size[3],
                             const float *chromext, const float *mzext,
                             const float *zext, const float *intensity)
{
    const uint64_t total = (uint64_t)size[2] * (uint64_t)size[0] * (uint64_t)size[1];
    uint64_t nonzero_count = 0;
    for (uint64_t i = 0; i < total; i++) {
        if (intensity[i] != 0) { nonzero_count++; }
    }

    FILE *file_ptr = fopen(filename, "wb");
    if (file_ptr == NULL) {
        fprintf(stderr, "Error opening %s\n", filename);
        return 0;
    }
    const uint32_t dimensions[3] = {
        (uint32_t)size[2], (uint32_t)size[0], (uint32_t)size[1]
    };
    if (fwrite(UCCD_MAGIC, 1, 8, file_ptr) != 8 ||
        fwrite(dimensions, sizeof(uint32_t), 3, file_ptr) != 3 ||
        fwrite(&nonzero_count, sizeof(uint64_t), 1, file_ptr) != 1 ||
        fwrite(chromext, sizeof(float), (size_t)size[2], file_ptr) != (size_t)size[2] ||
        fwrite(mzext, sizeof(float), (size_t)size[0], file_ptr) != (size_t)size[0] ||
        fwrite(zext, sizeof(float), (size_t)size[1], file_ptr) != (size_t)size[1]) {
        fprintf(stderr, "Error writing UCCD header or axes to %s\n", filename);
        fclose(file_ptr);
        return 0;
    }
    const size_t buffer_length = 65536;
    UCCDRecord *records = malloc(buffer_length * sizeof(UCCDRecord));
    if (records == NULL) {
        fprintf(stderr, "Error allocating UCCD sparse write buffer\n");
        fclose(file_ptr);
        return 0;
    }
    size_t record_count = 0;
    for (uint64_t i = 0; i < total; i++) {
        if (intensity[i] != 0) {
            records[record_count].index = (uint32_t)i;
            records[record_count].intensity = intensity[i];
            record_count++;
            if (record_count == buffer_length) {
                if (fwrite(records, sizeof(UCCDRecord), record_count, file_ptr) != record_count) {
                    fprintf(stderr, "Error writing UCCD sparse data to %s\n", filename);
                    free(records);
                    fclose(file_ptr);
                    return 0;
                }
                record_count = 0;
            }
        }
    }
    if (record_count > 0 &&
        fwrite(records, sizeof(UCCDRecord), record_count, file_ptr) != record_count) {
        fprintf(stderr, "Error writing UCCD sparse data to %s\n", filename);
        free(records);
        fclose(file_ptr);
        return 0;
    }
    free(records);
    fclose(file_ptr);
    printf("Wrote %llu nonzero UCCD values to: %s\n",
           (unsigned long long)nonzero_count, filename);
    return 1;
}


static int collect_nonzero_indices_UCCD(const float *intensity, const int length,
                                        uint32_t **nonzero_indices)
{
    int nonzero_count = 0;
    for (int i = 0; i < length; i++) {
        if (intensity[i] != 0) { nonzero_count++; }
    }
    *nonzero_indices = NULL;
    if (nonzero_count == 0) { return 0; }

    *nonzero_indices = malloc((size_t)nonzero_count * sizeof(uint32_t));
    if (*nonzero_indices == NULL) {
        fprintf(stderr, "Error allocating UCCD sparse index list\n");
        return -1;
    }
    int output_index = 0;
    for (int i = 0; i < length; i++) {
        if (intensity[i] != 0) {
            (*nonzero_indices)[output_index++] = (uint32_t)i;
        }
    }
    return nonzero_count;
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


static void softargmax_scan_UCCD(float *blur, float *scratch,
                                 const int lengthmz, const int numz,
                                 const float beta)
{
    const int length = lengthmz * numz;
    memcpy(scratch, blur, (size_t)length * sizeof(float));
    for (int i = 0; i < lengthmz; i++) {
        float sum2 = 0;
        float sum1 = 0;
        float factor = 0;
        float min2 = 1.0f;
        for (int j = 0; j < numz; j++) {
            const int index = index2D(numz, i, j);
            const float value = scratch[index];
            const float exponential = expf(beta * value);
            sum1 += value;
            if (exponential < min2) { min2 = exponential; }
            blur[index] = exponential;
            sum2 += exponential;
        }
        const float denominator = sum2 - min2 * (float)numz;
        if (denominator != 0) { factor = sum1 / denominator; }
        for (int j = 0; j < numz; j++) {
            const int index = index2D(numz, i, j);
            if (factor > 0) {
                blur[index] = (blur[index] - min2) * factor;
            } else {
                blur[index] = 0;
            }
        }
    }
}


static void point_smoothing_scan_UCCD(float *blur, float *scratch,
                                      const char *barr, const int lengthmz,
                                      const int numz, const int width)
{
    const int length = lengthmz * numz;
    memcpy(scratch, blur, (size_t)length * sizeof(float));
    const float denominator = 1.0f + 2.0f * (float)width;
    for (int i = 0; i < lengthmz; i++) {
        const int low = i - width > 0 ? i - width : 0;
        const int high = i + width + 1 < lengthmz ? i + width + 1 : lengthmz;
        for (int j = 0; j < numz; j++) {
            const int index = index2D(numz, i, j);
            if (barr[index] == 1) {
                float sum = 0;
                for (int k = low; k < high; k++) {
                    sum += scratch[index2D(numz, k, j)];
                }
                blur[index] = sum / denominator;
            }
        }
    }
}


void setup_blur_z_UCCD(int *zupind, int *zloind, const float *mzdat,
                       const float *zdat, const int scan_length,
                       const float adductmass, const float mzranges[4],
                       const int size[3])
{
    #pragma omp parallel for schedule(static) if(scan_length >= UCCD_OMP_MIN_LENGTH)
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
    #pragma omp parallel for schedule(static) if(scan_length >= UCCD_OMP_MIN_LENGTH)
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
    #pragma omp parallel for schedule(static) if(size[2] > 1 && scan_length * size[2] >= UCCD_OMP_MIN_LENGTH)
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


static int initialize_fft_UCCD(UCCDFFTContext *context, const int size[3],
                               const int batched_2d)
{
    memset(context, 0, sizeof(*context));
    const int scan_length = size[0] * size[1];
    context->real_length = scan_length * size[2];
    context->batched_2d = batched_2d;
    if (batched_2d) {
        context->kernel_length = size[0] * (size[1] / 2 + 1);
        context->spectrum_length = size[2] * context->kernel_length;
        context->normalization = scan_length;
    } else {
        context->kernel_length = size[2] * size[0] * (size[1] / 2 + 1);
        context->spectrum_length = context->kernel_length;
        context->normalization = context->real_length;
    }

    context->real_work = fftwf_malloc((size_t)context->real_length * sizeof(float));
    context->spectrum_work = fftwf_malloc(
        (size_t)context->spectrum_length * sizeof(fftwf_complex));
    if (context->real_work == NULL || context->spectrum_work == NULL) {
        fprintf(stderr, "Error allocating UCCD real FFT workspace\n");
        return 0;
    }

    if (batched_2d) {
        const int dimensions[2] = {size[0], size[1]};
        context->forward_plan = fftwf_plan_many_dft_r2c(
            2, dimensions, size[2], context->real_work, NULL, 1, scan_length,
            context->spectrum_work, NULL, 1, context->kernel_length,
            UCCD_FFT_PLAN_FLAGS);
        context->backward_plan = fftwf_plan_many_dft_c2r(
            2, dimensions, size[2], context->spectrum_work, NULL, 1,
            context->kernel_length, context->real_work, NULL, 1, scan_length,
            UCCD_FFT_PLAN_FLAGS);
    } else {
        context->forward_plan = fftwf_plan_dft_r2c_3d(
            size[2], size[0], size[1], context->real_work,
            context->spectrum_work, UCCD_FFT_PLAN_FLAGS);
        context->backward_plan = fftwf_plan_dft_c2r_3d(
            size[2], size[0], size[1], context->spectrum_work,
            context->real_work, UCCD_FFT_PLAN_FLAGS);
    }
    if (context->forward_plan == NULL || context->backward_plan == NULL) {
        fprintf(stderr, "Error creating UCCD real FFT plans\n");
        return 0;
    }
    return 1;
}


static void destroy_fft_UCCD(UCCDFFTContext *context)
{
    if (context->forward_plan != NULL) { fftwf_destroy_plan(context->forward_plan); }
    if (context->backward_plan != NULL) { fftwf_destroy_plan(context->backward_plan); }
    fftwf_free(context->real_work);
    fftwf_free(context->spectrum_work);
    memset(context, 0, sizeof(*context));
}


static void make_kernel_fft_UCCD(fftwf_complex *kernel_fft,
                                 UCCDFFTContext *context, const int size[3],
                                 const float *chromext, const float *mzext,
                                 const float *zext, const float chromsig,
                                 const float mzsig, const float zsig,
                                 const int psfun, const int zpsfun)
{
    if (context->batched_2d) {
        const int scan_length = size[0] * size[1];
        MakeKernel2D(context->real_work, size, mzext, zext,
                     mzsig, zsig, psfun, zpsfun);
        #pragma omp parallel for schedule(static) if(context->real_length >= UCCD_OMP_MIN_LENGTH)
        for (int scan = 1; scan < size[2]; scan++) {
            memcpy(context->real_work + scan * scan_length,
                   context->real_work, (size_t)scan_length * sizeof(float));
        }
    } else {
        make_kernel3D_UCCD(context->real_work, size, chromext, mzext, zext,
                           chromsig, mzsig, zsig, psfun, zpsfun);
    }
    fftwf_execute(context->forward_plan);
    memcpy(kernel_fft, context->spectrum_work,
           (size_t)context->kernel_length * sizeof(fftwf_complex));
}


static int normalize_kernel_fft_UCCD(fftwf_complex *kernel_fft,
                                     const int kernel_length)
{
    /* The zero-frequency component is the sum of the real-space kernel.
     * Dividing the full spectrum by it makes convolution conserve signal. */
    const float kernel_sum = kernel_fft[0][0];
    if (kernel_sum == 0 || !isfinite(kernel_sum)) {
        fprintf(stderr, "Invalid UCCD reconvolution kernel sum: %g\n", kernel_sum);
        return 0;
    }
    const float inverse_sum = 1.0f / kernel_sum;
    #pragma omp parallel for schedule(static) if(kernel_length >= UCCD_OMP_MIN_LENGTH)
    for (int i = 0; i < kernel_length; i++) {
        kernel_fft[i][0] *= inverse_sum;
        kernel_fft[i][1] *= inverse_sum;
    }
    return 1;
}


static void execute_convolution_UCCD(const fftwf_complex *kernel_fft,
                                     const int conjugate_kernel,
                                     UCCDFFTContext *context,
                                     const int scan_count)
{
    fftwf_execute(context->forward_plan);

    if (context->batched_2d) {
        #pragma omp parallel for schedule(static) if(context->spectrum_length >= UCCD_OMP_MIN_LENGTH)
        for (int scan = 0; scan < scan_count; scan++) {
            const int offset = scan * context->kernel_length;
            for (int i = 0; i < context->kernel_length; i++) {
                const float a = context->spectrum_work[offset + i][0];
                const float b = context->spectrum_work[offset + i][1];
                const float c = kernel_fft[i][0];
                const float d = conjugate_kernel ? -kernel_fft[i][1] : kernel_fft[i][1];
                context->spectrum_work[offset + i][0] = a * c - b * d;
                context->spectrum_work[offset + i][1] = b * c + a * d;
            }
        }
    } else {
        #pragma omp parallel for schedule(static) if(context->spectrum_length >= UCCD_OMP_MIN_LENGTH)
        for (int i = 0; i < context->spectrum_length; i++) {
            const float a = context->spectrum_work[i][0];
            const float b = context->spectrum_work[i][1];
            const float c = kernel_fft[i][0];
            const float d = conjugate_kernel ? -kernel_fft[i][1] : kernel_fft[i][1];
            context->spectrum_work[i][0] = a * c - b * d;
            context->spectrum_work[i][1] = b * c + a * d;
        }
    }
    fftwf_execute(context->backward_plan);
}


static void compute_convolution_UCCD(const float *list,
                                     const fftwf_complex *kernel_fft,
                                     const int conjugate_kernel,
                                     UCCDFFTContext *context,
                                     const int scan_count)
{
    #pragma omp parallel for schedule(static) if(context->real_length >= UCCD_OMP_MIN_LENGTH)
    for (int i = 0; i < context->real_length; i++) {
        context->real_work[i] = list[i];
    }
    execute_convolution_UCCD(kernel_fft, conjugate_kernel, context, scan_count);
}


static void compute_sparse_correction_UCCD(
    const float *data_int, const uint32_t *nonzero_indices,
    const int nonzero_count, const fftwf_complex *kernel_fft,
    UCCDFFTContext *context, const int scan_count, const float fft_scale,
    float *ratio_values)
{
    /* The Richardson-Lucy ratio is zero wherever the original input is zero.
     * Save only its nonzero candidates before reusing the FFT output buffer. */
    #pragma omp parallel for schedule(static) if(nonzero_count >= UCCD_OMP_MIN_LENGTH)
    for (int i = 0; i < nonzero_count; i++) {
        const uint32_t index = nonzero_indices[i];
        const float predicted = context->real_work[index] * fft_scale;
        ratio_values[i] = predicted != 0 ? data_int[index] / predicted : 0;
    }

    memset(context->real_work, 0,
           (size_t)context->real_length * sizeof(float));
    for (int i = 0; i < nonzero_count; i++) {
        context->real_work[nonzero_indices[i]] = ratio_values[i];
    }
    execute_convolution_UCCD(kernel_fft, 1, context, scan_count);
}


static void fftconvolve_precomputed_UCCD(float *output, const float *input,
                                         const fftwf_complex *kernel_fft,
                                         UCCDFFTContext *context,
                                         const int scan_count)
{
    compute_convolution_UCCD(input, kernel_fft, 0, context, scan_count);
    const float scale = 1.0f / (float)context->normalization;
    #pragma omp parallel for schedule(static) if(context->real_length >= UCCD_OMP_MIN_LENGTH)
    for (int i = 0; i < context->real_length; i++) {
        output[i] = context->real_work[i] * scale;
    }
}


int run_unidec_UCCD(int argc, char *argv[], Config config)
{
    const time_t starttime = time(NULL);
    char input_filename[550];
    char decon_filename[550];
    if (!make_filename_UCCD(input_filename, sizeof(input_filename),
                            config.outfile, UCCD_INPUT_SUFFIX) ||
        !make_filename_UCCD(decon_filename, sizeof(decon_filename),
                            config.outfile, UCCD_DECON_SUFFIX)) {
        return 1;
    }
    /* Fit output is temporarily disabled.
    char fit_filename[550];
    if (!make_filename_UCCD(fit_filename, sizeof(fit_filename),
                            config.outfile, UCCD_FIT_SUFFIX)) {
        return 1;
    }
    */

    printf("Opening sparse binary UCCD file: %s\n", input_filename);

    int size[3] = {0, 0, 0};
    float *chromext = NULL;
    float *mzext = NULL;
    float *zext = NULL;
    float *dataInt = NULL;
    if (!read_sparse_UCCD(input_filename, size, &chromext, &mzext, &zext, &dataInt)) {
        return 2;
    }

    const int scan_length = size[0] * size[1];
    const int output_size[3] = {size[0], size[1], size[2]};
    const int output_lines = scan_length * output_size[2];
    printf("Dimensions: %d chromatography by %d m/z by %d charge: %d total\n",
           output_size[2], size[0], size[1], output_lines);

    uint32_t *nonzero_indices = NULL;
    const int nonzero_count = collect_nonzero_indices_UCCD(
        dataInt, output_lines, &nonzero_indices);
    if (nonzero_count <= 0) {
        fprintf(stderr, nonzero_count == 0 ? "UCCD input contains no signal\n" :
                                             "Error collecting UCCD sparse indexes\n");
        fftwf_free(dataInt);
        free(chromext); free(mzext); free(zext);
        return 1;
    }

    const int batched_2d = config.dtsig == 0;
    int compacted_scans = 0;
    int *active_scans = NULL;
    if (batched_2d) {
        char *active_flags = calloc((size_t)output_size[2], sizeof(char));
        if (active_flags == NULL) {
            fprintf(stderr, "Error allocating UCCD active-scan flags\n");
            free(nonzero_indices); fftwf_free(dataInt);
            free(chromext); free(mzext); free(zext);
            return 1;
        }
        int active_count = 0;
        for (int i = 0; i < nonzero_count; i++) {
            const int scan = (int)(nonzero_indices[i] / (uint32_t)scan_length);
            if (!active_flags[scan]) {
                active_flags[scan] = 1;
                active_count++;
            }
        }

        if (active_count < output_size[2]) {
            active_scans = malloc((size_t)active_count * sizeof(int));
            int *scan_map = malloc((size_t)output_size[2] * sizeof(int));
            float *compact_data = fftwf_malloc(
                (size_t)active_count * (size_t)scan_length * sizeof(float));
            if (active_scans == NULL || scan_map == NULL || compact_data == NULL) {
                fprintf(stderr, "Error allocating compact UCCD scan arrays\n");
                free(active_flags); free(active_scans); free(scan_map);
                fftwf_free(compact_data); free(nonzero_indices); fftwf_free(dataInt);
                free(chromext); free(mzext); free(zext);
                return 1;
            }
            for (int scan = 0; scan < output_size[2]; scan++) { scan_map[scan] = -1; }
            int packed_scan = 0;
            for (int scan = 0; scan < output_size[2]; scan++) {
                if (active_flags[scan]) {
                    active_scans[packed_scan] = scan;
                    scan_map[scan] = packed_scan++;
                }
            }
            memset(compact_data, 0,
                   (size_t)active_count * (size_t)scan_length * sizeof(float));
            for (int i = 0; i < nonzero_count; i++) {
                const uint32_t old_index = nonzero_indices[i];
                const int original_scan = (int)(old_index / (uint32_t)scan_length);
                const uint32_t scan_index = old_index % (uint32_t)scan_length;
                const uint32_t compact_index =
                    (uint32_t)(scan_map[original_scan] * scan_length) + scan_index;
                compact_data[compact_index] = dataInt[old_index];
                nonzero_indices[i] = compact_index;
            }
            fftwf_free(dataInt);
            dataInt = compact_data;
            size[2] = active_count;
            compacted_scans = 1;
            free(scan_map);
        }
        printf("Active-scan FFT compaction: %d of %d scans (%d skipped)\n",
               size[2], output_size[2], output_size[2] - size[2]);
        free(active_flags);
    }

    const int lines = scan_length * size[2];

    /* The mass/charge smoothing setup only needs one dense m/z-charge scan. */
    float *mzdat = calloc((size_t)scan_length, sizeof(float));
    float *zdat = calloc((size_t)scan_length, sizeof(float));
    if (mzdat == NULL || zdat == NULL) {
        fprintf(stderr, "Error allocating UCCD coordinate grid\n");
        return 1;
    }
    for (int mzindex = 0; mzindex < size[0]; mzindex++) {
        for (int zindex = 0; zindex < size[1]; zindex++) {
            const int index = mzindex * size[1] + zindex;
            mzdat[index] = mzext[mzindex];
            zdat[index] = zext[zindex];
        }
    }
    const float mzranges[4] = {mzext[0], mzext[size[0] - 1],
                               zext[0], zext[size[1] - 1]};
    printf("Chromatography range: %f to %f\n",
           chromext[0], chromext[output_size[2] - 1]);
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

    const int explicit_thread_count = argc > 3 && strcmp(argv[2], "-nthreads") == 0;
#ifdef _OPENMP
    int thread_count = omp_get_max_threads();
    if (!explicit_thread_count && thread_count > 8) {
        thread_count = 8;
        omp_set_num_threads(thread_count);
    }
    omp_set_dynamic(0);
#else
    const int thread_count = 1;
#endif
    printf("UCCD worker threads: %d%s\n", thread_count,
           explicit_thread_count ? " (requested)" : "");

#ifdef UNIDEC_USE_MKL
    const int previous_mkl_threads = mkl_set_num_threads_local(thread_count);
#endif
#ifdef UNIDEC_USE_FFTW_THREADS
    const int fftw_threads_initialized = fftwf_init_threads();
    if (!fftw_threads_initialized) {
        fprintf(stderr, "Error initializing threaded FFTW for UCCD\n");
        return 1;
    }
    fftwf_plan_with_nthreads(thread_count);
#endif

    float *blur = fftwf_malloc((size_t)lines * sizeof(float));
    float *newblur = fftwf_malloc((size_t)lines * sizeof(float));
    float *oldblur = fftwf_malloc((size_t)lines * sizeof(float));
    float *ratio_values = malloc((size_t)nonzero_count * sizeof(float));
    /* Fit output is temporarily disabled.
    float *newblur2 = fftwf_malloc((size_t)lines * sizeof(float));
    */
    if (blur == NULL || newblur == NULL || oldblur == NULL || ratio_values == NULL) {
        fprintf(stderr, "Error allocating UCCD processing arrays\n");
        return 1;
    }

    UCCDFFTContext fft_context;
    if (!initialize_fft_UCCD(&fft_context, size, batched_2d)) {
        destroy_fft_UCCD(&fft_context);
        return 1;
    }
    fftwf_complex *kernel_fft = fftwf_malloc(
        (size_t)fft_context.kernel_length * sizeof(fftwf_complex));
    if (kernel_fft == NULL) {
        fprintf(stderr, "Error allocating UCCD kernel spectrum\n");
        destroy_fft_UCCD(&fft_context);
        return 1;
    }

    printf("Peak widths: chromatography %f, m/z %f, charge %f\n",
           config.dtsig, config.mzsig, config.csig);
    printf("FFT mode: %s real transforms\n",
           batched_2d ? "batched 2-D" : "coupled 3-D");
    make_kernel_fft_UCCD(kernel_fft, &fft_context, size, chromext, mzext, zext,
                         config.dtsig, config.mzsig, config.csig,
                         config.psfun, config.zpsfun);

    const size_t matsize = (size_t)lines * sizeof(float);
    memcpy(blur, dataInt, matsize);
    memcpy(oldblur, blur, matsize);
    printf("Iterating.");
    double conv = 0;
    int off = 0;
    for (int iteration = 0; iteration < config.numit; iteration++) {
        /* These operations intentionally do not cross chromatography scans. */
        if (config.beta > 0 || config.psig > 0 || config.zsig != 0 || config.msig != 0) {
            #pragma omp parallel for schedule(static) if(size[2] > 1 && lines >= UCCD_OMP_MIN_LENGTH)
            for (int scan = 0; scan < size[2]; scan++) {
                const int offset = scan * scan_length;
                float *const scratch = newblur + offset;
                if (config.beta > 0) {
                    softargmax_scan_UCCD(blur + offset, scratch, size[0], size[1],
                                         config.beta / betafactor);
                }
                if (config.psig > 0) {
                    point_smoothing_scan_UCCD(blur + offset, scratch, barr,
                                              size[0], size[1], abs((int)config.psig));
                }
                if (config.zsig != 0) {
                    blur_it_UCCD(scratch, blur + offset, zupind, zloind,
                                 scan_length, config.zsig * dmax);
                    memcpy(blur + offset, scratch,
                           (size_t)scan_length * sizeof(float));
                }
                if (config.msig != 0) {
                    blur_it_UCCD(scratch, blur + offset, mupind, mloind,
                                 scan_length, config.msig * dmax);
                    memcpy(blur + offset, scratch,
                           (size_t)scan_length * sizeof(float));
                }
            }
        }

        compute_convolution_UCCD(blur, kernel_fft, 0, &fft_context, size[2]);
        const float fft_scale = 1.0f / (float)fft_context.normalization;
        compute_sparse_correction_UCCD(
            dataInt, nonzero_indices, nonzero_count, kernel_fft, &fft_context,
            size[2], fft_scale, ratio_values);
        #pragma omp parallel for schedule(static) if(lines >= UCCD_OMP_MIN_LENGTH)
        for (int i = 0; i < lines; i++) {
            blur[i] *= fft_context.real_work[i] * fft_scale;
        }

        if (config.numit < 10 || iteration % 10 == 0 || iteration % 10 == 1 ||
            iteration > 0.9 * config.numit) {
            double diff = 0;
            double total = 0;
            #pragma omp parallel for reduction(+:diff,total) schedule(static) if(lines >= UCCD_OMP_MIN_LENGTH)
            for (int i = 0; i < lines; i++) {
                const double delta = (double)blur[i] - (double)oldblur[i];
                diff += delta * delta;
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
                    printf("Converged in %d iterations.\n", iteration + 1);
                    break;
                }
                off = 1;
            }
            memcpy(oldblur, blur, matsize);
        }
    }
    printf("Completed iterations\n");

    int result = 0;
    float *expanded_output = NULL;
    /* Fit calculation and _uccd_fit.bin output are temporarily disabled.
    fftconvolve_precomputed_UCCD(newblur2, blur, kernel_fft,
                                &fft_context, size[2]);
    if (config.datanorm == 1) {
        const float fitmax = Max(newblur2, lines);
        if (dmax != 0 && fitmax != 0) { Normalize(lines, newblur2, fitmax / dmax); }
    }
    ApplyCutoff(newblur2, 0, lines);
    if (!write_sparse_UCCD(fit_filename, size, chromext, mzext, zext, newblur2)) {
        result = 1;
        goto cleanup_processing_UCCD;
    }
    */

    if (config.rawflag == 0) {
        make_kernel_fft_UCCD(kernel_fft, &fft_context, size, chromext, mzext, zext,
                             config.dtsig, config.mzsig, 0,
                             config.psfun, config.zpsfun);
        if (!normalize_kernel_fft_UCCD(kernel_fft, fft_context.kernel_length)) {
            result = 1;
            goto cleanup_processing_UCCD;
        }
        fftconvolve_precomputed_UCCD(blur, blur, kernel_fft,
                                    &fft_context, size[2]);
        printf("Reconvolved with chromatography and m/z dimensions\n");
    }
    if (config.datanorm == 1) {
        const float blurmax = Max(blur, lines);
        if (dmax != 0 && blurmax != 0) { Normalize(lines, blur, blurmax / dmax); }
    }
    ApplyCutoff(blur, 0, lines);
    const float *output_data = blur;
    if (compacted_scans) {
        expanded_output = fftwf_malloc((size_t)output_lines * sizeof(float));
        if (expanded_output == NULL) {
            fprintf(stderr, "Error allocating expanded UCCD output\n");
            result = 1;
            goto cleanup_processing_UCCD;
        }
        memset(expanded_output, 0, (size_t)output_lines * sizeof(float));
        #pragma omp parallel for schedule(static) if(size[2] > 1 && lines >= UCCD_OMP_MIN_LENGTH)
        for (int packed_scan = 0; packed_scan < size[2]; packed_scan++) {
            memcpy(expanded_output + active_scans[packed_scan] * scan_length,
                   blur + packed_scan * scan_length,
                   (size_t)scan_length * sizeof(float));
        }
        output_data = expanded_output;
    }
    if (!write_sparse_UCCD(decon_filename, output_size, chromext, mzext, zext,
                           output_data)) {
        result = 1;
    }

cleanup_processing_UCCD:
    fftwf_free(expanded_output);
    fftwf_free(kernel_fft);
    destroy_fft_UCCD(&fft_context);
    fftwf_free(blur); fftwf_free(newblur); fftwf_free(oldblur);
    /* Fit output is temporarily disabled.
    fftwf_free(newblur2);
    */
    free(ratio_values);
#ifdef UNIDEC_USE_FFTW_THREADS
    fftwf_cleanup_threads();
#endif
#ifdef UNIDEC_USE_MKL
    mkl_set_num_threads_local(previous_mkl_threads);
#endif

    free(mzdat); free(zdat); fftwf_free(dataInt);
    free(chromext); free(mzext); free(zext);
    free(zupind); free(zloind); free(mupind); free(mloind); free(barr);
    free(nonzero_indices); free(active_scans);

    printf("Done in %ds!\n", (int)difftime(time(NULL), starttime));
    return result;
}
