/*
 * Coupled UniChrom deconvolution on ragged processed m/z axes.
 *
 * The direct response is stored in forward and transposed CSR form. This
 * preserves the Richardson-Lucy forward/adjoint pair without interpolating
 * observations onto a common axis before solving.
 */

#include "UniChrom_Nonlinear.h"
#include "UniChrom_Main.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#define UC_DIRECT_OMP_MIN_LENGTH 32768

typedef struct {
    int scan_count;
    int point_count;
    int *scan_offsets;
    /* Scan numbers now; this is the single axis to replace with acquisition
     * times when a time-based dtsig mode is added. */
    float *chrom_axis;
    float *mz;
    float *observed;
} UCRaggedData;

typedef struct {
    int length;
    size_t edge_count;
    size_t *offsets;
    int *indices;
    float *weights;
    size_t *transpose_offsets;
    int *transpose_indices;
    float *transpose_weights;
    float *sensitivity;
} UCDirectOperator;


static int lower_bound_UC(const float *axis, const int length, const float value)
{
    int low = 0;
    int high = length;
    while (low < high) {
        const int middle = low + (high - low) / 2;
        if (axis[middle] < value) { low = middle + 1; }
        else { high = middle; }
    }
    return low;
}


static int mirrored_scan_UC(int scan, const int scan_count)
{
    if (scan_count == 1) { return 0; }
    const int period = 2 * scan_count;
    scan %= period;
    if (scan < 0) { scan += period; }
    return scan < scan_count ? scan : period - scan - 1;
}


static float periodic_scan_peak_UC(const int length, const int index,
                                   const float width, const int shape)
{
    if (length == 1) { return 1; }
    return mzpeakshape(0, (float)index, width, shape) +
           mzpeakshape((float)length, (float)index, width, shape);
}


static int read_processed_UC(const Config config, UCRaggedData *data)
{
    memset(data, 0, sizeof(*data));
    data->scan_count = int_attr(config.file_id, "/ms_dataset", "num", 0);
    data->scan_offsets = calloc((size_t)data->scan_count + 1, sizeof(int));
    data->chrom_axis = calloc((size_t)data->scan_count, sizeof(float));
    if (data->scan_offsets == NULL || data->chrom_axis == NULL) { return 0; }

    int64_t total = 0;
    for (int scan = 0; scan < data->scan_count; scan++) {
        char path[1024];
        snprintf(path, sizeof(path), "/ms_dataset/%d/processed_data", scan);
        const int length = mh5getfilelength(config.file_id, path);
        if (length < 0 || total + length > INT_MAX) {
            fprintf(stderr, "Invalid nonlinear UniChrom processed spectrum %d\n", scan);
            return 0;
        }
        data->scan_offsets[scan] = (int)total;
        data->chrom_axis[scan] = (float)scan;
        total += length;
    }
    data->scan_offsets[data->scan_count] = (int)total;
    data->point_count = (int)total;
    if (data->point_count < 1) {
        fprintf(stderr, "Nonlinear UniChrom has no processed m/z points\n");
        return 0;
    }
    data->mz = calloc((size_t)data->point_count, sizeof(float));
    data->observed = calloc((size_t)data->point_count, sizeof(float));
    if (data->mz == NULL || data->observed == NULL) { return 0; }

    for (int scan = 0; scan < data->scan_count; scan++) {
        char path[1024];
        snprintf(path, sizeof(path), "/ms_dataset/%d/processed_data", scan);
        const int first = data->scan_offsets[scan];
        const int length = data->scan_offsets[scan + 1] - first;
        if (length == 0) { continue; }
        mh5readfile2d(config.file_id, path, length,
                      data->mz + first, data->observed + first);
        for (int i = first + 1; i < first + length; i++) {
            if (!(data->mz[i] > data->mz[i - 1])) {
                fprintf(stderr, "Nonlinear UniChrom m/z axis %d is not strictly increasing\n", scan);
                return 0;
            }
        }
    }
    return 1;
}


static void free_ragged_UC(UCRaggedData *data)
{
    free(data->scan_offsets);
    free(data->chrom_axis);
    free(data->mz);
    free(data->observed);
    memset(data, 0, sizeof(*data));
}


static int scan_for_point_UC(const UCRaggedData *data, const int point)
{
    int low = 0;
    int high = data->scan_count;
    while (low + 1 < high) {
        const int middle = low + (high - low) / 2;
        if (data->scan_offsets[middle] <= point) { low = middle; }
        else { high = middle; }
    }
    return low;
}


static int interpolation_neighbors_UC(const float *axis, const int length,
                                      const float value, int indexes[2],
                                      float weights[2])
{
    if (length < 1) { return 0; }
    if (value < axis[0] || value > axis[length - 1]) { return 0; }
    const int upper = lower_bound_UC(axis, length, value);
    if (upper == 0 || upper == length || axis[upper] == value) {
        indexes[0] = upper == length ? length - 1 : upper;
        weights[0] = 1;
        return 1;
    }
    const int lower = upper - 1;
    const float fraction = (value - axis[lower]) / (axis[upper] - axis[lower]);
    indexes[0] = lower;
    indexes[1] = upper;
    weights[0] = 1 - fraction;
    weights[1] = fraction;
    return 2;
}


static size_t count_mz_edges_UC(const UCRaggedData *data,
                                const int source_scan, const float target_mz,
                                const Config config)
{
    const int first = data->scan_offsets[source_scan];
    const int length = data->scan_offsets[source_scan + 1] - first;
    const float *axis = data->mz + first;
    if (config.mzsig == 0) {
        int indexes[2]; float weights[2];
        return (size_t)interpolation_neighbors_UC(axis, length, target_mz,
                                                   indexes, weights);
    }
    const float window = config.psmzthresh > 0 ? config.psmzthresh : 6 * fabsf(config.mzsig);
    const int begin = lower_bound_UC(axis, length, target_mz - window);
    const int end = lower_bound_UC(axis, length, target_mz + window);
    return end > begin ? (size_t)(end - begin) : 0;
}


static size_t write_mz_edges_UC(const UCRaggedData *data,
                                const int source_scan, const float target_mz,
                                const float time_weight, const Config config,
                                int *indices, float *weights)
{
    const int first = data->scan_offsets[source_scan];
    const int length = data->scan_offsets[source_scan + 1] - first;
    const float *axis = data->mz + first;
    if (config.mzsig == 0) {
        int local[2]; float interpolation[2];
        const int count = interpolation_neighbors_UC(axis, length, target_mz,
                                                       local, interpolation);
        for (int i = 0; i < count; i++) {
            indices[i] = first + local[i];
            weights[i] = time_weight * interpolation[i];
        }
        return (size_t)count;
    }
    const float window = config.psmzthresh > 0 ? config.psmzthresh : 6 * fabsf(config.mzsig);
    const int begin = lower_bound_UC(axis, length, target_mz - window);
    const int end = lower_bound_UC(axis, length, target_mz + window);
    size_t written = 0;
    for (int source = begin; source < end; source++) {
        const float weight = time_weight * mzpeakshape(
            axis[source], target_mz,
            fabsf(config.mzsig) * config.peakshapeinflate, config.psfun);
        if (weight > 0 && isfinite(weight)) {
            indices[written] = first + source;
            weights[written++] = weight;
        }
    }
    return written;
}


static int build_direct_operator_UC(const UCRaggedData *data,
                                    const Config config, const int zero_padding,
                                    UCDirectOperator *op)
{
    memset(op, 0, sizeof(*op));
    op->length = data->point_count;
    const double padding_value = ceil(3.0 * config.dtsig);
    if (padding_value > (INT_MAX - data->scan_count) / 2.0) { return 0; }
    const int padding = (int)padding_value;
    const int padded_count = data->scan_count + 2 * padding;
    float *time_weights = calloc((size_t)data->scan_count * data->scan_count,
                                 sizeof(float));
    op->offsets = calloc((size_t)op->length + 1, sizeof(size_t));
    if (time_weights == NULL || op->offsets == NULL) { free(time_weights); return 0; }

    /* Keep construction of the temporal matrix separate from the ragged m/z
     * response so a later mode can calculate these weights from chrom_axis. */
    for (int target = 0; target < data->scan_count; target++) {
        const int target_extended = target + padding;
        for (int extended = 0; extended < padded_count; extended++) {
            int source = extended - padding;
            if (zero_padding && (source < 0 || source >= data->scan_count)) { continue; }
            if (!zero_padding) { source = mirrored_scan_UC(source, data->scan_count); }
            int delta = target_extended - extended;
            if (delta < 0) { delta += padded_count; }
            time_weights[index2D(data->scan_count, target, source)] +=
                periodic_scan_peak_UC(padded_count, delta, config.dtsig, config.psfun);
        }
    }

    size_t edge_count = 0;
    for (int target_point = 0; target_point < op->length; target_point++) {
        const int target_scan = scan_for_point_UC(data, target_point);
        op->offsets[target_point] = edge_count;
        for (int source_scan = 0; source_scan < data->scan_count; source_scan++) {
            const float time_weight = time_weights[index2D(
                data->scan_count, target_scan, source_scan)];
            if (!(time_weight > 0)) { continue; }
            const size_t count = count_mz_edges_UC(
                data, source_scan, data->mz[target_point], config);
            if (SIZE_MAX - edge_count < count) { free(time_weights); return 0; }
            edge_count += count;
        }
    }
    op->offsets[op->length] = edge_count;
    op->edge_count = edge_count;
    op->indices = malloc(edge_count * sizeof(int));
    op->weights = malloc(edge_count * sizeof(float));
    if (op->indices == NULL || op->weights == NULL) {
        free(time_weights); return 0;
    }

    for (int target_point = 0; target_point < op->length; target_point++) {
        const int target_scan = scan_for_point_UC(data, target_point);
        size_t edge = op->offsets[target_point];
        for (int source_scan = 0; source_scan < data->scan_count; source_scan++) {
            const float time_weight = time_weights[index2D(
                data->scan_count, target_scan, source_scan)];
            if (!(time_weight > 0)) { continue; }
            edge += write_mz_edges_UC(data, source_scan, data->mz[target_point],
                                      time_weight, config,
                                      op->indices + edge, op->weights + edge);
        }
        op->offsets[target_point + 1] = edge;
    }
    op->edge_count = op->offsets[op->length];
    op->transpose_offsets = calloc((size_t)op->length + 1, sizeof(size_t));
    if (op->transpose_offsets == NULL) { free(time_weights); return 0; }
    for (size_t edge = 0; edge < op->edge_count; edge++) {
        op->transpose_offsets[op->indices[edge] + 1]++;
    }
    for (int i = 0; i < op->length; i++) {
        op->transpose_offsets[i + 1] += op->transpose_offsets[i];
    }
    op->transpose_indices = malloc(op->edge_count * sizeof(int));
    op->transpose_weights = malloc(op->edge_count * sizeof(float));
    size_t *positions = malloc((size_t)op->length * sizeof(size_t));
    op->sensitivity = calloc((size_t)op->length, sizeof(float));
    if (op->transpose_indices == NULL || op->transpose_weights == NULL ||
        positions == NULL || op->sensitivity == NULL) {
        free(time_weights); free(positions); return 0;
    }
    memcpy(positions, op->transpose_offsets, (size_t)op->length * sizeof(size_t));
    for (int row = 0; row < op->length; row++) {
        for (size_t edge = op->offsets[row]; edge < op->offsets[row + 1]; edge++) {
            const int column = op->indices[edge];
            const size_t position = positions[column]++;
            op->transpose_indices[position] = row;
            op->transpose_weights[position] = op->weights[edge];
            op->sensitivity[column] += op->weights[edge];
        }
    }
    free(positions);
    free(time_weights);
    for (int i = 0; i < op->length; i++) {
        if (!(op->sensitivity[i] > 0) || op->offsets[i] == op->offsets[i + 1]) {
            fprintf(stderr, "Nonlinear UniChrom response has an unsupported m/z point\n");
            return 0;
        }
    }
    printf("UniChrom nonlinear response: %d points, %zu direct edges\n",
           op->length, op->edge_count);
    return 1;
}


static void free_direct_operator_UC(UCDirectOperator *op)
{
    free(op->offsets); free(op->indices); free(op->weights);
    free(op->transpose_offsets); free(op->transpose_indices);
    free(op->transpose_weights); free(op->sensitivity);
    memset(op, 0, sizeof(*op));
}


static void forward_UC(const UCDirectOperator *op, const float *input, float *output)
{
    #pragma omp parallel for schedule(static) if(op->length >= UC_DIRECT_OMP_MIN_LENGTH)
    for (int row = 0; row < op->length; row++) {
        float value = 0;
        for (size_t edge = op->offsets[row]; edge < op->offsets[row + 1]; edge++) {
            value += op->weights[edge] * input[op->indices[edge]];
        }
        output[row] = value;
    }
}


static void adjoint_UC(const UCDirectOperator *op, const float *input, float *output)
{
    #pragma omp parallel for schedule(static) if(op->length >= UC_DIRECT_OMP_MIN_LENGTH)
    for (int row = 0; row < op->length; row++) {
        float value = 0;
        for (size_t edge = op->transpose_offsets[row];
             edge < op->transpose_offsets[row + 1]; edge++) {
            value += op->transpose_weights[edge] * input[op->transpose_indices[edge]];
        }
        output[row] = value;
    }
}


static int nearest_in_scan_UC(const UCRaggedData *data, const int scan,
                              const float value)
{
    const int first = data->scan_offsets[scan];
    const int length = data->scan_offsets[scan + 1] - first;
    if (length < 1) { return -1; }
    if (value < data->mz[first] || value > data->mz[first + length - 1]) { return -1; }
    return first + nearfast(data->mz + first, value, length);
}


static void make_regularizer_indexes_UC(const UCRaggedData *data,
                                        const Config config,
                                        int *zup, int *zdown,
                                        int *mup, int *mdown)
{
    const int charge_count = config.numz;
    #pragma omp parallel for schedule(static) if(data->point_count * charge_count >= UC_DIRECT_OMP_MIN_LENGTH)
    for (int point = 0; point < data->point_count; point++) {
        const int scan = scan_for_point_UC(data, point);
        const float mz = data->mz[point];
        for (int zi = 0; zi < charge_count; zi++) {
            const int index = point * charge_count + zi;
            const int charge = config.startz + zi;
            const float mass = calcmass(mz, charge, config.adductmass);
            if (zup != NULL) {
                const int mapped = nearest_in_scan_UC(
                    data, scan, calcmz(mass, charge + 1, config.adductmass));
                zup[index] = zi == charge_count - 1 || mapped < 0 ? index :
                    mapped * charge_count + zi + 1;
                const int lower = charge == 1 ? -1 : nearest_in_scan_UC(
                    data, scan, calcmz(mass, charge - 1, config.adductmass));
                zdown[index] = zi == 0 || lower < 0 ? index :
                    lower * charge_count + zi - 1;
            }
            if (mup != NULL) {
                const int upper = nearest_in_scan_UC(
                    data, scan, calcmz(mass + config.molig, charge, config.adductmass));
                const int lower = nearest_in_scan_UC(
                    data, scan, calcmz(mass - config.molig, charge, config.adductmass));
                mup[index] = upper < 0 ? index : upper * charge_count + zi;
                mdown[index] = lower < 0 ? index : lower * charge_count + zi;
            }
        }
    }
}


static void clear_cube_UC(float *cube, const char *allowed, const int length)
{
    #pragma omp parallel for schedule(static) if(length >= UC_DIRECT_OMP_MIN_LENGTH)
    for (int i = 0; i < length; i++) {
        if (!allowed[i] || !isfinite(cube[i]) || cube[i] < 0) { cube[i] = 0; }
    }
}


static void scale_sum_UC(float *sum, const int sum_length,
                         const float *grid, const int grid_length,
                         const int datanorm)
{
    if (datanorm == 1) { norm1d(sum, sum_length); return; }
    const float grid_max = Max(grid, grid_length);
    const float sum_max = Max(sum, sum_length);
    if (sum_max > 0) { Normalize(sum_length, sum, sum_max / grid_max); }
}


static int write_outputs_UC(const Config config, const UCRaggedData *data,
                            const UCDirectOperator *op, const float *cube,
                            const char *allowed)
{
    const int charge_count = config.numz;
    const int cube_length = data->point_count * charge_count;
    float *output = malloc((size_t)cube_length * sizeof(float));
    float *plane = calloc((size_t)data->point_count, sizeof(float));
    float *convolved = calloc((size_t)data->point_count, sizeof(float));
    if (output == NULL || plane == NULL || convolved == NULL) {
        free(output); free(plane); free(convolved); return 0;
    }
    memcpy(output, cube, (size_t)cube_length * sizeof(float));
    if (config.rawflag == 0 || config.rawflag == 2) {
        for (int charge = 0; charge < charge_count; charge++) {
            for (int point = 0; point < data->point_count; point++) {
                plane[point] = cube[(size_t)point * charge_count + charge];
            }
            forward_UC(op, plane, convolved);
            for (int point = 0; point < data->point_count; point++) {
                output[(size_t)point * charge_count + charge] =
                    convolved[point];
            }
        }
        clear_cube_UC(output, allowed, cube_length);
    }
    if (config.datanorm == 1) {
        const float maximum = Max(output, cube_length);
        const float data_max = Max(data->observed, data->point_count);
        if (maximum > 0 && data_max > 0) { Normalize(cube_length, output, maximum / data_max); }
    }
    ApplyCutoff(output, 0, cube_length);

    float mz_min = INFINITY;
    float mz_max = -INFINITY;
    float mz_step = INFINITY;
    for (int scan = 0; scan < data->scan_count; scan++) {
        const int first = data->scan_offsets[scan];
        const int count = data->scan_offsets[scan + 1] - first;
        if (count == 0) { continue; }
        if (data->mz[first] < mz_min) { mz_min = data->mz[first]; }
        if (data->mz[first + count - 1] > mz_max) { mz_max = data->mz[first + count - 1]; }
        if (count > 1) {
            const float step = (data->mz[first + count - 1] - data->mz[first]) /
                               (float)count;
            if (step > 0 && step < mz_step) { mz_step = step; }
        }
    }
    if (!isfinite(mz_step)) { mz_step = mz_max > mz_min ? mz_max - mz_min : 1; }
    const double mz_bins = ((double)mz_max - mz_min) / mz_step;
    if (!isfinite(mz_bins) || mz_bins > INT_MAX - 1 ||
        (int64_t)data->scan_count * (1 + (int)mz_bins) > INT_MAX) {
        fprintf(stderr, "Nonlinear UniChrom merged m/z grid is too large\n");
        free(output); free(plane); free(convolved);
        return 0;
    }
    const int mz_count = 1 + (int)mz_bins;
    float *mz_axis = calloc((size_t)mz_count, sizeof(float));
    float *mz_grid = calloc((size_t)data->scan_count * mz_count, sizeof(float));
    float *mz_sum = calloc((size_t)mz_count, sizeof(float));
    if (mz_axis == NULL || mz_grid == NULL || mz_sum == NULL) {
        free(output); free(plane); free(convolved);
        free(mz_axis); free(mz_grid); free(mz_sum); return 0;
    }
    for (int i = 0; i < mz_count; i++) { mz_axis[i] = mz_min + i * mz_step; }

    float cube_max = Max(output, cube_length);
    const float cutoff = cube_max * .000001f;
    float mass_min = config.fixedmassaxis ? config.masslb : config.massub;
    float mass_max = config.fixedmassaxis ? config.massub : config.masslb;
    if (!config.fixedmassaxis) {
        for (int point = 0; point < data->point_count; point++) {
            for (int zi = 0; zi < charge_count; zi++) {
                const int index = point * charge_count + zi;
                if (allowed[index] && output[index] > cutoff) {
                    const float mass = calcmass(data->mz[point], config.startz + zi,
                                                config.adductmass);
                    const float low = mass - config.psmzthresh * (config.startz + zi);
                    const float high = mass + config.psmzthresh * (config.startz + zi) + config.massbins;
                    if (low < mass_min) { mass_min = low; }
                    if (high > mass_max) { mass_max = high; }
                }
            }
        }
        if (mass_max <= mass_min) { mass_min = config.masslb; mass_max = config.massub; }
        mass_min = floorf(mass_min / config.massbins) * config.massbins;
        mass_max = ceilf(mass_max / config.massbins) * config.massbins;
    }
    const double mass_bins = ((double)mass_max - mass_min) / config.massbins;
    if (!isfinite(mass_bins) || mass_bins < 0 || mass_bins > INT_MAX - 1 ||
        (int64_t)data->scan_count * (1 + (int)mass_bins) > INT_MAX) {
        fprintf(stderr, "Nonlinear UniChrom merged mass grid is too large\n");
        free(output); free(plane); free(convolved);
        free(mz_axis); free(mz_grid); free(mz_sum);
        return 0;
    }
    const int mass_count = 1 + (int)mass_bins;
    float *mass_axis = calloc((size_t)mass_count, sizeof(float));
    float *mass_grid = calloc((size_t)data->scan_count * mass_count, sizeof(float));
    float *mass_sum = calloc((size_t)mass_count, sizeof(float));
    float *merged_cube = calloc((size_t)mz_count * charge_count, sizeof(float));
    float *merged_charge = calloc((size_t)mz_count, sizeof(float));
    if (mass_axis == NULL || mass_grid == NULL || mass_sum == NULL ||
        merged_cube == NULL || merged_charge == NULL) {
        free(output); free(plane); free(convolved);
        free(mz_axis); free(mz_grid); free(mz_sum);
        free(mass_axis); free(mass_grid); free(mass_sum);
        free(merged_cube); free(merged_charge);
        return 0;
    }
    for (int i = 0; i < mass_count; i++) { mass_axis[i] = mass_min + i * config.massbins; }

    for (int scan = 0; scan < data->scan_count; scan++) {
        const int first = data->scan_offsets[scan];
        const int point_count = data->scan_offsets[scan + 1] - first;
        float *mass_values = mass_grid + (size_t)scan * mass_count;
        float *mz_values = mz_grid + (size_t)scan * mz_count;
        if (point_count == 0) {
            char group[1024], path[1024];
            snprintf(group, sizeof(group), "/ms_dataset/%d", scan);
            snprintf(path, sizeof(path), "%s/mass_data", group);
            mh5writefile2d(config.file_id, path, mass_count, mass_axis, mass_values);
            write_attr_int(config.file_id, group, "length_mz", 0);
            write_attr_int(config.file_id, group, "length_mass", mass_count);
            write_attr_float(config.file_id, group, "mzsig", config.mzsig);
            write_attr_float(config.file_id, group, "dtsig", config.dtsig);
            continue;
        }
        memset(merged_cube, 0, (size_t)mz_count * charge_count * sizeof(float));
        for (int zi = 0; zi < charge_count; zi++) {
            for (int local = 0; local < point_count; local++) {
                const int index = (first + local) * charge_count + zi;
                plane[local] = allowed[index] ? output[index] : 0;
            }
            memset(merged_charge, 0, (size_t)mz_count * sizeof(float));
            interpolate_merge(mz_axis, merged_charge, data->mz + first, plane,
                              mz_count, point_count);
            for (int mz = 0; mz < mz_count; mz++) {
                const float value = merged_charge[mz];
                mz_values[mz] += value;
                merged_cube[(size_t)mz * charge_count + zi] = value;
            }
        }
        transform_mass_grid_UniChrom(config, merged_cube, mz_axis, 1, mz_count,
                                     charge_count, mass_min, mass_max, mass_axis,
                                     mass_count, mass_values);
        char group[1024], path[1024];
        snprintf(group, sizeof(group), "/ms_dataset/%d", scan);
        snprintf(path, sizeof(path), "%s/mass_data", group);
        mh5writefile2d(config.file_id, path, mass_count, mass_axis, mass_values);
        write_attr_int(config.file_id, group, "length_mz", point_count);
        write_attr_int(config.file_id, group, "length_mass", mass_count);
        write_attr_float(config.file_id, group, "mzsig", config.mzsig);
        write_attr_float(config.file_id, group, "dtsig", config.dtsig);
    }
    for (int scan = 0; scan < data->scan_count; scan++) {
        for (int mz = 0; mz < mz_count; mz++) {
            mz_sum[mz] += mz_grid[(size_t)scan * mz_count + mz];
        }
        for (int mass = 0; mass < mass_count; mass++) {
            mass_sum[mass] += mass_grid[(size_t)scan * mass_count + mass];
        }
    }
    scale_sum_UC(mz_sum, mz_count, mz_grid,
                 data->scan_count * mz_count, config.datanorm);
    scale_sum_UC(mass_sum, mass_count, mass_grid,
                 data->scan_count * mass_count, config.datanorm);
    mh5writefile1d(config.file_id, "/ms_dataset/mz_grid",
                   data->scan_count * mz_count, mz_grid);
    mh5writefile1d(config.file_id, "/ms_dataset/mz_axis", mz_count, mz_axis);
    mh5writefile1d(config.file_id, "/ms_dataset/mz_sum", mz_count, mz_sum);
    mh5writefile1d(config.file_id, "/ms_dataset/mass_grid",
                   data->scan_count * mass_count, mass_grid);
    mh5writefile1d(config.file_id, "/ms_dataset/mass_axis", mass_count, mass_axis);
    mh5writefile1d(config.file_id, "/ms_dataset/mass_sum", mass_count, mass_sum);
    set_got_grids(config.file_id);
    free(output); free(plane); free(convolved);
    free(mz_axis); free(mz_grid); free(mz_sum);
    free(mass_axis); free(mass_grid); free(mass_sum);
    free(merged_cube); free(merged_charge);
    return 1;
}


int run_chromatogram_nonlinear(int argc, char *argv[], Config config)
{
    const clock_t start = clock();
    int result = 1;
    UCRaggedData data;
    UCDirectOperator op;
    memset(&data, 0, sizeof(data));
    memset(&op, 0, sizeof(op));
    config.file_id = H5Fopen(argv[1], H5F_ACC_RDWR, H5P_DEFAULT);
    if (config.file_id < 0) { return 2; }
    const int zero_padding = int_attr(config.file_id, "/config", "unichromzeropad", 0);
    const int scan_count = int_attr(config.file_id, "/ms_dataset", "num", 0);
    if (scan_count < 1 || config.numz < 1 || config.massbins <= 0 ||
        (int64_t)scan_count * scan_count > INT_MAX ||
        !(config.dtsig > 0) || !isfinite(config.dtsig)) {
        fprintf(stderr, "Invalid nonlinear UniChrom dimensions, dtsig, or mass bin size\n");
        goto cleanup;
    }
    printf("UniChrom: processing %d nonlinear spectra\n", scan_count);
    for (int scan = 0; scan < scan_count; scan++) {
        config.metamode = scan;
        config.silent = 1;
        process_data(argc, argv, config);
    }
    if (!read_processed_UC(config, &data) ||
        !build_direct_operator_UC(&data, config, zero_padding, &op)) {
        fprintf(stderr, "Unable to initialize nonlinear UniChrom response\n");
        goto cleanup;
    }

    if ((int64_t)data.point_count * config.numz > INT_MAX) {
        fprintf(stderr, "Nonlinear UniChrom latent cube is too large\n");
        goto cleanup;
    }
    const int cube_length = data.point_count * config.numz;
    float *blur = calloc((size_t)cube_length, sizeof(float));
    float *scratch = calloc((size_t)cube_length, sizeof(float));
    float *oldblur = calloc((size_t)cube_length, sizeof(float));
    float *projection = calloc((size_t)data.point_count, sizeof(float));
    float *response = calloc((size_t)data.point_count, sizeof(float));
    float *ratio = calloc((size_t)data.point_count, sizeof(float));
    float *correction = calloc((size_t)data.point_count, sizeof(float));
    float *smoothing_sums = config.psig >= 1 ?
        calloc((size_t)data.scan_count * config.numz, sizeof(float)) : NULL;
    char *allowed = calloc((size_t)cube_length, sizeof(char));
    int *zup = config.zsig != 0 ? malloc((size_t)cube_length * sizeof(int)) : NULL;
    int *zdown = config.zsig != 0 ? malloc((size_t)cube_length * sizeof(int)) : NULL;
    int *mup = config.msig != 0 ? malloc((size_t)cube_length * sizeof(int)) : NULL;
    int *mdown = config.msig != 0 ? malloc((size_t)cube_length * sizeof(int)) : NULL;
    int *nztab = config.suppression_harmonic > 0 ?
        malloc((size_t)config.numz * sizeof(int)) : NULL;
    if (blur == NULL || scratch == NULL || oldblur == NULL || projection == NULL ||
        response == NULL || ratio == NULL || correction == NULL || allowed == NULL ||
        (config.psig >= 1 && smoothing_sums == NULL) ||
        (config.zsig != 0 && (zup == NULL || zdown == NULL)) ||
        (config.msig != 0 && (mup == NULL || mdown == NULL)) ||
        (config.suppression_harmonic > 0 && nztab == NULL)) {
        fprintf(stderr, "Unable to allocate nonlinear UniChrom iteration arrays\n");
        free(blur); free(scratch); free(oldblur); free(projection); free(response);
        free(ratio); free(correction); free(smoothing_sums); free(allowed);
        free(zup); free(zdown); free(mup); free(mdown); free(nztab);
        goto cleanup;
    }

    for (int point = 0; point < data.point_count; point++) {
        for (int zi = 0; zi < config.numz; zi++) {
            const int index = point * config.numz + zi;
            const int charge = config.startz + zi;
            const float mass = calcmass(data.mz[point], charge, config.adductmass);
            const float native_limit = nativecharge(mass, 0);
            allowed[index] = data.observed[point] > config.intthresh &&
                mass > config.masslb && mass < config.massub &&
                charge < native_limit + config.nativezub &&
                charge > native_limit + config.nativezlb;
            if (allowed[index]) {
                blur[index] = data.observed[point] / (config.numz + 2.0f);
            }
        }
    }
    memcpy(oldblur, blur, (size_t)cube_length * sizeof(float));
    make_regularizer_indexes_UC(&data, config, zup, zdown, mup, mdown);
    if (nztab != NULL) {
        for (int zi = 0; zi < config.numz; zi++) { nztab[zi] = config.startz + zi; }
    }
    const float data_max = Max(data.observed, data.point_count);
    const float beta_factor = data_max > 1 ? data_max : 1;
    double convergence = 0;
    int convergence_seen = 0;
    printf("UniChrom: iterating with ragged direct convolution.");
    for (int iteration = 0; iteration < abs(config.numit); iteration++) {
        if (config.beta > 0 && iteration > 0) {
            softargmax(blur, data.point_count, config.numz, config.beta / beta_factor);
        }
        if (config.psig >= 1 && iteration > 0) {
            #pragma omp parallel for schedule(static) if(data.scan_count > 1)
            for (int scan = 0; scan < data.scan_count; scan++) {
                const int first = data.scan_offsets[scan];
                const int count = data.scan_offsets[scan + 1] - first;
                point_smoothing_to_scratch(blur + (size_t)first * config.numz,
                    scratch + (size_t)first * config.numz,
                    smoothing_sums + (size_t)scan * config.numz,
                    allowed + (size_t)first * config.numz,
                    count, config.numz, abs((int)config.psig));
            }
            float *swap = blur; blur = scratch; scratch = swap;
        }
        if (iteration > config.suppression_startit &&
            (config.suppression_satellite > 0 || config.suppression_harmonic > 0 ||
             config.suppression_topn > 0 || config.suppression_topx > 0)) {
            apply_suppressions(blur, scratch, data.point_count, config.numz,
                               config.suppression_satellite, config.suppression_harmonic,
                               nztab, config.suppression_topn, config.suppression_topx,
                               config.suppression_percent);
        }
        if (config.zsig != 0) {
            blur_it_UCCD(blur, scratch, zup, zdown, cube_length, config.zsig * data_max);
        }
        if (config.msig != 0) {
            blur_it_UCCD(blur, scratch, mup, mdown, cube_length, config.msig * data_max);
        }
        clear_cube_UC(blur, allowed, cube_length);
        #pragma omp parallel for schedule(static) if(data.point_count >= UC_DIRECT_OMP_MIN_LENGTH)
        for (int point = 0; point < data.point_count; point++) {
            float value = 0;
            for (int zi = 0; zi < config.numz; zi++) {
                value += blur[(size_t)point * config.numz + zi];
            }
            projection[point] = value;
        }
        forward_UC(&op, projection, response);
        #pragma omp parallel for schedule(static) if(data.point_count >= UC_DIRECT_OMP_MIN_LENGTH)
        for (int point = 0; point < data.point_count; point++) {
            ratio[point] = response[point] > 0 ? data.observed[point] / response[point] : 0;
        }
        adjoint_UC(&op, ratio, correction);
        #pragma omp parallel for schedule(static) if(data.point_count >= UC_DIRECT_OMP_MIN_LENGTH)
        for (int point = 0; point < data.point_count; point++) {
            const float factor = correction[point] / op.sensitivity[point];
            for (int zi = 0; zi < config.numz; zi++) {
                const int index = point * config.numz + zi;
                const float value = blur[index] * factor;
                blur[index] = allowed[index] && isfinite(value) && value >= 0 ? value : 0;
            }
        }
        if (abs(config.numit) < 10 || iteration == 1 || iteration % 10 == 0 ||
            iteration >= 9 * abs(config.numit) / 10) {
            double difference = 0, total = 0;
            #pragma omp parallel for reduction(+:difference,total) schedule(static) if(cube_length >= UC_DIRECT_OMP_MIN_LENGTH)
            for (int i = 0; i < cube_length; i++) {
                const double delta = blur[i] - oldblur[i];
                difference += delta * delta;
                total += blur[i];
            }
            convergence = total > 0 ? difference / total : INFINITY;
            if (convergence < .000001) {
                if (convergence_seen && config.numit > 0) {
                    printf(" converged in %d iterations", iteration + 1);
                    break;
                }
                convergence_seen = 1;
            } else { convergence_seen = 0; }
            memcpy(oldblur, blur, (size_t)cube_length * sizeof(float));
        }
    }
    printf(" done (metric %.6g)\n", convergence);
    if (write_outputs_UC(config, &data, &op, blur, allowed)) { result = 0; }
    free(blur); free(scratch); free(oldblur); free(projection); free(response);
    free(ratio); free(correction); free(smoothing_sums); free(allowed);
    free(zup); free(zdown); free(mup); free(mdown); free(nztab);

cleanup:
    free_direct_operator_UC(&op);
    free_ragged_UC(&data);
    H5Fclose(config.file_id);
    printf("Nonlinear UniChrom finished in %.3f s\n",
           (float)(clock() - start) / CLOCKS_PER_SEC);
    return result;
}
