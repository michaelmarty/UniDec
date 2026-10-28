#include "isodeclib.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>

int match_group_centroids_batch(
    const double* mz, const double* intensity, int spectrum_count,
    const int* windows, const double* isotope_mz, const int* isotope_offsets,
    int peak_count, double ppm_tolerance,
    int* match_counts, int* match_indexes, int match_capacity)
{
    if (spectrum_count < 0 || peak_count < 0 || peak_count > INT_MAX / 2 || match_capacity < 0 ||
        !isfinite(ppm_tolerance) || ppm_tolerance < 0 ||
        (spectrum_count && (!mz || !intensity)) ||
        (peak_count && (!windows || !isotope_offsets || !match_counts)) ||
        (match_capacity && !match_indexes)) return -1;
    if (!peak_count) return 0;
    if (isotope_offsets[0] != 0) return -1;
    for (int i = 0; i < spectrum_count; i++) {
        if (!isfinite(mz[i]) || !isfinite(intensity[i]) || mz[i] <= 0 ||
            intensity[i] < 0 || (i && mz[i] < mz[i - 1])) return -1;
    }
    int output_count = 0;
    for (int peak = 0; peak < peak_count; peak++) {
        const int start = windows[2 * peak];
        const int end = windows[2 * peak + 1];
        const int iso_start = isotope_offsets[peak];
        const int iso_end = isotope_offsets[peak + 1];
        if (start < 0 || end > spectrum_count || end <= start ||
            iso_start < 0 || iso_end < iso_start || iso_end > match_capacity ||
            (iso_end > iso_start && !isotope_mz)) return -1;
        match_counts[peak] = 0;
        const double tolerance = mz[start] * ppm_tolerance * 1e-6;
        int lowest = start;
        for (int iso = iso_start; iso < iso_end; iso++) {
            const double query = isotope_mz[iso];
            if (!isfinite(query) || query <= 0) return -1;
            const double low = query - tolerance;
            const double high = query + tolerance;
            double top_intensity = 0;
            int top_index = -1;
            for (int index = lowest; index < end; index++) {
                if (mz[index] > high) break;
                if (mz[index] < low) {
                    lowest = index + 1;
                } else if (intensity[index] > top_intensity) {
                    top_intensity = intensity[index];
                    top_index = index - start;
                }
            }
            if (top_index >= 0) {
                if (output_count >= match_capacity) return -1;
                match_indexes[output_count++] = top_index;
                match_counts[peak]++;
            }
        }
    }
    return 0;
}

struct GroupNode {
    struct MassGroupResult value;
    double lookup_mass;
    int* scans;
    int scan_count;
    int float32_intensity, float32_massdist;
};

int mass_group_abi_version(void)
{
    return ISODEC_MASS_GROUP_ABI_VERSION;
}

static double* copy_doubles(const double* source, int count)
{
    if (count <= 0) return NULL;
    double* result = (double*)malloc((size_t)count * sizeof(double));
    if (result) memcpy(result, source, (size_t)count * sizeof(double));
    return result;
}

static int nearest_scalar(const double* values, int count, double target)
{
    if (count < 2) return 0;
    int lo = 0, hi = count - 1;
    while (hi - lo > 1) {
        int middle = lo + (hi - lo) / 2;
        if (target < values[middle]) hi = middle;
        else if (target > values[middle]) lo = middle;
        else return middle;
    }
    return fabs(target - values[lo]) >= fabs(target - values[hi]) ? hi : lo;
}

static int nearest_pairs(const double* values, int count, double target)
{
    if (count < 2) return 0;
    int lo = 0, hi = count - 1;
    while (hi - lo > 1) {
        int middle = lo + (hi - lo) / 2;
        if (target < values[2 * middle]) hi = middle;
        else if (target > values[2 * middle]) lo = middle;
        else return middle;
    }
    return fabs(target - values[2 * lo]) >= fabs(target - values[2 * hi]) ? hi : lo;
}

static int within_ppm(double theoretical, double experimental, double tolerance)
{
    return theoretical > 0 && fabs((theoretical - experimental) / theoretical) * 1e6 <= tolerance;
}

static double cosine_from_data(const double* centroids, int centroid_count,
                               const double* theory, int theory_count, double ppm)
{
    if (!centroids || !theory || centroid_count < 1 || theory_count < 1) return 0;
    double ab = 0, aa = 0, bb = 0;
    double tolerance = theory[0] * ppm * 1e-6;
    for (int i = 0; i < theory_count; i++) {
        double observed = 0;
        double low = theory[2 * i] - tolerance;
        double high = theory[2 * i] + tolerance;
        for (int j = 0; j < centroid_count; j++) {
            double mass = centroids[2 * j];
            if (mass > high) break;
            if (mass >= low && centroids[2 * j + 1] > observed)
                observed = centroids[2 * j + 1];
        }
        double predicted = theory[2 * i + 1];
        ab += observed * predicted;
        aa += observed * observed;
        bb += predicted * predicted;
    }
    if (ab == 0 || aa == 0 || bb == 0) return 0;
    return ab / sqrt(aa * bb);
}

static int pair_compare(const void* left, const void* right)
{
    const double* a = (const double*)left;
    const double* b = (const double*)right;
    return a[0] < b[0] ? -1 : a[0] > b[0] ? 1 : 0;
}

static int merge_centroids(struct MassGroupResult* group,
                           const struct MassGroupPeak* peak, double ppm)
{
    if (peak->centroid_count > INT_MAX / 2 - group->centroid_count) return -1;
    int count = group->centroid_count + peak->centroid_count;
    if ((size_t)count > SIZE_MAX / (2 * sizeof(double))) return -1;
    double* values = (double*)malloc((size_t)count * 2 * sizeof(double));
    if (!values) return -1;
    memcpy(values, group->decon_centroids,
           (size_t)group->centroid_count * 2 * sizeof(double));
    memcpy(values + 2 * group->centroid_count, peak->decon_centroids,
           (size_t)peak->centroid_count * 2 * sizeof(double));
    qsort(values, (size_t)count, 2 * sizeof(double), pair_compare);
    int output = 0;
    double mass = values[0], intensity = values[1];
    double weighted = mass * intensity;
    for (int i = 1; i < count; i++) {
        double next_mass = values[2 * i];
        double next_intensity = values[2 * i + 1];
        if (fabs(mass - next_mass) / mass * 1e6 <= ppm) {
            weighted += next_mass * next_intensity;
            intensity += next_intensity;
            if (intensity > 0) mass = weighted / intensity;
        } else {
            values[2 * output] = mass;
            values[2 * output + 1] = intensity;
            output++;
            mass = next_mass;
            intensity = next_intensity;
            weighted = mass * intensity;
        }
    }
    values[2 * output] = mass;
    values[2 * output + 1] = intensity;
    free(group->decon_centroids);
    group->decon_centroids = values;
    group->centroid_count = output + 1;
    return 0;
}

static void fit_massdist(struct MassGroupResult* group, double ppm, int float32_massdist)
{
    double maximum = 0, difference = 0, scale = 1;
    for (int i = 0; i < group->massdist_count; i++) {
        double predicted = group->massdist[2 * i + 1];
        if (predicted <= maximum) continue;
        int nearest = nearest_pairs(group->decon_centroids, group->centroid_count,
                                    group->massdist[2 * i]);
        double observed_mass = group->decon_centroids[2 * nearest];
        if (within_ppm(observed_mass, group->massdist[2 * i], ppm)) {
            maximum = predicted;
            difference = observed_mass - group->massdist[2 * i];
            scale = group->decon_centroids[2 * nearest + 1] / predicted;
        }
    }
    int shift = fabs(difference) < ppm * group->massdist[0] / 1e6;
    for (int i = 0; i < group->massdist_count; i++) {
        if (shift) {
            double updated = group->massdist[2 * i] + difference;
            group->massdist[2 * i] = float32_massdist ? (float)updated : updated;
        }
        if (maximum != 0) {
            double updated = group->massdist[2 * i + 1] * scale;
            group->massdist[2 * i + 1] = float32_massdist ? (float)updated : updated;
        }
    }
}

static int matches_group(const struct GroupNode* node,
                         const struct MassGroupPeak* peak, double ppm,
                         int maxshift, double mass_diff, double cosine,
                         int scan_tolerance)
{
    const struct MassGroupResult* group = &node->value;
    if (fabs(peak->monoiso - group->monoiso) > maxshift * mass_diff * 1.1 ||
        fabs((double)peak->scan - node->scans[node->scan_count - 1]) > scan_tolerance)
        return 0;
    int mass_match = 0;
    for (int shift = -maxshift; shift <= maxshift; shift++) {
        if (within_ppm(group->monoiso, peak->monoiso + shift * mass_diff, ppm)) {
            mass_match = 1;
            break;
        }
    }
    return mass_match && cosine_from_data(group->decon_centroids,
                                          group->centroid_count, peak->massdist,
                                          peak->massdist_count, ppm) >= cosine;
}

static struct GroupNode* new_group(const struct MassGroupPeak* peak, int seed)
{
    struct GroupNode* node = (struct GroupNode*)calloc(1, sizeof(*node));
    if (!node) return NULL;
    struct MassGroupResult* group = &node->value;
    group->monoisos = copy_doubles(peak->monoisos, peak->monoisos_count);
    group->massdist = copy_doubles(peak->massdist, 2 * peak->massdist_count);
    group->decon_centroids = copy_doubles(peak->decon_centroids, 2 * peak->centroid_count);
    node->scans = (int*)malloc(sizeof(int));
    if (!group->monoisos || !group->massdist || !group->decon_centroids || !node->scans) {
        free(group->monoisos);
        free(group->massdist);
        free(group->decon_centroids);
        free(node->scans);
        free(node);
        return NULL;
    }
    node->lookup_mass = peak->monoiso;
    node->scans[0] = peak->scan;
    node->scan_count = 1;
    node->float32_intensity = peak->float32_intensity;
    node->float32_massdist = peak->float32_massdist;
    group->monoiso = peak->monoiso;
    group->lookup_mass = peak->monoiso;
    group->apexintensity = peak->peakint;
    group->totalintensity = peak->matchedintensity;
    group->avgmass = peak->avgmass;
    group->minrt = group->maxrt = group->apexrt = peak->rt;
    group->minscan = group->maxscan = group->apexscan = peak->scan;
    group->seed_index = seed;
    group->totalpeaks = 1;
    group->monoisos_count = peak->monoisos_count;
    group->massdist_count = peak->massdist_count;
    group->centroid_count = peak->centroid_count;
    return node;
}

static int merge_peak(struct GroupNode* node, const struct MassGroupPeak* peak,
                      double ppm)
{
    struct MassGroupResult* group = &node->value;
    if (peak->peakint > group->apexintensity) {
        group->apexintensity = peak->peakint;
        group->apexscan = peak->scan;
        group->apexrt = peak->rt;
    }
    if (peak->scan > group->maxscan) {
        group->maxscan = peak->scan;
        group->maxrt = peak->rt;
    }
    if (peak->scan < group->minscan) {
        group->minscan = peak->scan;
        group->minrt = peak->rt;
    }
    if (merge_centroids(group, peak, ppm) != 0) return -1;
    for (int i = 0; i < peak->monoisos_count; i++) {
        double incoming = peak->monoisos[i];
        int nearest = nearest_scalar(group->monoisos, group->monoisos_count, incoming);
        if (within_ppm(incoming, group->monoisos[nearest], ppm)) {
            group->monoisos[nearest] =
                (group->monoisos[nearest] * group->totalintensity +
                 incoming * peak->matchedintensity) /
                (group->totalintensity + peak->matchedintensity);
        } else {
            if (group->monoisos_count >= INT_MAX / 2) return -1;
            double* extended = (double*)realloc(group->monoisos,
                (size_t)(group->monoisos_count + 1) * sizeof(double));
            if (!extended) return -1;
            group->monoisos = extended;
            group->monoisos[group->monoisos_count++] = incoming;
        }
    }
    if (within_ppm(group->monoiso, peak->monoiso, ppm)) {
        group->monoiso = (group->monoiso * group->totalintensity +
                          peak->monoiso * peak->matchedintensity) /
                         (group->totalintensity + peak->matchedintensity);
    } else {
        double new_css = cosine_from_data(group->decon_centroids, group->centroid_count,
                                          peak->massdist, peak->massdist_count, ppm);
        double old_css = cosine_from_data(group->decon_centroids, group->centroid_count,
                                          group->massdist, group->massdist_count, ppm);
        if (new_css > old_css) {
            double* replacement = copy_doubles(peak->massdist, 2 * peak->massdist_count);
            if (!replacement) return -1;
            free(group->massdist);
            group->massdist = replacement;
            group->massdist_count = peak->massdist_count;
            node->float32_massdist = peak->float32_massdist;
            group->monoiso = peak->monoiso;
        }
    }
    fit_massdist(group, ppm, node->float32_massdist);
    int known_scan = 0;
    for (int i = 0; i < node->scan_count; i++) {
        if (node->scans[i] == peak->scan) known_scan = 1;
    }
    if (!known_scan) {
        int* extended = (int*)realloc(node->scans,
                                      (size_t)(node->scan_count + 1) * sizeof(int));
        if (!extended) return -1;
        node->scans = extended;
        node->scans[node->scan_count++] = peak->scan;
    }
    double updated_total = group->totalintensity + peak->matchedintensity;
    group->totalintensity = node->float32_intensity ? (float)updated_total : updated_total;
    group->totalpeaks++;
    return 0;
}

void free_mass_groups(struct MassGroupResult* groups, int group_count, int* group_ids)
{
    if (groups) {
        for (int i = 0; i < group_count; i++) {
            free(groups[i].monoisos);
            free(groups[i].massdist);
            free(groups[i].decon_centroids);
        }
    }
    free(groups);
    free(group_ids);
}

int group_mass_peaks_batch(
    const struct MassGroupPeak* peaks, int peak_count,
    double ppm_tolerance, int maxshift, double mass_diff_c,
    double cosine_threshold, int scan_tolerance, int order,
    struct MassGroupResult** groups, int* group_count, int** group_ids)
{
    if (!groups || !group_count || !group_ids) return -1;
    *groups = NULL;
    *group_count = 0;
    *group_ids = NULL;
    if (peak_count < 0 || peak_count > INT_MAX / 2 ||
        (size_t)peak_count > SIZE_MAX / sizeof(struct GroupNode*) ||
        (peak_count && !peaks) || !isfinite(ppm_tolerance) || ppm_tolerance < 0 ||
        maxshift < 0 || maxshift > 64 ||
        !isfinite(mass_diff_c) || mass_diff_c <= 0 ||
        !isfinite(cosine_threshold) || scan_tolerance < 0 ||
        (order != 0 && order != 1)) return -1;
    if (!peak_count) return 0;
    struct GroupNode** nodes = (struct GroupNode**)calloc((size_t)peak_count, sizeof(*nodes));
    struct GroupNode** assigned = (struct GroupNode**)calloc((size_t)peak_count, sizeof(*assigned));
    int* processing = (int*)malloc((size_t)peak_count * sizeof(int));
    int count = 0, status = -1;
    if (!nodes || !assigned || !processing) goto cleanup;
    for (int i = 0; i < peak_count; i++) {
        const struct MassGroupPeak* peak = &peaks[i];
        if (!isfinite(peak->monoiso) || peak->monoiso <= 0 ||
            !isfinite(peak->matchedintensity) || peak->matchedintensity <= 0 ||
            !isfinite(peak->peakint) || !isfinite(peak->rt) ||
            !isfinite(peak->avgmass) ||
            peak->charge < 1 || peak->monoisos_count < 1 ||
            peak->massdist_count < 1 || peak->centroid_count < 1 ||
            peak->monoisos_count > INT_MAX / 2 ||
            peak->massdist_count > INT_MAX / 2 || peak->centroid_count > INT_MAX / 2 ||
            (size_t)peak->monoisos_count > SIZE_MAX / sizeof(double) ||
            (size_t)peak->massdist_count > SIZE_MAX / (2 * sizeof(double)) ||
            (size_t)peak->centroid_count > SIZE_MAX / (2 * sizeof(double)) ||
            !peak->monoisos || !peak->massdist || !peak->decon_centroids) goto cleanup;
        for (int j = 0; j < peak->monoisos_count; j++) {
            if (!isfinite(peak->monoisos[j]) || peak->monoisos[j] <= 0) goto cleanup;
        }
        for (int j = 0; j < peak->massdist_count; j++) {
            if (!isfinite(peak->massdist[2 * j]) || peak->massdist[2 * j] <= 0 ||
                !isfinite(peak->massdist[2 * j + 1]) || peak->massdist[2 * j + 1] < 0)
                goto cleanup;
        }
        for (int j = 0; j < peak->centroid_count; j++) {
            if (!isfinite(peak->decon_centroids[2 * j]) ||
                peak->decon_centroids[2 * j] <= 0 ||
                !isfinite(peak->decon_centroids[2 * j + 1]) ||
                peak->decon_centroids[2 * j + 1] < 0 ||
                (j && peak->decon_centroids[2 * j] < peak->decon_centroids[2 * (j - 1)]))
                goto cleanup;
        }
        processing[i] = i;
    }
    if (order == 1) {
        // Stable descending order without global comparator state.
        for (int i = 1; i < peak_count; i++) {
            int incoming = processing[i], j = i;
            while (j > 0 &&
                   peaks[processing[j - 1]].matchedintensity < peaks[incoming].matchedintensity) {
                processing[j] = processing[j - 1];
                j--;
            }
            processing[j] = incoming;
        }
    }
    for (int step = 0; step < peak_count; step++) {
        int index = processing[step];
        const struct MassGroupPeak* peak = &peaks[index];
        int nearest = 0, selected = -1;
        double closest_rt = INFINITY;
        if (count) {
            int lo = 0, hi = count - 1;
            while (hi - lo > 1) {
                int middle = lo + (hi - lo) / 2;
                if (peak->monoiso < nodes[middle]->lookup_mass) hi = middle;
                else if (peak->monoiso > nodes[middle]->lookup_mass) lo = middle;
                else { lo = hi = middle; break; }
            }
            nearest = lo == hi ? lo :
                (fabs(peak->monoiso - nodes[lo]->lookup_mass) >=
                 fabs(peak->monoiso - nodes[hi]->lookup_mass) ? hi : lo);
            if (fabs(nodes[nearest]->lookup_mass - peak->monoiso) <= maxshift * 1.25) {
                int upper = nearest, lower = nearest - 1;
                int first = 1;
                while (upper < count || lower >= 0) {
                    if (upper < count) {
                        if (fabs(nodes[upper]->lookup_mass - peak->monoiso) <= maxshift * 1.25) {
                            if (matches_group(nodes[upper], peak, ppm_tolerance,
                                              maxshift, mass_diff_c, cosine_threshold,
                                              scan_tolerance)) {
                                double rt_diff = fabs(peak->rt - nodes[upper]->value.apexrt);
                                if (rt_diff < closest_rt) { selected = upper; closest_rt = rt_diff; }
                            }
                            upper++;
                        } else upper = count;
                    }
                    if (first) { first = 0; continue; }
                    if (lower >= 0) {
                        if (fabs(nodes[lower]->lookup_mass - peak->monoiso) <= maxshift * 1.25) {
                            if (matches_group(nodes[lower], peak, ppm_tolerance,
                                              maxshift, mass_diff_c, cosine_threshold,
                                              scan_tolerance)) {
                                double rt_diff = fabs(peak->rt - nodes[lower]->value.apexrt);
                                if (rt_diff < closest_rt) { selected = lower; closest_rt = rt_diff; }
                            }
                            lower--;
                        } else lower = -1;
                    }
                }
            }
        }
        if (selected >= 0) {
            if (merge_peak(nodes[selected], peak, ppm_tolerance) != 0) goto cleanup;
            assigned[index] = nodes[selected];
        } else {
            struct GroupNode* node = new_group(peak, index);
            if (!node) goto cleanup;
            int insertion = count ? nearest + (peak->monoiso > nodes[nearest]->lookup_mass) : 0;
            memmove(nodes + insertion + 1, nodes + insertion,
                    (size_t)(count - insertion) * sizeof(*nodes));
            nodes[insertion] = node;
            count++;
            assigned[index] = node;
        }
    }
    *groups = (struct MassGroupResult*)calloc((size_t)count, sizeof(**groups));
    *group_ids = (int*)malloc((size_t)peak_count * sizeof(int));
    if (!*groups || !*group_ids) goto cleanup;
    for (int i = 0; i < count; i++) {
        (*groups)[i] = nodes[i]->value;
        nodes[i]->value.monoisos = NULL;
        nodes[i]->value.massdist = NULL;
        nodes[i]->value.decon_centroids = NULL;
        for (int j = 0; j < peak_count; j++) {
            if (assigned[j] == nodes[i]) (*group_ids)[j] = i;
        }
    }
    *group_count = count;
    status = 0;
cleanup:
    for (int i = 0; i < count; i++) {
        if (nodes && nodes[i]) {
            free(nodes[i]->value.monoisos);
            free(nodes[i]->value.massdist);
            free(nodes[i]->value.decon_centroids);
            free(nodes[i]->scans);
            free(nodes[i]);
        }
    }
    free(nodes);
    free(assigned);
    free(processing);
    if (status != 0) {
        free_mass_groups(*groups, *group_count, *group_ids);
        *groups = NULL;
        *group_ids = NULL;
        *group_count = 0;
    }
    return status;
}
