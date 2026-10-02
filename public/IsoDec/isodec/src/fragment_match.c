#include "isodeclib.h"

#include <math.h>
#include <limits.h>
#include <stdlib.h>

static int lower_bound(const double* values, int count, double query)
{
    int left = 0, right = count;
    while (left < right) {
        int middle = left + (right - left) / 2;
        if (values[middle] < query) left = middle + 1;
        else right = middle;
    }
    return left;
}

static int upper_bound(const double* values, int count, double query)
{
    int left = 0, right = count;
    while (left < right) {
        int middle = left + (right - left) / 2;
        if (values[middle] <= query) left = middle + 1;
        else right = middle;
    }
    return left;
}

static void insert_top(double* top, int count, double value)
{
    if (value <= top[0]) return;
    top[0] = value;
    for (int i = 1; i < count && top[i - 1] > top[i]; i++) {
        double temporary = top[i - 1];
        top[i - 1] = top[i];
        top[i] = temporary;
    }
}

void free_fragment_hits(struct FragmentHit* hits)
{
    free(hits);
}

int match_fragment_batch(
    const double* mz, const double* intensity, int spectrum_count,
    const double* masses, const float* envelopes, int fragment_count, int isolen,
    int max_charge, double adduct, double isotope_threshold, double ppm_tolerance,
    int min_peaks, double cosine_threshold, double min_area, int minus_one_as_zero,
    struct FragmentHit** output, int* output_count)
{
    // Retained in the ABI; no preceding isotope is present in these windows.
    (void)minus_one_as_zero;
    if (!output || !output_count || spectrum_count < 0 || fragment_count < 0 ||
        isolen < 1 || isolen > ISODEC_FRAGMENT_ISOTOPES || min_peaks < 1 ||
        min_peaks > isolen || max_charge < 0 || !isfinite(adduct) ||
        !isfinite(isotope_threshold) || isotope_threshold < 0 ||
        !isfinite(ppm_tolerance) || ppm_tolerance < 0 ||
        !isfinite(cosine_threshold) || !isfinite(min_area) ||
        (spectrum_count && (!mz || !intensity)) ||
        (fragment_count && (!masses || !envelopes))) return -1;
    *output = NULL;
    *output_count = 0;
    if (!spectrum_count || !fragment_count) return 0;
    for (int i = 0; i < spectrum_count; i++) {
        if (!isfinite(mz[i]) || !isfinite(intensity[i]) || mz[i] <= 0 ||
            intensity[i] < 0 || (i && mz[i] < mz[i - 1])) return -1;
    }
    if (mz[0] <= adduct && max_charge == 0) return -1;

    struct FragmentHit* hits = NULL;
    int capacity = 0, count = 0;
    for (int fragment = 0; fragment < fragment_count; fragment++) {
        double mass = masses[fragment];
        if (!isfinite(mass) || mass <= 0) continue;
        const float* envelope = envelopes + (size_t)fragment * (size_t)isolen;
        double maximum = 0;
        for (int j = 0; j < isolen; j++) {
            if (!isfinite(envelope[j]) || envelope[j] < 0) goto failure;
            if (envelope[j] > maximum) maximum = envelope[j];
        }
        if (maximum <= 0) continue;
        int n_iso = 0;
        double theory[ISODEC_FRAGMENT_ISOTOPES], neutral[ISODEC_FRAGMENT_ISOTOPES];
        for (int j = 0; j < isolen; j++) {
            if (envelope[j] > maximum * isotope_threshold) {
                theory[n_iso] = envelope[j];
                neutral[n_iso] = mass + j * 1.0033;
                n_iso++;
            }
        }
        if (n_iso < min_peaks || mz[spectrum_count - 1] <= adduct) continue;
        double first_charge = ceil(neutral[0] / (mz[spectrum_count - 1] - adduct));
        if (first_charge >= INT_MAX) goto failure;
        int first_z = (int)first_charge;
        if (first_z < 1) first_z = 1;
        double last_charge = mz[0] > adduct ? floor(neutral[n_iso - 1] / (mz[0] - adduct)) : max_charge;
        if (max_charge && last_charge > max_charge) last_charge = max_charge;
        if (last_charge >= INT_MAX) goto failure;
        int last_z = (int)last_charge;
        if (max_charge && last_z > max_charge) last_z = max_charge;
        for (int z = first_z; z <= last_z; z++) {
            double first_mz = neutral[0] / z + adduct;
            double last_mz = neutral[n_iso - 1] / z + adduct;
            double margin = last_mz * ppm_tolerance * 1e-6;
            int left = lower_bound(mz, spectrum_count, first_mz - margin);
            int right = upper_bound(mz, spectrum_count, last_mz + margin);
            if (right - left < min_peaks) continue;
            double observed[ISODEC_FRAGMENT_ISOTOPES] = {0};
            int matched_iso[ISODEC_FRAGMENT_ISOTOPES], matched_centroid[ISODEC_FRAGMENT_ISOTOPES];
            int n_match = 0, lowest = left, unique = 0;
            double diff = mz[left] * ppm_tolerance * 1e-6;
            for (int k = 0; k < n_iso; k++) {
                double query = neutral[k] / z + adduct;
                double low = query - diff, high = query + diff;
                double best_intensity = 0;
                int best = -1;
                for (int p = lowest; p < right; p++) {
                    if (mz[p] > high) break;
                    if (mz[p] < low) lowest = p + 1;
                    else if (intensity[p] > best_intensity) {
                        best_intensity = intensity[p];
                        best = p;
                    }
                }
                if (best >= 0) {
                    matched_iso[n_match] = k;
                    matched_centroid[n_match] = best - left;
                    observed[k] = intensity[best];
                    int seen = 0;
                    for (int t = 0; t < n_match; t++) if (matched_centroid[t] == best - left) seen = 1;
                    if (!seen) unique++;
                    n_match++;
                }
            }
            if (unique < min_peaks) continue;
            double ab = 0, a2 = 0, b2 = 0, observed_max = 0;
            // This window starts at the first retained isotope. There is no
            // preceding sample to penalize; never wrap around to the last one.
            for (int k = 0; k < n_iso; k++) {
                ab += observed[k] * theory[k];
                a2 += observed[k] * observed[k];
                b2 += theory[k] * theory[k];
                if (observed[k] > observed_max) observed_max = observed[k];
            }
            double score = ab && a2 && b2 ? ab / sqrt(a2 * b2) : 0;
            if (score < cosine_threshold) continue;
            double scale = observed_max / maximum;
            double area = 0, local_area = 0;
            for (int k = 0; k < n_match; k++) area += theory[matched_iso[k]] * scale;
            for (int p = left; p < right; p++) local_area += intensity[p];
            double top_match[ISODEC_FRAGMENT_ISOTOPES] = {0};
            double top_local[ISODEC_FRAGMENT_ISOTOPES] = {0};
            for (int k = 0; k < n_match; k++) insert_top(top_match, min_peaks, intensity[left + matched_centroid[k]]);
            for (int p = left; p < right; p++) insert_top(top_local, min_peaks, intensity[p]);
            int top_equal = 1;
            for (int k = 0; k < min_peaks; k++) if (top_match[k] != top_local[k]) top_equal = 0;
            if ((local_area == 0 || area / local_area <= min_area) && !top_equal) continue;
            if (count == capacity) {
                int next = capacity ? capacity * 2 : 64;
                struct FragmentHit* enlarged = realloc(hits, (size_t)next * sizeof(*hits));
                if (!enlarged) goto failure;
                hits = enlarged;
                capacity = next;
            }
            struct FragmentHit* hit = &hits[count++];
            hit->fragment_index = fragment;
            hit->charge = z;
            hit->left = left;
            hit->right = right;
            hit->isotope_count = n_iso;
            hit->match_count = n_match;
            hit->score = score;
            hit->scale = scale;
            for (int k = 0; k < n_match; k++) {
                hit->isotope_indexes[k] = matched_iso[k];
                hit->centroid_indexes[k] = matched_centroid[k];
            }
        }
    }
    *output = hits;
    *output_count = count;
    return 0;
failure:
    free(hits);
    return -1;
}
