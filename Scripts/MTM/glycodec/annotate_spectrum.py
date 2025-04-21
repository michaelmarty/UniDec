from collections import namedtuple, defaultdict

import glypy
import numpy as np

from matplotlib import pyplot as plt, font_manager, transforms as mtransform

from glypy.io import glycoct
from glypy.composition.composition_transform import derivatize
from glypy.structure.fragment import latex_symbol_map
from glypy import plot

from pyteomics import mgf


def _make_tex_name(self):
    buff = []
    for c in self.name:
        if c in latex_symbol_map:
            buff.append("_{}".format(latex_symbol_map[c]))
        else:
            buff.append(c)
    return ''.join(buff)


FragmentGroup = namedtuple("FragmentGroup", ('name', ))
font_options = font_manager.FontProperties(family='sans serif')

default_ion_series_to_color = {
    'B': 'blue',
    'Y': 'red',
}

water = glypy.Composition("H2O")
water_mass = water.mass


Peak = namedtuple('Peak', ("mz", "intensity", "charge", "neutral_mass"))


def neutral_mass(mz, charge):
    return (mz * charge) - (charge * 1.00727646677)


class PeakSet(object):
    def __init__(self, peaks):
        self.peaks = sorted(peaks, key=lambda x: x.neutral_mass)

    def __getitem__(self, i):
        return self.peaks[i]

    def __len__(self):
        return len(self.peaks)

    def all_peaks_for(self, mass, error_tolerance=2e-5):
        n = len(self.peaks)
        lo = 0
        hi = n
        if n == 0:
            return []
        while hi != lo:
            mid = (hi + lo) // 2
            peak = self[mid]
            error = (peak.neutral_mass - mass) / mass
            if (abs(error) <= error_tolerance) or (hi - 1 == lo):
                i = mid
                while i >= 0:
                    peak = self[i]
                    if abs((peak.neutral_mass - mass) / mass) >= error_tolerance:
                        break
                    i -= 1
                start = i + 1
                i = mid
                while i < n:
                    peak = self[i]
                    if abs((peak.neutral_mass - mass) / mass) >= error_tolerance:
                        break
                    i += 1
                end = i
                return self[start:end]
            elif peak.neutral_mass > mass:
                hi = mid
            else:
                lo = mid
        return []

    @classmethod
    def from_mgf(cls, mgf_scan):
        mz_array = mgf_scan['m/z array']
        intensity_array = mgf_scan['intensity array']
        charge_array = mgf_scan['charge array']

        peaks = []
        for i, mz in enumerate(mz_array):
            peaks.append(
                Peak(mz, intensity_array[i], charge_array[i],
                     neutral_mass(mz, charge_array[i])))
        return cls(peaks)


class GlycanSpectrumMatcher(object):

    def __init__(self, spectrum, target):
        self.target = target
        self.fragment_matches = set()
        self.spectrum = spectrum

    def _match_fragments(self, series, error_tolerance=2e-5, max_cleavages=2, include_neutral_losses=None):
        for fragment in self.target.fragments(series, max_cleavages=max_cleavages):
            peaks = self.spectrum.all_peaks_for(fragment.mass, error_tolerance)
            for peak in peaks:
                self.fragment_matches.add((peak, fragment))
            if include_neutral_losses:
                for loss, label in include_neutral_losses:
                    peaks = self.spectrum.all_peaks_for(fragment.mass - loss.mass, error_tolerance)
                    if peaks:
                        fragment_loss = fragment.copy()
                        fragment_loss.mass -= loss.mass
                        fragment_loss.composition -= loss
                        fragment_loss.name += label
                    for peak in peaks:
                        self.fragment_matches.add((peak, fragment_loss))

    def match(self, error_tolerance=2e-5, *args, **kwargs):
        include_neutral_losses = kwargs.get("include_neutral_losses", False)
        max_cleavages = kwargs.get("max_cleavages", 2)
        if include_neutral_losses and isinstance(include_neutral_losses, (int, bool)):
            include_neutral_losses = [(glypy.Composition("H2O"), "-H2O")]

        self._match_fragments(
            "BY", error_tolerance, max_cleavages=max_cleavages,
            include_neutral_losses=include_neutral_losses)

        return self

    @classmethod
    def evaluate(cls, scan, target, *args, **kwargs):
        return cls(scan, target).match(*args, **kwargs)


class SpectrumAnnotator(object):
    def __init__(self, spectrum_match, ax=None, clip_labels=True):
        if ax is None:
            _, ax = plt.subplots(1)
        self.spectrum_match = spectrum_match
        self.ax = ax
        self.clip_labels = clip_labels
        self.upper = max(
            spectrum_match.spectrum, key=lambda x: x.intensity
        ).intensity * 1.35
        self.peak_labels = []

    def draw_all_peaks(self, color='black', alpha=0.5, **kwargs):
        draw_peaklist(
            self.spectrum_match.spectrum,
            alpha=0.3, color='grey', ax=self.ax, **kwargs)

    def label_peak(self, fragment, peak, fontsize=12, rotation=90, **kw):
        label = "%s" % fragment.name
        if peak.charge > 1:
            label += "$^{%d}$" % peak.charge
        y = peak.intensity
        y = min(y + 100, self.upper * 0.95)

        kw.setdefault("clip_on", self.clip_labels)
        clip_on = kw['clip_on']

        text = self.ax.text(
            peak.mz, y, label, rotation=rotation, va='bottom',
            ha='center', fontsize=fontsize, fontproperties=font_options,
            clip_on=clip_on)
        self.peak_labels.append(text)
        return text

    def format_axes(self):
        draw_peaklist([], self.ax, pretty=True)
        self.ax.set_ylim(0, self.upper)

    def draw_matched_peaks(self, color='red', alpha=0.8, fontsize=12, ion_series_to_color=None, **kwargs):
        if ion_series_to_color is None:
            ion_series_to_color = {}

        by_peak = defaultdict(list)
        for peak, fragment in self.spectrum_match.fragment_matches:
            by_peak[peak].append(fragment)

        for peak, fragments in by_peak.items():
            fragment = fragments[0]
            label = '\n'.join([r"$\overline{%s}$" % _make_tex_name(f) if i != 0 else "$%s$" % _make_tex_name(f)
                               for i, f in enumerate(fragments)])
            try:
                peak_color = ion_series_to_color.get(fragment.series, color)
            except AttributeError:
                peak_color = ion_series_to_color.get(fragment.kind, color)
            draw_peaklist([peak], alpha=alpha, ax=self.ax, color=peak_color)
            self.label_peak(FragmentGroup(label), peak, fontsize=fontsize, **kwargs)

    def draw(self, **kwargs):
        fontsize = kwargs.pop('fontsize', 9)
        rotation = kwargs.pop("rotation", 90)
        clip_labels = kwargs.pop("clip_labels", self.clip_labels)
        self.clip_labels = clip_labels
        ion_series_to_color = kwargs.pop("ion_series_to_color", default_ion_series_to_color)
        self.draw_all_peaks(**kwargs)
        self.draw_matched_peaks(
            fontsize=fontsize, ion_series_to_color=ion_series_to_color,
            rotation=rotation, **kwargs)
        self.format_axes()
        return self


def peaklist_to_vector(peaklist, width=0.000001):
    mzs = []
    intensities = []
    for peak in sorted(peaklist, key=lambda x: x.mz):
        mzs.append(peak.mz - width)
        intensities.append(0.)
        mzs.append(peak.mz)
        intensities.append(peak.intensity)
        mzs.append(peak.mz + width)
        intensities.append(0.)
    return np.array(mzs), np.array(intensities)


def draw_peaklist(peaklist, ax=None, normalize=False, **kwargs):
    pretty = kwargs.pop("pretty", True)
    if ax is None:
        _, ax = plt.subplots(1)
    mz_array, intensity_array = peaklist_to_vector(peaklist)
    if normalize:
        intensity_array = intensity_array / intensity_array.max() * 100.0
    ax.plot(mz_array, intensity_array, **kwargs)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Relative Intensity")
    if pretty:
        if intensity_array.shape[0] > 0:
            set_ylim = intensity_array.min() >= 0
        else:
            set_ylim = True
        _beautify_axes(ax, set_ylim)
    return ax


def _beautify_axes(ax, set_ylim=True):
    ax.axes.spines['right'].set_visible(False)
    ax.axes.spines['top'].set_visible(False)
    ax.yaxis.tick_left()
    ax.xaxis.tick_bottom()
    ax.xaxis.set_ticks_position('none')
    ax.get_xaxis().get_major_formatter().set_useOffset(False)
    if set_ylim:
        ax.set_ylim(0, max(ax.get_ylim()))
    return ax


def main():
    # Load glycan structure and derivatize it
    biant = glycoct.load(open("biantennary.glycoct", "rt"))
    biant.reducing_end = glypy.ReducedEnd("H[2]1H1")
    derivatize(biant, 'methyl')
    print(biant)
    exit()

    # Load spectrum
    scan = next(mgf.read("example.mgf"))
    peaks = PeakSet.from_mgf(scan)

    # Match fragments
    match = GlycanSpectrumMatcher.evaluate(
        peaks, biant,
        include_neutral_losses=[(water, "-H2O"),  # Water Loss
                                (glypy.Composition("CH2"), '-CH2'),  # Methyl Loss
                                (glypy.Composition("CH2H2O"), '-CH4O'),  # Water + Methyl Loss
                                ])

    # Configure Canvas
    fig, ax = plt.subplots(1)
    fig.set_figwidth(12)
    fig.set_figheight(8)
    fig.dpi = 120

    art = SpectrumAnnotator(match, ax=ax)
    art.draw()
    art.format_axes()
    art.ax.set_xlim(0, 2000)
    tax = ax.twinx()
    tax.axis('off')
    # stretch out the structure to leave more spce between branches for annotations
    dt, _ = plot.plot(biant, ax=tax, layout_args={"transform": mtransform.Affine2D().scale(2, 1)})
    dt.transform.scale(2e2, 0.125).translate(1800, 0.55)
    seen = set()
    fragment_matches = sorted(match.fragment_matches, key=lambda peak_fragment: peak_fragment[1].break_count,
                              reverse=True)
    for peak, fragment in fragment_matches:
        if fragment in seen:
            continue
        else:
            seen.add(fragment)
        if fragment.is_internal():
            dt.annotate_fragment(fragment, label=False, lw=2, line_length_scale=0.5, color='navy')
        else:
            dt.annotate_fragment(fragment, label="-" not in fragment.name, lw=2, line_length_scale=0.75,
                                 fontsize=0.33)

    ymax = max(match.spectrum, key=lambda x: x.intensity).intensity
    ticks = np.linspace(0, ymax, 5)
    ax.set_yticks(ticks)
    ax.set_yticklabels("%0.3g%%" % i for i in ticks / ticks.max() * 100.0)
    # plt.savefig("annotated_spectrum.pdf", bbox_inches='tight')
    plt.show()


if __name__ == '__main__':
    main()
