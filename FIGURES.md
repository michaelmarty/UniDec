# Scientific Figure Guide

This document defines the default visual style for scientific figures produced
by the Marty Lab. It is written for both people and AI agents. Its goals are to
make figures consistent across projects, readable at publication size, and
reproducible from source data.

The initial style is based on the figures in
`Scripts/MTM/ArtifactSuppression/`. Those figures use compact journal-sized
layouts, dark-blue data traces, normalized intensity axes, sparse decoration,
direct panel labels, and matched scales for comparisons.

## How to use this guide

The words **MUST**, **SHOULD**, and **MAY** indicate requirement strength:

- **MUST**: follow unless a journal requirement or the scientific meaning of
  the data makes it impossible.
- **SHOULD**: use by default; a clear figure-specific reason may override it.
- **MAY**: optional.

When a journal supplies figure specifications, follow the journal's dimensions,
file format, resolution, and color-mode requirements. Keep the visual choices in
this guide wherever they do not conflict.

## Agent workflow

An AI agent creating or revising a figure MUST:

1. Read this file and the complete script that generates the figure.
2. Identify the comparison being made, the intended panel order, the final
   display width, and any journal requirements.
3. Preserve scientific meaning. Do not smooth, normalize, crop, omit, or
   transform data merely to improve appearance.
4. Define shared styles once and apply them to every panel.
5. Use the same axis limits and normalization for panels intended for direct
   comparison. If different scales are scientifically necessary, label them
   clearly.
6. Save the figure from the script in both editable/vector and review/raster
   forms.
7. Inspect the rendered output at its final physical size. Check labels,
   overlaps, clipping, line visibility, and color interpretation.
8. Report any intentional departure from this guide when handing off the work.

Do not rewrite analysis code solely to standardize a figure. Make the smallest
change that produces a consistent and reproducible result.

## Core visual specification

| Element | Lab default |
|---|---|
| Plotting library | Matplotlib |
| Font | Arial, 11 pt |
| Background | White |
| Primary data color | `darkblue` (`#00008B`) |
| Primary line width | 1.5 pt |
| Panel letters | Uppercase, reading order, upper-left inside each panel |
| Panel description | Upper-right inside each panel |
| Axes | Left and bottom spines only for ordinary 1D plots |
| Relative intensity | Normalize to 0–100; ticks at 0, 50, and 100 |
| Standard two-column width | 6.5 in |
| Standard single-column width | 3.25 in |
| Raster output | PNG at 300 dpi or higher |
| Vector output | PDF |

Arial is the target publication font. If it is unavailable, use a metrically
similar sans-serif font for drafts and disclose the substitution. Do not allow
silent font substitution in a final figure.

## Figure size and layout

- Choose the final physical width before adjusting fonts or spacing. Use 6.5 in
  for a normal two-column figure and 3.25 in for a single-column figure unless
  the journal specifies another width.
- Figures MUST not exceed the journal's maximum width. If a figure is too wide, reduce the number
  of panels or split it into multiple figures or give a warning.
- Set height from the content. As reference points, the Artifact Suppression
  figures use approximately 6.5 x 6.5 in for ten compact spectra, 6.5 x 8 in
  for ten taller spectra, and 6.5 x 8 in for paired 2D and 1D plots.
- Prefer a regular row-by-column grid. The order should be logical but could be either by column or by row. Avoid irregular layouts unless the scientific question requires it.
- Align axes and plotting regions across comparisons. Keep column widths and row
  heights consistent unless the plot types differ.
- Use `GridSpec` or `subplots` with explicit `hspace` and `wspace`. Use
  `layout="constrained"` when it gives stable results; otherwise set margins
  explicitly after checking the rendered file.
- Use whitespace to separate panels. Do not add boxes, background shading, or
  grid lines merely to fill empty space.

## Typography and text

- Use title case for titles and panel descriptions: `Default Settings`, not
  `DEFAULT SETTINGS`.
- Keep text at 11 pt by default. Small secondary annotations MAY use 8–9 pt but
  MUST remain legible at final size.
- Use short, descriptive labels. Put experimental conditions or method names in
  the upper-right of each panel when this avoids a legend.
- Place panel letters at `(0.03, 0.95)` in axes coordinates, aligned left/top.
  Place panel descriptions at `(0.95, 0.95)`, aligned right/top. Adjust only to
  prevent collision with data.
- Letter panels `A`, `B`, `C`, ... in visual reading order: left-to-right, then
  top-to-bottom. Letters MAY be lower case when journal requires it.
- Use a centered title above a column only when it identifies the whole column,
  such as a sample name or molecular mass. Same with rows if grouped by row.
- Use standard scientific notation and units. Put units in parentheses:
  `Mass (kDa)`, `Time (min)`, and `Collision voltage (V)`.
- Format mass-to-charge ratio as italic `m/z` with Matplotlib math text:
  `r'$\mathit{m/z}$'`.
- Put charge-state labels next to the corresponding peaks for raw mass spec data when available, for example `27+`.
  Offset labels so they do not touch the trace or frame.
- Avoid figure-wide titles in publication figures unless the destination calls
  for one. The caption supplies the broader explanation.

## Axes and scales

- Remove the top and right spines from ordinary 1D plots.
- Keep the left and bottom spines, tick marks, and tick labels black.
- Do not show a grid by default.
- For spectra compared by shape, normalize each trace to its maximum and show
  relative intensity from 0 to 100. Use ticks at 0, 50, and 100, displaying the
  middle tick as `%` when the surrounding text or caption makes the convention
  clear.
- A normalization function MUST handle a zero maximum without division by zero.
  It MUST document whether normalization is per trace, per panel, per condition,
  or global.
- Use shared limits for direct comparisons. Do not use independent automatic
  limits when they could exaggerate or hide a difference.
- Show numeric y-axis labels only on the leftmost panel of a repeated row. Show
  numeric x-axis labels only on the bottom panel of a repeated column. Retain
  tick marks where they help alignment.
- Label the x-axis on a non-bottom panel when its quantity differs from the
  panels below it, as in an `m/z` spectrum above mass distributions.
- Use kDa for large molecular masses when it keeps tick labels compact, and
  transform the plotted values explicitly (`mass / 1000`).
- Start quantitative bar axes at zero. For other plot types, choose limits based
  on the scientific question and disclose any axis break.

## Lines, colors, and annotations

Use `darkblue` (`#00008B`) for the primary spectrum, chromatogram, or fitted
trace. Keep its meaning consistent throughout a figure. Use a 1.5 pt line by
default; increase it only when reduction to final size makes the trace hard to
see.

Accent colors are for scientific annotation, not decoration. The Artifact
Suppression figures use red, magenta, orange, and goldenrod to mark distinct
artifact regions or expected masses. These colors do not have universal lab
meanings. A new figure MUST define their meaning in the caption or legend and
use the same mapping in every panel.

- Prefer direct labels over legends when there are only a few unambiguous
  elements.
- If a legend is necessary, keep it frameless and place it where it does not
  cover data.
- Do not rely on color alone. Combine color with a label, marker, line style, or
  position so the annotation survives grayscale printing and color-vision
  differences.
- Use saturated accent colors sparingly. The data should remain visually
  dominant.
- Reference lines SHOULD be thinner and lighter than primary data. The example
  expected-mass markers use dashed goldenrod lines, about 1 pt wide and 50%
  opacity.
- Use horizontal bars or shaded spans to identify ranges. A span is generally
  clearer when the interval itself matters; a short bar is appropriate when the
  annotation should not obscure data.
- Avoid arbitrary rainbow palettes for categorical data.

## Heatmaps, images, and colorbars

- Use a perceptually uniform sequential colormap for nonnegative intensity data
  unless a domain-specific colormap is required. Use a diverging colormap only
  when a meaningful center value exists.
- Set color limits explicitly and keep them identical across heatmaps being
  compared. Relative-intensity heatmaps normally use 0–100.
- Include a colorbar for each distinct scale. Repeated panels with the same scale
  MAY share one colorbar if the association is unmistakable.
- Label colorbar units or relative intensity. Use a small number of meaningful
  ticks, commonly 0, 50, and 100.
- Choose white or black annotation text according to local background contrast.
- Rasterize dense heatmap or contour artists inside a vector PDF when needed to
  keep file size manageable; keep text and axes as vectors.
- Prefer viridis color map. Other colormaps to use are ocean_r and nipy_spectral.

## Data integrity and reproducibility

- Figure scripts MUST run from documented inputs. These should be to locations on the shared data server, "Z://".
- Build paths with `pathlib.Path`, relative to the repository, an explicit input
  directory, or an explicit output directory.
- Separate analysis parameters from appearance settings. Record scientific
  processing choices in named constants or configuration objects.
- Seed stochastic simulations explicitly. The Artifact Suppression simulations
  use `np.random.seed(42)`.
- Do not manually edit plotted values or move data marks after export.
- Do not clip unexplained data, hide inconvenient points, or use normalization
  that changes the apparent comparison.
- Captions or methods MUST state material transformations such as smoothing,
  baseline subtraction, normalization, aggregation, fitting, and excluded data.
- Store the script, environment information, and source-data location needed to
  regenerate the final figure. Do not commit private or licensed source data
  merely to make a script self-contained.

## Export

Every final figure script MUST save:

1. A PDF for publication, editing, and archival use.
2. A PNG at 300 dpi or higher for review and presentation use.

Use a transparent background when requested. Otherwise save on white. Avoid
JPEG for plots because its compression introduces artifacts around lines and
text. Give files stable, descriptive names such as `Figure2.pdf` and
`Figure2.png`. Do not encode a developer name or local path in the filename.

Call `savefig` before `show`, because some backends clear or close a figure after
display. Prefer `WxAgg` or `TkAgg` interactive window for standalone figures.

## Matplotlib starting point

This is a starting point, not a shared abstraction that every script must
import. Keeping a short style block in a figure script is often more
reproducible than depending on an unrelated helper module.

```python
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


OUTPUT_DIR = Path(__file__).resolve().parent / "figures"

plt.rcParams.update({
    "font.family": "Arial",
    "font.size": 11,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.5,
    "pdf.fonttype": 42,
    "savefig.facecolor": "white",
})


def normalize_to_percent(intensity):
    """Return intensity normalized per trace to a maximum of 100."""
    intensity = np.asarray(intensity)
    maximum = np.max(intensity)
    return np.zeros_like(intensity) if maximum == 0 else intensity / maximum * 100


def style_spectrum_axis(axis):
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)
    axis.set_ylim(0, 100)
    axis.set_yticks([0, 50, 100], labels=["0", "%", "100"])


fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.8), layout="constrained")

for panel_index, axis in enumerate(axes.flat):
    # Replace x and intensity with the data for this panel.
    axis.plot(x, normalize_to_percent(intensity), color="darkblue")
    axis.text(0.03, 0.95, chr(ord("A") + panel_index),
              transform=axis.transAxes, ha="left", va="top")
    axis.text(0.95, 0.95, panel_descriptions[panel_index],
              transform=axis.transAxes, ha="right", va="top")
    style_spectrum_axis(axis)

for axis in axes[:, 1]:
    axis.tick_params(axis="y", labelleft=False)
for axis in axes[:-1, :].flat:
    axis.tick_params(axis="x", labelbottom=False)
for axis in axes[-1, :]:
    axis.set_xlabel("Mass (kDa)")

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
for suffix, options in {
    ".pdf": {},
    ".png": {"dpi": 300},
}.items():
    fig.savefig(OUTPUT_DIR / f"Figure1{suffix}", **options)
```

Replace the placeholder data and descriptions, then set explicit limits shared
by comparable panels. If the script is intended for headless automation, set a
noninteractive Matplotlib backend before importing `pyplot` or allow the runtime
environment to select one.

## Final review checklist

Before considering a figure complete, inspect both the PDF and PNG and confirm:

- [ ] The figure answers one clear scientific question or comparison.
- [ ] Panel order follows the visual reading order and matches the caption.
- [ ] Font family and apparent size are consistent across every panel.
- [ ] Text remains readable at the final publication width.
- [ ] Comparable panels use the same limits, units, and normalization.
- [ ] The normalization scope and all material data transformations are stated.
- [ ] Repeated tick labels and legends have been removed without creating
      ambiguity.
- [ ] Panel letters, descriptions, titles, and annotations do not overlap data.
- [ ] Colors have consistent meanings and the figure remains interpretable
      without color alone.
- [ ] Lines, markers, and reference annotations remain visible after reduction.
- [ ] No labels, tick marks, or data are clipped at the canvas edge.
- [ ] PDF text and axes remain vector elements where practical.
- [ ] The PNG is at least 300 dpi at final size.
- [ ] The script has no machine-specific paths and regenerates both outputs.
- [ ] The saved outputs were created before any interactive `show` call.
