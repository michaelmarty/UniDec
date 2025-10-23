# python
from pathlib import Path
import re
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# file paths
chrom_path = Path(r"C:\Users\bht442\Desktop\20251021\20251021_BHT_0.1mgmL_carbonicanhydrase_6inj_0.1gas_50inj_zoom_CDMS_C4_1_2025-10-21-09-53-49_unidecfiles\TIC_chrom.txt")
decon_path = Path(r"C:\Users\bht442\Desktop\20251021\20251021_BHT_0.1mgmL_carbonicanhydrase_6inj_0.1gas_50inj_zoom_CDMS_C4_1_2025-10-21-09-53-49_unidecfiles\TIC_decon.txt")

num_re = re.compile(r'[-+]?\d*\.\d+(?:[eE][-+]?\d+)?|[-+]?\d+(?:[eE][-+]?\d+)?')

def read_two_columns(p: Path):
    t = []
    y = []
    with p.open('r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            nums = num_re.findall(line)
            if len(nums) >= 2:
                t.append(float(nums[0]))
                y.append(float(nums[1]))
    if not t:
        raise SystemExit(f"No numeric data found in {p}")
    return np.array(t), np.array(y)

# read data
time1, inten1 = read_two_columns(chrom_path)
time2, inten2 = read_two_columns(decon_path)

# Apply Mass Distribution plot styling
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 28
# make major tick marks bigger and a bit thicker (applies globally)
plt.rcParams['xtick.major.size'] = 12
plt.rcParams['xtick.major.width'] = 2
plt.rcParams['ytick.major.size'] = 12
plt.rcParams['ytick.major.width'] = 2
# move tick labels a bit further away from the tick marks
plt.rcParams['xtick.major.pad'] = 8
plt.rcParams['ytick.major.pad'] = 8
# increase default spacing between axis labels and axes, and between title and axes
plt.rcParams['axes.labelpad'] = 18
plt.rcParams['axes.titlepad'] = 18
# ensure axes edge (spine) width matches major tick width
plt.rcParams['axes.linewidth'] = plt.rcParams['xtick.major.width']

plt.figure(figsize=(16, 9))

# choose a more contrasting but still muted color scheme (colorblind-friendly)
raw_color = '#001f4d'      # changed to a darker blue (deep navy)
demux_color = '#6A4C93'    # muted purple

# Use colorblind-friendly, high-contrast colors for comparison
plt.plot(time1, inten1, '-', color=raw_color, lw=4, label='Raw TIC')
plt.plot(time2, inten2, '-', color=demux_color, lw=4, label='Demultiplexed TIC')

plt.xlabel('Time (min)', labelpad=18)
plt.ylabel('Intensity (# of ions)', labelpad=18)   # changed from 'Ion Count' to 'Intensity (# of ions)'
plt.title('Total Ion Chromatogram', pad=18)

ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
# make visible spines match tick thickness
ax.spines['left'].set_linewidth(plt.rcParams['xtick.major.width'])
ax.spines['bottom'].set_linewidth(plt.rcParams['xtick.major.width'])
# ensure axis applies the larger major ticks (no minor ticks)
ax.tick_params(which='major', length=12, width=2, direction='out', pad=plt.rcParams['xtick.major.pad'])

# compute y-limits from raw data with a small padding
ymin = min(np.min(inten1), np.min(inten2))
ymax = max(np.max(inten1), np.max(inten2))

# if all intensities are non-negative, set bottom to 0
if ymin >= 0:
    ymin = 0.0

# add padding and force top at 60; ensure 60 is a y-tick
if np.isfinite(ymax) and np.isfinite(ymin):
    span = ymax - ymin
    pad = span * 0.02 if span > 0 else abs(ymax) * 0.1 + 1e-9
    top = 60.0                                  # force top to 60 so the 60 tick is visible
    bottom = ymin - pad
    # guard: ensure bottom is below top
    if bottom >= top:
        bottom = max(0.0, top - abs(pad) - 1.0)
    ax.set_ylim(bottom, top)
    yt = list(ax.get_yticks())
    if 60.0 not in yt:
        yt.append(60.0)
        yt = sorted(set(yt))
        ax.set_yticks(yt)

# x-limits spanning both traces; ensure 0 is shown on the unzoomed plot
xmin = min(np.min(time1), np.min(time2))
xmax = max(np.max(time1), np.max(time2))
if xmin > 0:
    ax.set_xlim(0.0, xmax)
else:
    ax.set_xlim(xmin, xmax)

# ensure 0 is present as a tick on the unzoomed axis
xt = list(ax.get_xticks())
if 0.0 not in xt:
    xt.append(0.0)
    xt = sorted(set(xt))
    ax.set_xticks(xt)

formatter = ScalarFormatter(useOffset=False, useMathText=False)
formatter.set_scientific(False)
ax.xaxis.set_major_formatter(formatter)

# add legend without a surrounding box
ax.legend(frameon=False)

plt.tight_layout()

# save full-range plot, then close it so we can create a separate zoom figure
out_png = chrom_path.with_name(f"{chrom_path.stem}_with_decon.png")
plt.savefig(out_png, dpi=300, transparent=True)
plt.close()
print(f"Saved plot to: {out_png}")

# Zoom window (minutes)
zoom_min, zoom_max = 10.0, 14.0

# build masks for each trace inside the zoom range
mask1 = (time1 >= zoom_min) & (time1 <= zoom_max)
mask2 = (time2 >= zoom_min) & (time2 <= zoom_max)

if not (mask1.any() or mask2.any()):
    print(f"No data in zoom range {zoom_min}-{zoom_max} min; skipping zoom plot.")
else:
    plt.figure(figsize=(16, 9))
    plt.plot(time1[mask1], inten1[mask1], '-', color=raw_color, lw=4, label='Raw TIC')
    plt.plot(time2[mask2], inten2[mask2], '-', color=demux_color, lw=4, label='Demultiplexed TIC')

    plt.xlabel('Time (min)', labelpad=18)
    plt.ylabel('Intensity (# of ions)', labelpad=18)   # changed from 'Ion Count' to 'Intensity (# of ions)'
    plt.title(f'Total Ion Chromatogram (Zoom {zoom_min:.0f}-{zoom_max:.0f} min)', pad=18)

    axz = plt.gca()
    axz.spines['top'].set_visible(False)
    axz.spines['right'].set_visible(False)
    # make visible spines match tick thickness on zoom axis
    axz.spines['left'].set_linewidth(plt.rcParams['xtick.major.width'])
    axz.spines['bottom'].set_linewidth(plt.rcParams['xtick.major.width'])

    # y-limits for zoom calculated from masked data
    y_vals = np.concatenate([inten1[mask1], inten2[mask2]]) if (mask1.any() or mask2.any()) else np.array([0.0])
    ymin_z = np.min(y_vals)
    ymax_z = np.max(y_vals)
    if ymin_z >= 0:
        ymin_z = 0.0
    if np.isfinite(ymin_z) and np.isfinite(ymax_z):
        span = ymax_z - ymin_z
        pad = span * 0.02 if span > 0 else abs(ymax_z) * 0.1 + 1e-9
        top_z = 60.0                             # force top to 60 for zoom so the 60 tick is visible
        bottom_z = ymin_z - pad
        if bottom_z >= top_z:
            bottom_z = max(0.0, top_z - abs(pad) - 1.0)
        axz.set_ylim(bottom_z, top_z)
        ytz = list(axz.get_yticks())
        if 60.0 not in ytz:
            ytz.append(60.0)
            ytz = sorted(set(ytz))
            axz.set_yticks(ytz)

    axz.set_xlim(zoom_min, zoom_max)
    axz.xaxis.set_major_formatter(formatter)
    # apply the same major-tick styling to the zoomed axis (no minor ticks)
    axz.tick_params(which='major', length=12, width=2, direction='out', pad=plt.rcParams['xtick.major.pad'])

    # add legend to zoom plot without a surrounding box
    axz.legend(frameon=False)

    plt.tight_layout()

    out_png_zoom = chrom_path.with_name(f"{chrom_path.stem}_with_decon_zoom_{int(zoom_min)}-{int(zoom_max)}min.png")
    plt.savefig(out_png_zoom, dpi=300, transparent=True)
    plt.close()
    print(f"Saved zoomed plot to: {out_png_zoom}")
