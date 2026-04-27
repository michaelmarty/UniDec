import os
import re
import glob
import shutil
import colorsys
from unidec.tools import integrate
import numpy as np
import pandas as pd
import rainbow as rb
import matplotlib.pyplot as plt
from matplotlib import colormaps, rcParams
import matplotlib.colors as mcolors
from matplotlib.colorbar import ColorbarBase

'''======================================================================'''
G_FILENAME = "UIB1B_gdata.txt"
R_FILENAME = "UIB1B_rdata.txt"
NPZ_DIR = "NPZ files"
TXT_DIR = "Text files"
FIG_DIR = "Figures"

'''======================================================================'''
# Global color tuning (HSL Modifiers)
global_LUM = 0.8  # brightness (1 = normal, <1 = darker, >1 = brighter)
global_SAT = 1.4  # saturation (1 = normal, >1 = more vivid)
global_HUE = 0  # hue shift (0–1; cycles colors)

'''======================================================================'''
# Global colormap and or color dictionary associated with each channel related data
color_scheme = {
    "alpha": "winter_r",
    "beta": "cool",
    "gamma": "Wistia",
    "gdata": "winter",
    "rdata": "spring",
    "overlay": "cool",
}


def get_filtered_colors(cmap_name, num_colors=256, colorblind_safe=False, **kwargs):
    # Automatic toggle
    if colorblind_safe:
        cmap_name = "cividis" if "rdata" in cmap_name else "viridis"
        LUM_val = 1.0
        SAT_val = 1.0
    else:
        LUM_val = kwargs.get('LUM_factor', global_LUM)
        SAT_val = kwargs.get('SAT_factor', global_SAT)

    HUE_val = kwargs.get('HUE_shift', global_HUE)
    cmap = colormaps.get_cmap(cmap_name)
    processed_colors = []

    for i in range(num_colors):
        r, g, b, a = cmap(i / max(1, num_colors - 1))
        h, l, s = colorsys.rgb_to_hls(r, g, b)
        h = (h + HUE_val) % 1.0
        l = max(0, min(1, l * LUM_val))
        s = max(0, min(1, s * SAT_val))

        r2, g2, b2 = colorsys.hls_to_rgb(h, l, s)
        processed_colors.append((r2, g2, b2, a))

    return processed_colors

def rc_params_setup():
    rcParams['ps.useafm'] = True
    rcParams['ps.fonttype'] = 42
    rcParams['pdf.fonttype'] = 42
    rcParams['lines.linewidth'] = 1
    rcParams['errorbar.capsize'] = 3
    rcParams['patch.force_edgecolor'] = True
    rcParams['patch.facecolor'] = 'b'
    rcParams['lines.markersize'] = 7
    rcParams['font.size'] = 18
    rcParams['font.weight'] = 500
    rcParams['font.sans-serif'] = "Arial"
    rcParams['axes.titleweight'] = 500
    rcParams['axes.titlesize'] = 22
    rcParams['axes.linewidth'] = 3
    rcParams['xtick.major.width'] = 1
    rcParams['ytick.major.width'] = 1
    rcParams['xtick.labelsize'] = rcParams.get('xtick.labelsize', rcParams['font.size'])
    rcParams['ytick.labelsize'] = rcParams.get('ytick.labelsize', rcParams['font.size'])

    try:
        tick_size = float(rcParams['xtick.labelsize'])
    except Exception:
        tick_size = float(rcParams.get('font.size', 12))

    legend_size = max(1.0, tick_size - 3.0)

    rcParams['legend.fontsize'] = legend_size
    rcParams['legend.title_fontsize'] = legend_size

    globals()['LEGEND_FONT_WEIGHT'] = 'bold'

def extract_all_files(files):
    for path in files:
        try:
            datadir = rb.read(path)
            txt_out = os.path.join(path, TXT_DIR)
            npz_out = os.path.join(path, NPZ_DIR)

            os.makedirs(txt_out, exist_ok=True)
            os.makedirs(npz_out, exist_ok=True)

            for f in os.listdir(path):
                if f.lower().endswith('.ch'):
                    print(f"Extracting {f}...")
                    df = datadir.get_file(f)
                    stacked = np.column_stack((df.xlabels, df.data))
                    name = os.path.splitext(f)[0] + "_rawdata"
                    np.savez_compressed(os.path.join(npz_out, f"{name}.npz"), data=stacked)
                    np.savetxt(os.path.join(txt_out, f"{name}.txt"), X=stacked, delimiter="\t")
        except Exception as e:
            print(f" [EXTRACT ERROR] Failed on {os.path.basename(path)}: {e}")

def nonlinearize(datatop, num_compressed, gap=20):
    if num_compressed == 0:
        return datatop
    else:
        num_compressed = int(num_compressed)
        gap = int(gap)
        out = np.array([np.mean(datatop[index + gap:index + num_compressed - gap], axis=0) for index in
                        range(0, len(datatop), num_compressed)])
        # remove nans (if any)
        out = out[~np.isnan(out).any(axis=1)]
        return out

def calc_lock_index(data, lock_time =0.003):
    times = data[:,0]
    diffs = np.abs(times - lock_time)
    lock_index = np.argmin(diffs)
    print("Lock Index:", lock_index, "at time:", data[lock_index,0])
    return lock_index

def parse_lock_data(data, num_compressed, gap=20):
    start_idx = calc_lock_index(data, lock_time = 0.003)

    parseddata = nonlinearize(data[start_idx:], num_compressed, gap=gap)

    pdata = parseddata[0::2] # Positive data
    ndata = parseddata[1::2] # Negative data
    #invert ndata
    ndata[:, 1] = ndata[:, 1] * -1  # Flip the negative channel to positive values
    rdata = ndata
    gdata = pdata
    rdata[:, 1] = rdata[:, 1] - np.min(rdata[:, 1])  # Subtract the minimum to avoid negative values
    gdata[:, 1] = gdata[:, 1] - np.min(gdata[:, 1])  # Subtract the minimum to avoid negative values


    return rdata, gdata

def calc_percents(data, range1=[2.5, 3.5], range2=[10, 12]):
    """
    Calculates the percentage of signal area in two defined windows.
    Returns: (percent1, percent2)
    """
    try:
        # data is an N x 2 array: [time, intensity]
        s1, _ = integrate(data, range1[0], range1[1])  # Area of peak 1
        s2, _ = integrate(data, range2[0], range2[1])  # Area of peak 2

        total = s1 + s2

        if total == 0:
            return 0.0, 0.0

        percent1 = (s1 / total) * 100
        percent2 = (s2 / total) * 100

        return s1, s2, percent1, percent2
    except Exception as e:
        # Returns 0,0 if integration fails (e.g. empty data)
        # to prevent the dashboard from crashing.
        return 0.0, 0.0

def plot_subtraction_verification(path, original, blank, subtracted, channel_name, interp_blank=None):
    plt.figure(figsize=(12, 5))
    plt.axhline(0, color='black', lw=0.5, alpha=0.3)

    # Raw Data
    plt.plot(original[:, 0], original[:, 1], label="Raw Tracer", alpha=0.3, color='gray')

    # Interpolated Blank (This is the one that was subtracted)
    if interp_blank is not None:
        plt.plot(original[:, 0], interp_blank, label="Interpolated Blank", alpha=0.5, color='tab:orange',
                 linestyle='--')

    # Final Result
    plt.plot(original[:, 0], subtracted[:, 1], label="Subtracted Result", color='tab:green', lw=1.5)

    plt.title(f" {os.path.basename(path)} | {channel_name}")
    plt.xlabel("Time (min)");
    plt.ylabel("Intensity")
    plt.legend(loc='upper right', fontsize='x-small')
    plt.tight_layout()

    save_dir = topdir if topdir else path
    fname = f"{os.path.basename(path)}_Verify_Subtraction_{channel_name}.png" if topdir else f"Verify_Subtraction_{channel_name}.png"
    plt.savefig(os.path.join(save_dir, fname), dpi=300)
    plt.close()

def extract_numeric_key(filename):
    """Extracts XX from t_XX or rep_XX for chronological sorting."""
    match = re.search(r'(?:t|rep)_?(\d+)', filename, re.IGNORECASE)
    return int(match.group(1)) if match else 9999

def map_samples_to_blanks(files):
    """Maps each sample to the nearest blank that occurred BEFORE it."""
    sorted_files = sorted(files, key=lambda x: (extract_numeric_key(os.path.basename(x)), x))
    blanks = [f for f in sorted_files if "blank" in os.path.basename(f).lower()]
    samples = [f for f in sorted_files if "blank" not in os.path.basename(f).lower()]

    if not blanks:
        print("[ERROR] No blanks found. Background subtraction impossible.")
        return None

    mapping = {}
    for s_path in samples:
        prior_blanks = [b for b in blanks if sorted_files.index(b) < sorted_files.index(s_path)]
        mapping[s_path] = prior_blanks[-1] if prior_blanks else blanks[0]
    return mapping

def analyze_UIB_samples(files, topdir, num_compressed=50):
    analysis_results = {}
    summary_data = []

    sample_map = map_samples_to_blanks(files)
    if not sample_map:
        return {}

    print(f"[PROCESS] Interpolating and subtracting blanks based on preceding reference...")

    for current_path, blank_path in sample_map.items():
        base_name = os.path.basename(current_path)

        #Metadata Extraction
        time_match = re.search(r't_?(\d+)', base_name, re.IGNORECASE)
        time_hours = int(time_match.group(1)) if time_match else 0
        rep_match = re.search(r'rep_?(\d+)', base_name, re.IGNORECASE)
        rep_number = int(rep_match.group(1)) if rep_match else 0

        #File Setup
        txt_out_dir = os.path.join(current_path, TXT_DIR)
        os.makedirs(txt_out_dir, exist_ok=True)

        try:
            t_raw_data = np.loadtxt(os.path.join(txt_out_dir, "UIB1B_rawdata.txt"))
            b_raw_data = np.loadtxt(os.path.join(blank_path, TXT_DIR, "UIB1B_rawdata.txt"))
        except Exception as e:
            print(f"[SKIP] Could not load data for {base_name}: {e}")
            continue

        #Filter and Demultiplex
        t_filtered = t_raw_data[t_raw_data[:, 0] > 0.005]
        b_filtered = b_raw_data[b_raw_data[:, 0] > 0.005]

        tr_raw, tg_raw = parse_lock_data(t_filtered, num_compressed, gap=20)
        br, bg = parse_lock_data(b_filtered, num_compressed, gap=20)

        #Shift Verification Plot
        # Two-row subplot: row 1 = normalized green & red, row 2 = raw multiplexed
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

        # Row 1: normalized demultiplexed traces
        ax1.plot(tg_raw[:, 0], tg_raw[:, 1] / np.max(tg_raw[:, 1]), label="Green (Norm)", alpha=0.8)
        ax1.plot(tr_raw[:, 0], tr_raw[:, 1] / np.max(tr_raw[:, 1]), label="Red (Norm)", alpha=0.8)

        ax1.set_ylabel("Normalized Intensity")
        ax1.legend(loc='upper right', fontsize='small')

        # Row 2: raw multiplexed trace
        ax2.plot(tg_raw[:, 0], tg_raw[:, 1], label="Green(raw)", alpha=0.8)
        ax2.plot(tr_raw[:, 0], tr_raw[:, 1], label="Red (raw)", alpha=0.8)
        ax2.plot(t_filtered[:, 0], t_filtered[:, 1], label="Raw Multiplexed", color='gray',alpha=0.1)
        ax2.set_xlabel("Time (min)")
        ax2.set_ylabel("Intensity")
        ax2.legend(loc='upper right', fontsize='small')

        # Title and save
        fig.suptitle(f"Demux Verification: {base_name}\n(Using Blank: {os.path.basename(blank_path)})", fontsize=12)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.savefig(os.path.join(current_path, f"{base_name}_UIB1B_shift_diag.png"))
        plt.close(fig)

        #Interpolation and Subtraction Logic
        def subtract_background(sample_data, blank_data):
            interp_blank = np.interp(sample_data[:, 0], blank_data[:, 0], blank_data[:, 1])
            subbed = sample_data.copy()
            subbed[:, 1] = np.maximum(sample_data[:, 1] - interp_blank, 0)

            return subbed, interp_blank

        tg_sub, bg_interp = subtract_background(tg_raw, bg)
        tr_sub, br_interp = subtract_background(tr_raw, br)

        #Ratio Calculation
        g_r1, g_r2, g_p1, g_p2 = calc_percents(tg_sub)
        r_r1, r_r2, r_p1, r_p2 = calc_percents(tr_sub)

        #Store Data for Summary CSV
        summary_data.append({
            "BaseName": base_name,
            "Time_Point": time_hours,
            "Rep_Number": rep_number, "Green_1_raw":g_r1, "Green_2_raw":g_r2,
            "Green_P1": g_p1, "Green_P2": g_p2, "Red_1_raw":r_r1, "Red_2_raw":r_r2,
            "Red_P1": r_p1, "Red_P2": r_p2,
            "Used_Blank": os.path.basename(blank_path)
        })

        #Verification Plots & Saving Text Files
        plot_subtraction_verification(current_path, tg_raw, bg, tg_sub, "Green", interp_blank=bg_interp)
        plot_subtraction_verification(current_path, tr_raw, br, tr_sub, "Red", interp_blank=br_interp)

        np.savetxt(os.path.join(txt_out_dir, "UIB1B_gdata.txt"), tg_sub, header="time\tintensity")
        np.savetxt(os.path.join(txt_out_dir, "UIB1B_rdata.txt"), tr_sub, header="time\tintensity")

        #Store for in-memory access
        analysis_results[current_path] = {
            'tg_raw':tg_raw, 'tr_raw':tr_raw, 'tg_sub': tg_sub, 'tr_sub': tr_sub,
            'hours': time_hours, 'rep': rep_number,
            'g_p1': g_p1, 'g_p2': g_p2,'r_p1': r_p1, 'r_p2': r_p2
        }

    #Export Master Summary (OUTSIDE the loop)
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df = summary_df.sort_values(by=["Time_Point", "Rep_Number"])
        csv_path = os.path.join(topdir, "Dimension_One_Summary.csv")
        summary_df.to_csv(csv_path, index=False)
        print(f"[SUCCESS] Summary report saved to: {csv_path}")

    return analysis_results

def plot_subtraction_verification(path, original, blank, subtracted, channel_name, interp_blank=None):
    print("[PLOTTING] Generating subtraction verification plot...")
    plt.figure(figsize=(12, 5))
    plt.axhline(0, color='black', lw=0.5, alpha=0.3)

    # Raw Data
    plt.plot(original[:, 0], original[:, 1], label="Raw Tracer", alpha=0.5, color='skyblue')

    # Interpolated Blank (This is the one that was subtracted)
    if interp_blank is not None:
        plt.plot(original[:, 0], interp_blank, label="Interpolated Blank", alpha=0.3, color='magenta')

    # Final Result
    plt.plot(original[:, 0], subtracted[:, 1], label="Subtracted Result", color='deepskyblue', lw=1.5)

    plt.title(f" {os.path.basename(path)} | {channel_name}")
    plt.xlabel("Time (min)");
    plt.ylabel("Intensity")
    plt.legend(loc='upper right', fontsize='x-small')
    plt.tight_layout()

    plt.savefig(os.path.join(path, f"Verify_Subtraction_{channel_name}.png"), dpi=300)
    plt.close()

def plot_UIB_dashboard(analysis_results, topdir, sort_by='reps', kin_xaxis='reps', figsize=(20, 18)):

    print(f"[PLOTTING] Generating UIB Dashboard (X-Axis: {kin_xaxis})...")


    if sort_by == 'reps':
        sort_key = lambda x: (analysis_results[x]['hours'], analysis_results[x].get('rep', 0))
    else:
        sort_key = lambda x: (analysis_results[x].get('rep', 0), analysis_results[x]['hours'])

    sorted_paths = sorted(analysis_results.keys(), key=sort_key)


    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.2])

    ax_raw_g, ax_raw_r = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    ax_sub_g, ax_sub_r = fig.add_subplot(gs[1, 0]), fig.add_subplot(gs[1, 1])
    ax_kin_g, ax_kin_r = fig.add_subplot(gs[2, 0]), fig.add_subplot(gs[2, 1])

    g_colors = get_filtered_colors(color_scheme['gdata'], len(sorted_paths))
    r_colors = get_filtered_colors(color_scheme['rdata'], len(sorted_paths))


    for i, p in enumerate(sorted_paths):
        res = analysis_results[p]
        h_label = f"t={res['hours']}h (rep {res.get('rep', 0)})"

        ax_raw_g.plot(res['tg_raw'][:, 0], res['tg_raw'][:, 1], color=g_colors[i], label= h_label, lw=1.2)
        ax_raw_r.plot(res['tr_raw'][:, 0], res['tr_raw'][:, 1], color=r_colors[i], label=h_label,  lw=1.2)

        def plot_norm(ax, data, color):
            mx = np.amax(data[:, 1]) if data.size > 0 else 0
            if mx > 0: ax.plot(data[:, 0], data[:, 1] / mx, color=color, lw=1.2)

        plot_norm(ax_sub_g, res['tg_sub'], g_colors[i])
        plot_norm(ax_sub_r, res['tr_sub'], r_colors[i])

    # Kinetics Data Preparation
    if kin_xaxis == 'hours':
        x_data = [analysis_results[p]['hours'] for p in sorted_paths]
        x_label = "Experimental Time (h)"
    else:
        x_data = [analysis_results[p].get('rep', 0) for p in sorted_paths]
        x_label = "Replicate Number"

    gp2 = [analysis_results[p]['g_p2'] for p in sorted_paths]
    rp1 = [analysis_results[p]['r_p1'] for p in sorted_paths]

    # Plot Green Kinetics
    ax_kin_g.plot(x_data, gp2, color='lime', ls='--', alpha=0.3)
    ax_kin_g.scatter(x_data, gp2, color='lime', s=120, edgecolors='k', label='Bound (P2) %', zorder=3)

    # Plot Red Kinetics
    ax_kin_r.plot(x_data, rp1, color='magenta', ls='--', alpha=0.3)
    ax_kin_r.scatter(x_data, rp1, color='magenta', s=120, edgecolors='k', label='UnBound (P1) %', zorder=3)

    # --- Formatting ---
    titles = ["Green: Raw Multiplexed", "Red: Raw Multiplexed",
              "Green: Subtracted (Normalized)", "Red: Subtracted (Normalized)",
              f"Green Binding vs {kin_xaxis.capitalize()}", f"Red Binding vs {kin_xaxis.capitalize()}"]

    for ax, title in zip([ax_raw_g, ax_raw_r, ax_sub_g, ax_sub_r, ax_kin_g, ax_kin_r], titles):
        ax.set_title(title, fontweight='bold', fontsize=14)

        # Standard axis labels
    for ax in [ax_raw_g, ax_raw_r, ax_sub_g, ax_sub_r]:
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Intensity")

    plt.subplots_adjust(left=0.06, right=0.88, wspace=0.6, hspace=0.4, top=0.92, bottom=0.07)

    # Dynamic x-axis labels for Kinetics
    for ax in [ax_kin_g, ax_kin_r]:
        ax.set_xlabel(x_label)
        ax.set_ylabel("Percent (%)")
        ax.set_ylim(-5, 105)
        ax.grid(True, which='both', linestyle=':', alpha=0.5)
        # If plotting reps, you might want discrete integer ticks
        if kin_xaxis == 'reps':
            from matplotlib.ticker import MaxNLocator
            ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    # 1. Setup Scales and Colors
    norm = mcolors.Normalize(vmin=0, vmax=len(sorted_paths) - 1)
    cmap_g = mcolors.ListedColormap(g_colors)
    cmap_r = mcolors.ListedColormap(r_colors)

    fig.canvas.draw()

    # 3. Setup Scales
    norm = mcolors.Normalize(vmin=0, vmax=len(sorted_paths) - 1)
    cmap_g = mcolors.ListedColormap(g_colors)
    cmap_r = mcolors.ListedColormap(r_colors)

    # 2. Helper to create aligned colorbar
    def add_aligned_cbar(parent_ax, cmap, title):
        pos = parent_ax.get_position()

        # Define cax using the actual rendered position of the subplot
        # x0=pos.x1 + 0.005 puts it just barely to the right of the plot line
        cax = fig.add_axes([pos.x1 + 0.005, pos.y0, 0.01, pos.height])

        cb = ColorbarBase(cax, cmap=cmap, norm=norm)
        cax.set_title(title, fontsize=10, fontweight='bold', pad=10)

        indices = [0, len(sorted_paths) - 1]
        labels = [f"t={analysis_results[sorted_paths[idx]]['hours']}h" for idx in indices]

        cb.set_ticks(indices)
        cb.set_ticklabels(labels, fontsize=9)
        return cb

    # 5. Attach the 4 bars
    add_aligned_cbar(ax_raw_g, cmap_g, " " )
    add_aligned_cbar(ax_raw_r, cmap_r, " " )
    add_aligned_cbar(ax_sub_g, cmap_g, " " )
    add_aligned_cbar(ax_sub_r, cmap_r, " " )

    # 4. Final Layout Adjustment
    # Increase wspace to prevent the G-bars from hitting the R-plot Y-labels
    plt.savefig(os.path.join(topdir, f"UIB_Dashboard_{kin_xaxis}.png"), dpi=300)
    plt.savefig(os.path.join(topdir, f"UIB_Dashboard_{kin_xaxis}.pdf"), dpi=300)
    plt.show()

def parse_array_to_xy(arr):
    """
    Standardizes input arrays into X and Y components.
    Assumes X is column 0 (Time) and Y is column 1 (Intensity).
    """
    arr = np.asarray(arr)

    # Check if the array is 2D and has at least 2 columns
    if arr.ndim == 2 and arr.shape[1] >= 2:
        x = arr[:, 0]
        y = arr[:, 1]
        return x, y

    # Fallback for 1D arrays or unexpected shapes
    print(f" [!] Warning: Unexpected array shape {arr.shape}. Cannot parse to XY.")
    return None, None

def convert_acq_to_utf8(files):
    for d in files:
        for root, _, fnames in os.walk(d):
            for fname in fnames:
                if fname.lower() in ("acq.txt", "acqmethod.txt"):
                    src = os.path.join(root, fname)
                    dst = os.path.join(root, "acq_utf8.txt")
                    if os.path.exists(dst): continue
                    with open(src, "rb") as f:
                        raw = f.read()
                    text = raw.decode("utf-16-le" if b"\x00" in raw else "utf-8", errors="ignore")
                    with open(dst, "w", encoding="utf-8") as f:
                        f.write(re.sub(r"[\x00\xa0]", "", text))

def parse_instrument_method(files):
    name_re = re.compile(r"^[^\S\r\n]*[A-Za-z].*\(G\d{4}[A-Z]\)\s*$")
    for d in files:
        acq = os.path.join(d, "acq_utf8.txt")
        if not os.path.exists(acq): continue
        with open(acq, "r", encoding="utf-8") as f:
            lines = f.read().splitlines()[7:]
        name_idx = [(i, " ".join(line.strip().split())) for i, line in enumerate(lines) if name_re.match(line)]
        out_dir = os.path.join(d, "Method Timetables")
        os.makedirs(out_dir, exist_ok=True)
        for k, (start, mod_name) in enumerate(name_idx):
            end = name_idx[k + 1][0] if k + 1 < len(name_idx) else len(lines)
            body = "\n".join(lines[start + 1:end])
            m = re.search(r"(?i)\bTimetable\b", body)
            if not m: continue
            after = body[m.end():]
            cut_pts = [m.start() for m in [name_re.search(after), re.search(r"(?im)^\s*Instrument\s+Curves\b", after),
                                           re.search(r"\n=+\s*\n", after)] if m]
            tt = after[:min(cut_pts)].strip() if cut_pts else after.strip()
            if len([l for l in tt.splitlines() if l.strip()]) >= 3:
                safe = re.sub(r"[^\w\s\(\)\-\.]", "", mod_name).strip()
                with open(os.path.join(out_dir, f"{safe}.txt"), "w", encoding="utf-8") as f: f.write(tt)

def plot_instrument_timetables(files):
    COLOR_SCHEMES = color_scheme
    print("[PLOTTING] Generating Instrument Method Timetable Visualizations...")

    def ffill(vals):
        last = None
        out = []
        for v in vals:
            if v is None or v == "---":
                out.append(last)
            else:
                try:
                    last = float(v)
                    out.append(last)
                except:
                    out.append(last)
        return out

    for d in files:
        tt_dir = os.path.join(d, "Method Timetables")
        if not os.path.isdir(tt_dir): continue

        files_map = {k: next((f for f in os.listdir(tt_dir) if k in f.lower()), None)
                     for k in ["alpha", "beta", "gamma", "g7116a", "g1316a"]}

        if not any(files_map.values()): continue

        fig, axs = plt.subplots(4, 1, figsize=(10, 12))
        fig.subplots_adjust(hspace=0.4)

        valve_times = []

        def parse_quat(fn):
            with open(os.path.join(tt_dir, fn)) as f:
                lines = [l.strip() for l in f if l.strip()]
            start = next((i for i, l in enumerate(lines) if re.match(r"^-{3,}$", l)), None)
            data = lines[start + 1:] if start else lines
            t, A, B, C, D, F = [], [], [], [], [], []
            for ln in data:
                p = re.split(r"\s+", ln)
                tm = 0.0 if p[0].lower().startswith("start") else float(p[0])
                v = [None if x == "---" else x for x in p[1:6]]
                if len(v) < 5: continue
                t.append(tm);
                A.append(v[0]);
                B.append(v[1]);
                C.append(v[2]);
                D.append(v[3]);
                F.append(v[4])
            return np.array(t), ffill(A), ffill(B), ffill(C), ffill(D), ffill(F)

        # 1 & 2: Alpha and Beta Pumps (Independent Smart Limits)
        for i, key in enumerate(["alpha", "beta"]):
            if files_map[key]:
                t, A, B, C, D, F = parse_quat(files_map[key])
                cmap = plt.get_cmap(COLOR_SCHEMES[key])
                colors = [cmap(x) for x in np.linspace(0.1, 0.9, 4)]
                for j in range(len(t) - 1):
                    dom = np.argmax([float(A[j] or 0), float(B[j] or 0), float(C[j] or 0), float(D[j] or 0)])
                    axs[i].plot(t[j:j + 2], F[j:j + 2], color=colors[dom], lw=3)

                # Independent Smart Limit
                if len(t) > 0:
                    axs[i].set_xlim(0, t.max() + 0.5)
                axs[i].set_title(f"{key.capitalize()} Pump")
        # 2. Extract G1316A Valve Switches (BEFORE Gamma plotting)
        if files_map["g1316a"]:
            with open(os.path.join(tt_dir, files_map["g1316a"])) as f:
                for line in f:
                    # RegEx to find a float at the very start, even if text follows immediately
                    # Captures '12.50' from '12.50Left'
                    match = re.match(r"^(\d+\.?\d*)", line.strip())
                    if match:
                        valve_times.append(float(match.group(1)))

        if files_map["gamma"]:
            with open(os.path.join(tt_dir, files_map["gamma"])) as f:
                lines = [l.strip() for l in f if l.strip()]
            start = next((i for i, l in enumerate(lines) if re.match(r"^-{3,}$", l)), None)
            t_g, B_g, F_g = [], [], []
            for ln in lines[start + 1:]:
                p = re.split(r"\s+", ln)
                t_g.append(0.0 if p[0].lower().startswith("start") else float(p[0]))
                B_g.append(p[2]);
                F_g.append(p[3])

            t_g = np.array(t_g)
            # --- Left Y-Axis: Flow ---
            ax_left = axs[2]
            line_flow = ax_left.plot(t_g, ffill(F_g), 'k', lw=1.5, label="Flow", zorder=3)
            ax_left.set_ylabel("Flow (mL/min)", color='black', fontweight='bold')
            ax_left.tick_params(axis='y', labelcolor='black')

            # --- Right Y-Axis: %B ---
            ax_right = ax_left.twinx()
            line_grad = ax_right.plot(t_g, ffill(B_g), color='tab:orange', ls='--', label="%B", zorder=3)
            ax_right.set_ylabel("Solvent B (%)", color='tab:orange', fontweight='bold')
            ax_right.tick_params(axis='y', labelcolor='tab:orange')
            # Match the spine color to the data
            ax_right.spines['right'].set_color('tab:orange')

            # --- Valve Overlays (Vertical Lines) ---
            for v_time in valve_times:
                ax_left.axvline(x=v_time, color='deepskyblue', linestyle=':', alpha=1, lw=1.5, zorder=1)
                # Position text slightly above the top of the plot
                ax_left.text(v_time, ax_left.get_ylim()[1] * 1.02, f'{v_time}',
                             color='red', fontsize=8, ha='center', weight='bold')

            # Dynamic Scaling
            if len(t_g) > 0:
                max_x = max(t_g.max(), (max(valve_times) if valve_times else 0))
                ax_left.set_xlim(0, max_x + 1)

            ax_left.set_title("Gamma Pump (Flow & Gradient Profiles)", fontweight='bold')

            # Unified Legend for both axes
            lines = line_flow + line_grad
            labels = [l.get_label() for l in lines]
            ax_left.legend(lines, labels, loc='upper left', fontsize='small', framealpha=0.8)

            # 4. Column Compartments (Visualizing the dots)
        axs[3].set_title("Column Compartment Events")
        if valve_times:
            axs[3].scatter(valve_times, [1] * len(valve_times), color='tab:orange', s=60, label="G1316A Switch")
            axs[3].set_xlim(0, max(valve_times) + 1)
            axs[3].set_yticks([])  # Hide Y axis as it's just a timeline
            axs[3].legend(fontsize=8)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        plt.savefig(os.path.join(d, "instrument_method.png"), dpi=200)
        plt.close()

def plot_individual_samples(files, figsize=(12, 10), save=True):
    """
    Standardizes plotting:
    - No normalization for UIB subplot.
    - Blank mapping/subtraction for non-UIB (FLD) channels.
    - Peak integration for FLD (5-7m and 11-13m).
    - Saves verification diagnostics and subtracted data.
    """
    sample_map = map_samples_to_blanks(files)
    if not sample_map:
        return

    # Define integration windows
    windows = [(6, 8), (14, 16)]
    results_log = []

    for s_path, b_path in sample_map.items():
        sample_name = os.path.basename(s_path)

        txt_out_dir = os.path.join(s_path, TXT_DIR)
        os.makedirs(txt_out_dir, exist_ok=True)

        targets = {
            "UIB1B_gdata.txt": ("Green Channel (UIB)", 0, False, None),
            "UIB1B_rdata.txt": ("Red Channel (UIB)", 0, False, None),
            "FLD1A_rawdata.txt": ("FLD Channel", 1, True, "FLD1A_subdata.txt")
        }

        fig, axes = plt.subplots(2, 1, figsize=figsize, sharex=True)
        data_found = False

        for fname, (label, row_idx, subtract, export_name) in targets.items():
            s_file = os.path.join(txt_out_dir, fname)
            b_file = os.path.join(b_path, TXT_DIR, fname)

            if not os.path.exists(s_file):
                s_file = os.path.join(s_path, fname)

            if not os.path.exists(s_file):
                continue

            try:
                s_raw = np.loadtxt(s_file)
                sx, sy = s_raw[:, 0], s_raw[:, 1]

                if subtract and os.path.exists(b_file):
                    b_raw = np.loadtxt(b_file)
                    bx, by = b_raw[:, 0], b_raw[:, 1]
                    interp_by = np.interp(sx, bx, by)
                    sy = np.maximum(sy - interp_by, 0)

                    plot_subtraction_verification(s_path, s_raw, b_raw,
                                                  np.column_stack((sx, sy)),
                                                  label, interp_blank=interp_by)

                    if export_name:
                        np.savetxt(os.path.join(txt_out_dir, export_name),
                                   np.column_stack((sx, sy)), header="time\tintensity")
                    label += " (Sub)"

                data_found = True
                ax = axes[row_idx]
                color = "green" if "Green" in label else ("red" if "Red" in label else "darkorchid")

                # Plotting: Raw scale for UIB (row 0), Normalization optional for FLD (row 1)
                ax.plot(sx, sy, label=label, color=color, linewidth=1.5)
                ax.fill_between(sx, sy, color=color, alpha=0.1)

                # --- Integration Logic for FLD (Row 1) ---
                # --- Integration Logic for FLD (Row 1) ---
                if row_idx == 1:
                    integration_text = "Integrals:\n"
                    # Prepare a dictionary for this specific sample
                    sample_results = {"Sample": sample_name}

                    for start, end in windows:
                        mask = (sx >= start) & (sx <= end)
                        try:
                            area = np.trapezoid(sy[mask], sx[mask])
                        except Exception:
                            area = np.trapz(sy[mask], sx[mask])

                        ax.fill_between(sx[mask], sy[mask], color=color, alpha=0.5)

                        # Store data for plot and for CSV
                        window_label = f"{start}-{end}m"
                        integration_text += f"{window_label}: {area:.2f}\n"
                        sample_results[window_label] = area

                    # Add text box to plot
                    ax.text(0.98, 0.5, integration_text, transform=ax.transAxes,
                            verticalalignment='center', horizontalalignment='right',
                            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8), fontsize=9)

                    # Append this sample's results to our master log
                    results_log.append(sample_results)

            except Exception as e:
                print(f" [!] Error processing {fname} for {sample_name}: {e}")


        if not data_found:
            plt.close(fig)
            continue

        # Formatting
        axes[0].set_title(f"Sample: {sample_name} | Dimension 1", fontweight='bold')
        axes[0].set_ylabel("Voltage (mV)")
        axes[0].legend(loc='upper right', fontsize='small')

        axes[1].set_title(f"Sample: {sample_name} | Dimension 2", fontweight='bold')
        axes[1].set_ylabel("Intensity (LU)")
        axes[1].set_xlabel("Time (min)")

        plt.tight_layout()

        if save:
            save_path = os.path.join(s_path, f"{sample_name}_report.png")
            save_path2 = save_path.replace(".png", ".pdf")
            plt.savefig(save_path2, dpi=300)
            plt.savefig(save_path, dpi=300)

            # Master Gallery copy
            # master_dir = os.path.join(topdir, "Report_Gallery")
            # os.makedirs(master_dir, exist_ok=True)
            # plt.savefig(os.path.join(master_dir, f"{sample_name}_report.png"), dpi=300)

        plt.close(fig)

    if results_log:
        import pandas as pd
        df = pd.DataFrame(results_log)
        # Reorder columns so 'Sample' is first
        cols = ['Sample'] + [c for c in df.columns if c != 'Sample']
        df = df[cols]

        csv_path = os.path.join(os.path.dirname(files[0]), "Dimension_Two_Summary.csv")
        df.to_csv(csv_path, index=False)
        print(f" Exported integration results to: {csv_path}")

def plot_multidetector_overlay(files, topdir, read_files=None, figsize=(16, 7), save=True):
    if not read_files:
        print("[plot_overlay] No files to read; aborting.")
        return


    effective_read_files = list(read_files)
    if any("UIB1B_rawdata" in f for f in effective_read_files):
        for dep in ["UIB1B_gdata.txt", "UIB1B_rdata.txt"]:
            if dep not in effective_read_files:
                effective_read_files.append(dep)

    data_dict = {rf: [] for rf in effective_read_files}
    names_dict = {rf: [] for rf in effective_read_files}


    for rf in effective_read_files:
        ext = os.path.splitext(rf)[1].lower()
        is_numpy = ext in (".npz", ".npy")
        for path in files:
            pref_dir = NPZ_DIR if is_numpy else TXT_DIR
            target = os.path.join(path, pref_dir, rf)
            if not os.path.exists(target): target = os.path.join(path, rf)
            if not os.path.exists(target): continue

            try:
                if is_numpy:
                    loaded = np.load(target, allow_pickle=True)
                    if hasattr(loaded, "files"):
                        for k in loaded.files: data_dict[rf].append(np.asarray(loaded[k]))
                    else:
                        data_dict[rf].append(np.asarray(loaded))
                else:
                    data_dict[rf].append(np.loadtxt(target))
                names_dict[rf].append(os.path.basename(path))
            except Exception as e:
                print(f" [!] Error loading {target}: {e}")


    processed_traces = {}
    for rf in effective_read_files:
        traces = []
        for arr, nm in zip(data_dict[rf], names_dict[rf]):
            x, y = parse_array_to_xy(arr)
            if x is not None:
                traces.append((x, y, nm))
        processed_traces[rf] = traces

    if not any(processed_traces.values()):
        print("[plot_overlay] No valid data found; aborting.")
        return


    n = len(read_files)
    fig, axes = plt.subplots(n, 1, figsize=(figsize[0] + 2, (figsize[1] * max(1, n)) + 2), squeeze=False)

    for i, rf in enumerate(read_files):
        ax = axes[i, 0]
        # Pull the main traces for this subplot
        current_traces = processed_traces.get(rf, [])

        # Plot Raw/Base Traces
        overlay_colors = get_filtered_colors(color_scheme.get("overlay", "viridis"), num_colors=len(current_traces))
        for idx, (x, y, label) in enumerate(current_traces):
            ax.plot(x, y, linewidth=0.8, alpha=0.3, label=f"{label} (Raw)", color=overlay_colors[idx], zorder=1)
            ax.fill_between(x, y, color=overlay_colors[idx], alpha=0.5, zorder=1)

        #Handle UIB Overlays (Accessing the "hidden" loaded dependencies)
        if "UIB1B_rawdata" in rf:
            channel_configs = [("UIB1B_gdata.txt", "gdata", "Green"), ("UIB1B_rdata.txt", "rdata", "Red")]
            for channel_key, scheme_key, suffix in channel_configs:
                sub_traces = processed_traces.get(channel_key, [])
                if sub_traces:
                    sub_colors = get_filtered_colors(color_scheme.get(scheme_key, "viridis"),
                                                     num_colors=len(sub_traces))
                    for j, (cx, cy, cnm) in enumerate(sub_traces):
                        ax.plot(cx, cy, linewidth=1.5, color=sub_colors[j], label=f"{cnm} ({suffix})", zorder=3,
                                alpha=0.8)

        # Formatting
        ax.set_title(f"Overlay: {rf}", fontweight='bold')
        ax.set_xlabel("Time (min)" if 'UIB' not in rf.upper() else "Time")
        ax.set_ylabel("Intensity" if 'UIB' not in rf.upper() else "Voltage (mV)")
        if current_traces:
            ax.legend(fontsize='x-small', loc='upper right', ncol=2, framealpha=0.8)

    plt.tight_layout()
    if save:
        outname = "overlay_" + "_".join([os.path.splitext(os.path.basename(rf))[0] for rf in read_files]) + ".png"
        outname2 = outname.replace(".png", ".pdf")
        fig.savefig(os.path.join(topdir, outname2), dpi=300)
        fig.savefig(os.path.join(topdir, outname), dpi=300)
    plt.show()


def plot_UIB_ridgelines_combined(analysis_results, topdir, sort_by='reps', v_step=1.0, scale_factor=0.8):
    print(f"[PLOTTING] Generating Combined UIB Ridgeline Plot (Sorted by: {sort_by})...")

    # 1. Sorting Logic
    if sort_by == 'reps':
        sort_key = lambda x: (analysis_results[x]['hours'], analysis_results[x].get('rep', 0))
    else:
        sort_key = lambda x: (analysis_results[x].get('rep', 0), analysis_results[x]['hours'])

    sorted_paths = sorted(analysis_results.keys(), key=sort_key)
    if not sorted_paths: return

    fig, ax = plt.subplots(figsize=(14, 10))

    # Colors
    g_colors = get_filtered_colors(color_scheme['gdata'], len(sorted_paths))
    r_colors = get_filtered_colors(color_scheme['rdata'], len(sorted_paths))

    # 2. Plotting Loop
    for i, p in enumerate(sorted_paths):
        baseline = i * v_step
        current_z = i
        res = analysis_results[p]

        # Green Channel
        if res['tg_sub'].size > 0:
            y_max = np.max(res['tg_sub'][:, 1])
            y_norm = (res['tg_sub'][:, 1] / y_max * scale_factor) if y_max > 0 else 0
            ax.plot(res['tg_sub'][:, 0], y_norm + baseline, color=g_colors[i], lw=1.5, zorder=current_z)
            ax.fill_between(res['tg_sub'][:, 0], y_norm + baseline, baseline, color=g_colors[i], alpha=0.3,
                            zorder=current_z)

        # Red Channel
        if res['tr_sub'].size > 0:
            y_max_r = np.max(res['tr_sub'][:, 1])
            y_norm_r = (res['tr_sub'][:, 1] / y_max_r * scale_factor) if y_max_r > 0 else 0
            ax.plot(res['tr_sub'][:, 0], y_norm_r + baseline, color=r_colors[i], lw=1.5, zorder=current_z)
            ax.fill_between(res['tr_sub'][:, 0], y_norm_r + baseline, baseline, color=r_colors[i], alpha=0.3,
                            zorder=current_z)

    # 3. Initial Formatting
    ax.set_title(f"UIB Ridgeline Analysis (Sorted by {sort_by.capitalize()})", fontweight='bold', fontsize=16)
    ax.set_xlabel("Time (min)", fontsize=12)
    ax.set_yticks([])  # Hide standard Y-ticks; colorbar will provide the scale
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)

    # 4. LOCK LAYOUT
    plt.subplots_adjust(left=0.08, right=0.82, top=0.92, bottom=0.1)
    fig.canvas.draw()

    # 5. Add Aligned Dual Colorbars
    pos = ax.get_position()
    norm = mcolors.Normalize(vmin=0, vmax=len(sorted_paths) - 1)

    # Green Bar
    cax_g = fig.add_axes([pos.x1 + 0.01, pos.y0, 0.012, pos.height])
    cb_g = ColorbarBase(cax_g, cmap=mcolors.ListedColormap(g_colors), norm=norm)
    cax_g.set_title("G", fontsize=10, fontweight='bold')
    cb_g.set_ticks([])  # USE THE COLORBAR OBJECT 'cb_g' HERE

    # Red Bar
    cax_r = fig.add_axes([pos.x1 + 0.06, pos.y0, 0.012, pos.height])
    cb_r = ColorbarBase(cax_r, cmap=mcolors.ListedColormap(r_colors), norm=norm)
    cax_r.set_title("R", fontsize=10, fontweight='bold')

    # Annotations: First and Last
    indices = [0, len(sorted_paths) - 1]
    labels = [
        f"t={analysis_results[sorted_paths[idx]]['hours']}h (r{analysis_results[sorted_paths[idx]].get('rep', 0)})" for
        idx in indices]

    cb_r.set_ticks(indices)  # USE THE COLORBAR OBJECT 'cb_r' HERE
    cb_r.set_ticklabels(labels, fontsize=9)

    # 6. Final Save
    plt.savefig(os.path.join(topdir, f"UIB_Ridgelines_{sort_by}.png"), dpi=300)
    plt.savefig(os.path.join(topdir, f"UIB_Ridgelines_{sort_by}.pdf"), dpi=300)
    plt.show()

def move_figures(files, topdir):
    print("[CLEANUP] Moving figures to subdirectories...")
    for d in files + [topdir]:
        f_dir = os.path.join(d, FIG_DIR)
        os.makedirs(f_dir, exist_ok=True)
        for ftype in ['*.png', '*.pdf']:
            for f_path in glob.glob(os.path.join(d, ftype)):
                if FIG_DIR not in f_path:
                    shutil.move(f_path, os.path.join(f_dir, os.path.basename(f_path)))

def main():
    rc_params_setup()
    topdir = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\Kinetics Replicates\20260218 POPC_POPG_LX_FLMSP rep3"
    files = glob.glob(os.path.join(topdir, "*.D"))
    detector_files = ["UIB1B_rawdata.npz", "FLD1A_rawdata.npz"]

    extract_all_files(files)
    convert_acq_to_utf8(files)
    parse_instrument_method(files)
    plot_instrument_timetables(files)

    analysis_results = analyze_UIB_samples(files, topdir, num_compressed=50)

    # Visualization & cleanup using results
    if analysis_results:
        plot_UIB_dashboard(analysis_results, topdir, sort_by='hours', kin_xaxis='hours', figsize=(20, 18)) #it is either 'reps' or 'hours'
        plot_UIB_ridgelines_combined(analysis_results, topdir, sort_by='hours', v_step=1.0, scale_factor=0.8)
    else:
        print("[MAIN] No analysis results generated; skipping dashboard/ridgeline plots.")

    plot_multidetector_overlay(files, topdir, read_files=detector_files, figsize=(16, 7), save=True)
    plot_individual_samples(files, figsize=(12, 10), save=True)
    move_figures(files, topdir)

    print("[Main] Pipeline Complete!")

if __name__ == "__main__":
    main()
