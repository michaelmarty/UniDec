import os
import numpy as np
import rainbow as rb
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import colormaps
from unidec.tools import integrate
from matplotlib import rcParams
import re
import colorsys
import glob

topdir = r"Z:\Group Share\Dash\Current Data\Agilent LC Data\LXEngine\20251030 LX\20251029 Control 2025-10-29 13-01-55"
files = glob.glob(os.path.join(topdir, "*.D"))


# Define the channel indices for multi-channel extraction only. Numbering starts from 0, so 3 and 4 correspond to the 4th and 5th channels in the data files.
channel_indices = [3, 4]

NPZ_DIR = "NPZ files"
TXT_DIR = "Text files"


#Read one of DAD1A.npz - DAD1B.npz - DAD1C - DAD1D - FLD1A.npz - UIB1B.npz
read_files=["FLD1A_data.npz", "UIB1B_data.npz"]

def get_filtered_colors(cmap_name='cool', num_colors=256, luminosity_factor=0.8, saturation_factor=2, hue_shift=0):
    '''

    cmap_name: Name of the matplotlib colormap to use (default is 'cool').
    num_colors: Number of colors to sample from the colormap.
    luminosity_factor: Multiplier for brightness (values < 1 darken, > 1 brighten).
    saturation_factor: Multiplier for color intensity (values < 1 desaturate, > 1 saturate).
    hue_shift: Value to shift the hue (color shade).
    '''

    viridis = colormaps.get_cmap(cmap_name)
    colors = []

    for i in range(num_colors):
        color = viridis(i / num_colors)  # Sample uniformly
        r, g, b, _ = color  # Unpack the color and ignore the alpha value
        h, l, s = colorsys.rgb_to_hls(r, g, b)


        new_hue = (h + hue_shift) % 1.0
        new_light = max(0, l * luminosity_factor)  # Ensure it doesn't go below 0
        new_sat = min(1, s * saturation_factor)  # Ensure it doesn't go below 0
        new_color = colorsys.hls_to_rgb(new_hue, new_light, new_sat)

        colors.append((new_color[0], new_color[1], new_color[2], color[3]))  # Retain the original alpha

    return colors

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

    # Ensure explicit tick label sizes (default to global font.size if unset)
    rcParams['xtick.labelsize'] = rcParams.get('xtick.labelsize', rcParams['font.size'])
    rcParams['ytick.labelsize'] = rcParams.get('ytick.labelsize', rcParams['font.size'])

    # Set legend font size to be 3 points smaller than tick label size and mark as bold.
    # Store bold weight in a small global helper so callers can apply it when creating legends.
    try:
        tick_size = float(rcParams['xtick.labelsize'])
    except Exception:
        tick_size = float(rcParams.get('font.size', 12))
    legend_size = max(1.0, tick_size - 3.0)
    rcParams['legend.fontsize'] = legend_size
    rcParams['legend.title_fontsize'] = legend_size
    globals()['LEGEND_FONT_WEIGHT'] = 'bold'

def get_just_UIB1B_files(files):
    """
    Save text and npz versions into `Text files` and `NPZ files` subdirs under each sample path.
    Sample name stored remains the basename of `path` (not the subdir).
    """
    for path in files:
        print(f"Processing: {path}")
        datadir = rb.read(path)
        datafile = datadir.get_file("UIB1B.ch")
        time = np.array(datafile.xlabels)
        voltages = np.array(datafile.data)

        # Ensure subdirs
        txt_outdir = os.path.join(path, TXT_DIR)
        npz_outdir = os.path.join(path, NPZ_DIR)
        os.makedirs(txt_outdir, exist_ok=True)
        os.makedirs(npz_outdir, exist_ok=True)

        # Save text
        txt_path = os.path.join(txt_outdir, "UIB1B_rawdata.txt")
        np.savetxt(txt_path, np.column_stack((time, voltages)), delimiter="\t")

        # Save npz with named array for consistency
        npz_path = os.path.join(npz_outdir, "UIB1B_rawdata.npz")
        np.savez_compressed(npz_path, data=np.column_stack((time, voltages)))

        # Save a quick preview png in the sample root (optional)
        plt.figure(figsize=(16, 9))
        plt.plot(time, voltages, marker='o', markersize=2, linestyle='-', linewidth=0.5, color='b', label='UIB1B Data')
        plt.savefig(os.path.join(path, "UIB1B_data_raw.png"))
        plt.close()

def extract_all_files(files):
    """
    Extract all .ch channel files and save text into `Text files` and npz into `NPZ files`.
    Skips `UIB1B.ch`.
    """
    skip_name = "UIB1B.ch"
    for path in files:
        print(f"Processing: {path}")
        datadir = rb.read(path)
        ch_files = [f for f in os.listdir(path) if f.lower().endswith('.ch') and f.lower() != skip_name]
        if not ch_files:
            continue

        txt_outdir = os.path.join(path, TXT_DIR)
        npz_outdir = os.path.join(path, NPZ_DIR)
        os.makedirs(txt_outdir, exist_ok=True)
        os.makedirs(npz_outdir, exist_ok=True)

        for ch_file in ch_files:
            print(f"Processing channel file: {ch_file}")
            try:
                datafile = datadir.get_file(ch_file)
            except Exception as e:
                print(f"  Failed to read channel file {ch_file}: {e}")
                continue

            time = np.array(datafile.xlabels)
            voltages = np.array(datafile.data)

            base_name = ch_file.replace('.ch', '_data')
            txt_path = os.path.join(txt_outdir, f"{base_name}.txt")
            np.savetxt(txt_path, np.column_stack((time, voltages)), delimiter="\t")

            npz_path = os.path.join(npz_outdir, f"{base_name}.npz")
            # Save with explicit key 'data' for predictable reading
            np.savez_compressed(npz_path, data=np.column_stack((time, voltages)))

def nonlinearize(datatop, num_compressed, gap=5):
    """
    Compress the data in a simple and nonlinear way.
    :param datatop: Data array (N x 2)
    :param num_compressed:
    :return:
    """
    if num_compressed == 0:
        return datatop
    else:
        num_compressed = int(num_compressed)
        gap = int(gap)
        return np.array([np.mean(datatop[index:index + num_compressed-gap], axis=0) for index in
                         range(0, len(datatop), num_compressed)])

def calc_start_index(data, numcompressed):
    adata = np.abs(data[:, 1])
    maxy = np.max(adata)
    ymaxindex = np.argmax(adata)
    threshold = 0.95 * maxy
    while ymaxindex > 0 and np.abs(data[ymaxindex, 1]) > threshold:
        ymaxindex -= 1
    return (ymaxindex+1) % numcompressed

def parse_rg_data(data, numcompressed, gap=5, negcolor="green", ramp=0):
    startindex = calc_start_index(data, numcompressed)
    print("Start Index:", startindex)

    parseddata = nonlinearize(data[startindex:], numcompressed, gap=gap)

    maxarg = np.argmax(parseddata[:, 1])
    minarg = np.argmin(parseddata[:, 1])

    print("Max Arg:", maxarg, "Min Arg:", minarg)

    nindex = minarg % 2
    pindex = maxarg % 2

    if nindex == pindex:
        raise ValueError("Negative and positive indices are the same, check the data.")

    ndata = parseddata[nindex::2]  # Negative data
    pdata = parseddata[pindex::2]  # Positive data

    #invert ndata
    ndata[:, 1] = ndata[:, 1] * -1  # Flip the negative channel to positive values

    if negcolor == "green":
        gdata = ndata
        rdata = pdata
    else:
        gdata = pdata
        rdata = ndata

    rdata[:, 1] = rdata[:, 1] - np.min(rdata[:, 1])  # Subtract the minimum to avoid negative values
    gdata[:, 1] = gdata[:, 1] - np.min(gdata[:, 1])  # Subtract the minimum to avoid negative values

    if ramp > 0:
        rdata[:, 1] = rdata[:, 1] * ramp

    print (rdata, gdata, parseddata)

    return rdata, gdata, parseddata

def process_UIB_data(files, numcompressed=50, gap=5, negcolor="magenta", ramp=0):
    for d in files:
        txt_outdir = os.path.join(d, TXT_DIR)
        os.makedirs(txt_outdir, exist_ok=True)

        txt_path = os.path.join(txt_outdir, "UIB1B_rawdata.txt")
        if not os.path.exists(txt_path):
            txt_path = os.path.join(d, "UIB1B_rawdata.txt")  # fallback to sample root
        if not os.path.exists(txt_path):
            print(f"No UIB1B_rawdata.txt found in {d}, skipping.")
            continue

        print("Processing file:", txt_path)
        data = np.loadtxt(txt_path)
        rdata, gdata, pdata = parse_rg_data(data, numcompressed, gap=gap, negcolor=negcolor, ramp=ramp)

        # Save plotted data into Text Files subdirectory
        g_out = os.path.join(txt_outdir, "UIB1B_gdata.txt")
        r_out = os.path.join(txt_outdir, "UIB1B_rdata.txt")

        header = "time\tintensity"
        if gdata is not None and gdata.size > 0:
            np.savetxt(g_out, gdata, delimiter="\t", header=header, comments='')
        else:
            with open(g_out, "w", encoding="utf-8") as f:
                f.write(header + "\n")

        if rdata is not None and rdata.size > 0:
            np.savetxt(r_out, rdata, delimiter="\t", header=header, comments='')
        else:
            with open(r_out, "w", encoding="utf-8") as f:
                f.write(header + "\n")

        print(f"Saved plotted data to `UIB1B_gdata.txt` and `UIB1B_rdata.txt` in `{txt_outdir}`")

        r1 = [2.5, 5]
        r2 = [8, 12]

        rps = np.round(calc_percents(gdata, r1, r2))
        gps = np.round(calc_percents(rdata, r1, r2))

        print("Red Channel Percent Unbound:", rps[0], "Bound:", rps[1])
        print("Green Channel Percent Unbound:", gps[0], "Bound:", gps[1])

        # Plot only when data present
        if rdata is not None and rdata.size > 0:
            plt.plot(rdata[:, 0], rdata[:, 1], label='Inverted Channel', color='magenta')
        if gdata is not None and gdata.size > 0:
            plt.plot(gdata[:, 0], gdata[:, 1], label='Uninverted Channel', color='lightseagreen')

        # Add vertical lines for the ranges
        plt.axvline(r1[0], color='grey', linestyle='--', label=f'Red Range {r1[0]}-{r1[1]}')
        plt.axvline(r1[1], color='grey', linestyle='--')
        plt.axvline(r2[0], color='blue', linestyle='--', label=f'Blue Range {r2[0]}-{r2[1]}')
        plt.axvline(r2[1], color='blue', linestyle='--')

        #plt.legend()
        # Save figure files
        plt.savefig(os.path.join(d, "UIB1B_rawdata.png"))
        plt.savefig( os.path.join(d, "UIB1B_rawdata.pdf"))
        #plt.show()
        plt.close()

def calc_percents(data, range1 = [2.2, 5], range2=[6.5,9]):
    """
    Calculate the percentage of data in the specified ranges.
    :param data: Data array (N x 2)
    :param range1: First range to calculate percentage
    :param range2: Second range to calculate percentage
    :return: Tuple of percentages in the two ranges
    """

    s1, _ = integrate(data, range1[0], range1[1])
    s2, _ = integrate(data, range2[0], range2[1])

    total = s1 + s2

    if total==0:
        percent1 = 0
        percent2 = 0
    else:
        percent1 = s1 / total * 100
        percent2 = s2 / total * 100

    return percent1, percent2

def read_data(files, read_files):
    """
    For each `rf` in `read_files` decide which subdir to search based on extension:
      - text-like -> `Text files` then sample root
      - npy/npz -> `NPZ files` then sample root
      - no/unknown ext -> try `Text files`, `NPZ files`, then sample root

    Always record sample name as basename(path) (not the subdir).
    """
    print(f"[read_data] start: read_files={read_files}, total_paths={len(files)}")
    if isinstance(read_files, str):
        read_files = [read_files]
    if read_files is None:
        read_files = []
    # filter out falsy
    read_files = [rf for rf in read_files if rf]

    data = {rf: [] for rf in read_files}
    file_names = {rf: [] for rf in read_files}

    text_exts = (".txt", ".csv", ".dat", ".asc", ".tsv")
    np_exts = (".npz", ".npy")

    for path in files:
        print(f"[read_data] checking path: {path}")
        for rf in read_files:
            rf_ext = os.path.splitext(rf)[1].lower()
            # choose candidate directories in preferred order
            if rf_ext in np_exts:
                candidate_dirs = [NPZ_DIR, ""]   # prefer NPZ subdir then root
            elif rf_ext in text_exts:
                candidate_dirs = [TXT_DIR, ""]   # prefer Text subdir then root
            else:
                # unknown extension: try text, npz subdirs first, then root
                candidate_dirs = [TXT_DIR, NPZ_DIR, ""]

            found = False
            for cand_dir in candidate_dirs:
                if cand_dir:
                    target_file = os.path.join(path, cand_dir, rf)
                else:
                    target_file = os.path.join(path, rf)

                print(f"[read_data] looking for {rf} -> {target_file}")
                if not os.path.exists(target_file):
                    continue

                try:
                    ext = os.path.splitext(target_file)[1].lower()
                    if ext in (".npy", ".npz"):
                        print(f"[read_data] loading binary numpy file: {target_file}")
                        arr = np.load(target_file, allow_pickle=True)
                        if hasattr(arr, 'files'):  # npz: may contain named arrays
                            for k in arr.files:
                                load_data = np.asarray(arr[k])
                                data[rf].append(load_data)
                                file_names[rf].append(os.path.basename(path))
                                print(f"    loaded array '{k}' shape={getattr(load_data, 'shape', None)}")
                        else:
                            load_data = np.asarray(arr)
                            data[rf].append(load_data)
                            file_names[rf].append(os.path.basename(path))
                            print(f"    loaded simple array shape={getattr(load_data, 'shape', None)}")
                    else:
                        # treat as text
                        print(f"[read_data] loading text file: {target_file}")
                        load_data = np.loadtxt(target_file)
                        load_data = np.asarray(load_data)
                        data[rf].append(load_data)
                        file_names[rf].append(os.path.basename(path))
                        print(f"    loaded text array shape={getattr(load_data, 'shape', None)}")
                    found = True
                    break  # stop searching candidate dirs for this rf/path
                except Exception as e:
                    import traceback
                    tb = traceback.format_exc()
                    print(f"[read_data] Failed to read {target_file}: {e}\n{tb}")
                    # continue to next candidate dir

            if not found:
                print(f"[read_data] no file found for {rf} in {path}")

    print(f"[read_data] finished. counts: " + ", ".join(f"{rf}={len(data[rf])}" for rf in read_files))
    return data, file_names

def parse_array_to_xy(arr):
    arr = np.asarray(arr)
    print(f"[parse_array_to_xy] input type={type(arr)} dtype={arr.dtype if hasattr(arr, 'dtype') else 'N/A'} ndim={arr.ndim} shape={getattr(arr, 'shape', None)}")
    if arr.ndim == 2 and arr.shape[1] == 2:
        x, y = arr[:, 0], arr[:, 1]
        print(f"[parse_array_to_xy] returning x.shape={x.shape} y.shape={y.shape}")
        return x, y
    print("[parse_array_to_xy] input is not a 2-column 2D array; returning (None, None)")
    return None, None

def plot_overlay(files, topdir, read_files=None, figsize=(16, 7), save=True):
    """
    Load arrays named in `read_files` (or defaults) using `read_data`, parse them
    into (x, y) and plot overlays with a dynamic number of subplots equal to
    the length of `read_files`.
    """
    if read_files is None:
        # get module-level default variable `read_files` if it exists
        module_default = globals().get('read_files')
        if module_default:
            read_files = module_default
            print("[plot_overlay] using module-level default `read_files`")
        else:
            print("[plot_overlay] no `read_files` provided and no module default; aborting")
            return
    elif isinstance(read_files, str):
        read_files = [read_files]

    # Filter out any falsy entries (None, empty strings)
    read_files = [rf for rf in read_files if rf]
    if not read_files:
        print("[plot_overlay] read_files is empty after normalization; aborting")
        return

    print(f"[plot_overlay] loading read_files={read_files} files_count={len(files)} topdir={topdir}")
    data_dict, names_dict = read_data(files, read_files)
    print("[plot_overlay] loaded data counts: " + ", ".join(
        f"{rf}={len(data_dict.get(rf, []))}" for rf in read_files))

    # Build traces per requested file type
    traces_per_file = []
    for rf in read_files:
        arrs = data_dict.get(rf, [])
        nms = names_dict.get(rf, [])
        traces = []
        for arr, nm in zip(arrs, nms):
            x, y = parse_array_to_xy(arr)
            if x is not None and y is not None:
                traces.append((x, y, nm))
        traces_per_file.append(traces)

    if all(len(t) == 0 for t in traces_per_file):
        print(f"No data found for {read_files} in the provided `files`.")
        return

    n = len(read_files)
    base_w, base_h = figsize
    fig_w = base_w
    fig_h = base_h * max(1, n)
    fig, axes = plt.subplots(n, 1, figsize=(fig_w, fig_h), squeeze=False)

    for i, rf in enumerate(read_files):
        ax = axes[i, 0]
        traces = traces_per_file[i]
        for x, y, label in traces:
            ax.plot(x, y, linewidth=0.9, alpha=0.8, label=label)

        #ax.set_title(f"Overlay: `{rf}`")
        rf_up = rf.upper()
        if 'UIB' in rf_up:
            ax.set_xlabel("Time")
            ax.set_ylabel("Voltage (mV)")
        else:
            ax.set_xlabel("Retention Time (min)")
            ax.set_ylabel("Intensity")

        if traces:
            ax.legend(fontsize='small', loc='upper right')

    plt.tight_layout()
    if save:
        outname = "overlay_" + "_".join([os.path.splitext(os.path.basename(rf))[0] for rf in read_files]) + ".png"
        outpath = os.path.join(topdir, outname)
        fig.savefig(outpath, dpi=300)
        print(f"Saved overlay figure to ` {outpath} `")
    plt.show()
    plt.close(fig)

def convert_acq_to_utf8(files):
    """
    Converts each acq.txt/acqmethod.txt from UTF-16-LE to UTF-8
    and saves as 'acq_utf8.txt' in the same folder.
    Returns list of new UTF-8 file paths.
    """
    utf8_files = []
    for d_path in files:
        for root, _, fnames in os.walk(d_path):
            for fname in fnames:
                if fname.lower() in ("acq.txt", "acqmethod.txt"):
                    src_path = os.path.join(root, fname)
                    dst_path = os.path.join(root, "acq_utf8.txt")
                    with open(src_path, "rb") as f:
                        raw = f.read()
                    # detect UTF-16-LE
                    if b"\x00" in raw:
                        text = raw.decode("utf-16-le", errors="ignore")
                    else:
                        text = raw.decode("utf-8", errors="ignore")
                    # remove any stray nulls or NBSPs
                    text = re.sub(r"[\x00\xa0]", "", text)
                    with open(dst_path, "w", encoding="utf-8") as f:
                        f.write(text)
                    utf8_files.append(dst_path)
                    print(f"Converted {fname} → {os.path.basename(dst_path)}")
    return utf8_files

def parse_instrument_method(files, debug=True):
    """
    Parses 'acq_utf8.txt' from each .D directory (no fallback).
    Extracts Timetable sections from each module and saves them into
    'Method Timetables/<Module Name>.txt'.
    """

    def read_text(path):
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            text = f.read()
        # remove any stray nulls or NBSPs
        text = re.sub(r"[\x00\xa0]", "", text)
        return text

    name_line_re = re.compile(r"^[^\S\r\n]*[A-Za-z].*\(G\d{4}[A-Z]\)\s*$")

    for d_path in files:
        print(f"\n Parsing instrument method in: {d_path}")

        #  Require acq_utf8.txt
        acq_file = None
        for root, _, fnames in os.walk(d_path):
            for fname in fnames:
                if fname.lower() == "acq_utf8.txt":
                    acq_file = os.path.join(root, fname)
                    break
            if acq_file:
                break

        if not acq_file:
            print(f" No acq_utf8.txt found in {d_path}, skipping.")
            continue

        print(f" Reading: {os.path.basename(acq_file)}")
        text = read_text(acq_file)
        lines = text.splitlines()

        # Skip metadata
        if len(lines) > 7:
            lines = lines[7:]

        #  Find module name lines
        name_idx = [(i, " ".join(line.strip().split()))
                    for i, line in enumerate(lines)
                    if name_line_re.match(line)]

        print(f" Detected {len(name_idx)} module(s)")
        if not name_idx:
            continue

        out_dir = os.path.join(d_path, "Method Timetables")
        os.makedirs(out_dir, exist_ok=True)

        saved = 0
        for k, (start_i, module_name) in enumerate(name_idx):
            end_i = name_idx[k + 1][0] if k + 1 < len(name_idx) else len(lines)
            body = "\n".join(lines[start_i + 1:end_i])

            #  Find Timetable section
            m = re.search(r"(?i)\bTimetable\b", body)
            if not m:
                if debug:
                    print(f"    No timetable in {module_name}")
                continue

            after = body[m.end():]

            # Find end boundary
            cut_points = []
            inner_mod = name_line_re.search(after)
            if inner_mod:
                cut_points.append(inner_mod.start())
            inst_curves = re.search(r"(?im)^\s*Instrument\s+Curves\b", after)
            if inst_curves:
                cut_points.append(inst_curves.start())
            next_sep = re.search(r"\n=+\s*\n", after)
            if next_sep:
                cut_points.append(next_sep.start())

            timetable = after[:min(cut_points)].strip() if cut_points else after.strip()

            lines_tt = [ln for ln in timetable.splitlines() if ln.strip()]
            if len(lines_tt) < 3:
                if debug:
                    print(f"    Empty timetable for {module_name}")
                continue

            #  Save cleanly
            safe_name = re.sub(r"[^\w\s\(\)\-\.]", "", module_name).strip()
            out_path = os.path.join(out_dir, f"{safe_name}.txt")
            with open(out_path, "w", encoding="utf-8") as f:
                f.write(timetable)
            saved += 1
            print(f" Saved {os.path.basename(out_path)} ({len(lines_tt)} lines)")

        print(f" Done parsing {os.path.basename(acq_file)} — {saved} timetable(s) written\n")

def plot_instrument_timetables(files):
    """
    Creates four fixed subplots:
      1. Alpha Pump (color-coded by dominant solvent)
      2. Beta Pump  (color-coded by dominant solvent)
      3. Gamma Pump (Flow + %B)
      4. Column Compartments (G7116A + G1316A)
    """

    def forward_fill(values):
        """Replace None or '---' with the last known numeric value."""
        last = None
        filled = []
        for v in values:
            if v is None or (isinstance(v, str) and v.strip() == "---"):
                filled.append(last)
            else:
                try:
                    val = float(v)
                    filled.append(val)
                    last = val
                except ValueError:
                    filled.append(last)
        return filled

    cmap_map = {
        "alpha": "winter",
        "beta": "cool",
        "gamma": "Wistia"
    }

    for d_path in files:
        timetable_dir = os.path.join(d_path, "Method Timetables")
        if not os.path.isdir(timetable_dir):
            print(f"⚠️ No Method Timetables found in {d_path}")
            continue

        # Explicit file selection
        alpha_file = next((f for f in os.listdir(timetable_dir) if "alpha" in f.lower()), None)
        beta_file  = next((f for f in os.listdir(timetable_dir) if "beta" in f.lower()), None)
        gamma_file = next((f for f in os.listdir(timetable_dir) if "gamma" in f.lower()), None)
        col1_file  = next((f for f in os.listdir(timetable_dir) if "g7116a" in f.lower()), None)
        col2_file  = next((f for f in os.listdir(timetable_dir) if "g1316a" in f.lower()), None)

        fig, axs = plt.subplots(4, 1, figsize=(10, 12), sharex=False)
        fig.subplots_adjust(hspace=0.4)
        fig.suptitle("Instrument Method Overview", fontsize=14)

        # ---------------- Helper for Quaternary Pumps ----------------
        def parse_quat_timetable(fpath):
            with open(fpath, "r", encoding="utf-8") as f:
                lines = [ln.strip() for ln in f if ln.strip()]
            data_start = next((i for i, l in enumerate(lines) if re.match(r"^-{3,}$", l)), None)
            data_lines = lines[data_start + 1:] if data_start else lines
            t, A, B, C, D, Flow = [], [], [], [], [], []
            for line in data_lines:
                parts = re.split(r"\s+", line.strip())
                if parts[0].lower().startswith("start"):
                    time = 0.0
                else:
                    try:
                        time = float(parts[0])
                    except ValueError:
                        continue
                vals = [None if p == "---" else p for p in parts[1:6]]
                if len(vals) < 5:
                    continue
                t.append(time)
                A.append(vals[0])
                B.append(vals[1])
                C.append(vals[2])
                D.append(vals[3])
                Flow.append(vals[4])
            A, B, C, D, Flow = map(forward_fill, [A, B, C, D, Flow])
            return np.array(t), np.array(A, float), np.array(B, float), np.array(C, float), np.array(D, float), np.array(Flow, float)

        # ---------------- Alpha Pump ----------------
        if alpha_file:
            t, A, B, C, D, Flow = parse_quat_timetable(os.path.join(timetable_dir, alpha_file))
            ax = axs[0]
            cmap = plt.get_cmap(cmap_map["alpha"])
            colors = [cmap(x) for x in np.linspace(0.1, 0.9, 4)]
            for i in range(len(t) - 1):
                vals = [A[i], B[i], C[i], D[i]]
                dominant = np.argmax(vals)
                ax.plot(t[i:i+2], Flow[i:i+2], color=colors[dominant], linewidth=3)
            ax.set_xlim(0, 15)
            ax.set_ylim(0, 1.0)
            ax.set_title("Alpha Pump", fontsize=10)
            ax.set_ylabel("Flow (mL/min)")
            ax.grid(False)

        # ---------------- Beta Pump ----------------
        if beta_file:
            t, A, B, C, D, Flow = parse_quat_timetable(os.path.join(timetable_dir, beta_file))
            ax = axs[1]
            cmap = plt.get_cmap(cmap_map["beta"])
            colors = [cmap(x) for x in np.linspace(0.1, 0.9, 4)]
            for i in range(len(t) - 1):
                vals = [A[i], B[i], C[i], D[i]]
                dominant = np.argmax(vals)
                ax.plot(t[i:i+2], Flow[i:i+2], color=colors[dominant], linewidth=3)
            ax.set_xlim(0, 15)
            ax.set_ylim(0, 1.0)
            ax.set_title("Beta Pump", fontsize=10)
            ax.set_ylabel("Flow (mL/min)")
            ax.grid(False)

        # ---------------- Gamma Pump ----------------
        if gamma_file:
            with open(os.path.join(timetable_dir, gamma_file), "r", encoding="utf-8") as f:
                lines = [ln.strip() for ln in f if ln.strip()]
            data_start = next((i for i, l in enumerate(lines) if re.match(r"^-{3,}$", l)), None)
            data_lines = lines[data_start + 1:] if data_start else lines
            t, A, B, Flow = [], [], [], []
            for line in data_lines:
                parts = re.split(r"\s+", line.strip())
                if parts[0].lower().startswith("start"):
                    time = 0.0
                else:
                    try:
                        time = float(parts[0])
                    except ValueError:
                        continue
                vals = [None if p == "---" else p for p in parts[1:4]]
                if len(vals) < 3:
                    continue
                t.append(time)
                A.append(vals[0])
                B.append(vals[1])
                Flow.append(vals[2])
            B = forward_fill(B)
            Flow = forward_fill(Flow)
            ax = axs[2]
            ax.plot(t, Flow, color="black", linewidth=2, label="Flow")
            axr = ax.twinx()
            axr.plot(t, B, color="tab:orange", linestyle="--", linewidth=2, label="%B")
            ax.set_xlim(0, 15)
            ax.set_ylim(0, 1.5)
            ax.set_ylabel("Flow (mL/min)")
            axr.set_ylabel("%B")
            ax.set_title("Gamma Pump", fontsize=10)
            ax.grid(False)

        # ---------------- Column Compartments ----------------
        ax = axs[3]
        ax.set_title("Column Compartments", fontsize=10)
        colors = {"G7116A": "tab:blue", "G1316A": "tab:orange"}
        for fname, color in [(col1_file, colors["G7116A"]), (col2_file, colors["G1316A"])]:
            if not fname:
                continue
            with open(os.path.join(timetable_dir, fname), "r", encoding="utf-8") as f:
                lines = [ln.strip() for ln in f if ln.strip()]
            times = []
            # find where the timetable data begins (after dashed line)
            data_start = next((i for i, l in enumerate(lines) if re.match(r"^-{3,}$", l)), None)
            data_lines = lines[data_start + 1:] if data_start is not None else lines

            for ln in data_lines:
                ln = ln.strip()
                if not ln:
                    continue  # skip empty lines
                ln = re.sub(r"(\d)([A-Za-z])", r"\1 \2", ln)
                parts = re.split(r"\s+", ln)
                # must start with a float (time)
                try:
                    t = float(parts[0])
                    times.append(t)
                except ValueError:
                    continue

            if times:
                ax.scatter(times, [1]*len(times), color=color, s=60, edgecolor='k', label=fname.replace('.txt',''))
        ax.set_xlim(0, 15)
        ax.set_yticks([])
        ax.set_xlabel("Time (min)")
        ax.set_ylabel("Column Valve Events")
        ax.legend(fontsize=8)
        ax.grid(False)

        plt.tight_layout(rect=[0, 0, 1, 0.96])
        out_path = os.path.join(d_path, "instrument_method.png")
        plt.savefig(out_path, dpi=200)
        plt.close(fig)
        print(f" Saved 4-panel instrument method plot: {out_path}")

def plot_UIB_overlay(files, topdir, figsize=(14, 8), save=True):
    """
    Reads all 'UIB1B_gdata.txt' and 'UIB1B_rdata.txt' files from each
    'Text files' subdirectory in the provided .D folders, then plots
    a figure with 2 subplots:
      - Top: overlay of all gdata traces (summer colormap)
      - Bottom: overlay of all rdata traces (cool colormap)
    Sorted alphabetically by the substring after the last '-' and before '.D',
    which is also used for the legend labels.
    """
    # --- Extract substring after the last '-' and before '.D' ---
    def extract_key_and_sort(path):
        base = os.path.basename(path)
        m = re.search(r'-([^-\\/.]+)\.D$', base, re.IGNORECASE)
        key = m.group(1).lower() if m else base.lower()
        return key

    # Build sorted list and keep mapping
    keys = {d: extract_key_and_sort(d) for d in files}
    files_sorted = sorted(files, key=lambda p: keys[p])

    g_all, r_all, labels = [], [], []

    for d in files_sorted:
        txtdir = os.path.join(d, "Text files")
        gpath = os.path.join(txtdir, "UIB1B_gdata.txt")
        rpath = os.path.join(txtdir, "UIB1B_rdata.txt")
        label = keys[d]  # legend label = extracted key

        def safe_loadtxt(path):
            if not os.path.exists(path):
                return None
            try:
                # Skip header if it starts with 'time'
                with open(path, 'r', encoding='utf-8') as f:
                    first = f.readline().strip().lower()
                skip = 1 if first.startswith("time") else 0
                return np.loadtxt(path, delimiter="\t", skiprows=skip)
            except Exception as e:
                print(f"Failed to read {path}: {e}")
                return None

        gdata = safe_loadtxt(gpath)
        rdata = safe_loadtxt(rpath)

        if gdata is not None and gdata.size > 0:
            g_all.append(gdata)
            labels.append(label)
        if rdata is not None and rdata.size > 0:
            r_all.append(rdata)

    if not g_all and not r_all:
        print("No gdata or rdata files found — nothing to plot.")
        return

    fig, (axg, axr) = plt.subplots(2, 1, figsize=figsize, sharex=False)
    fig.suptitle("UIB1B", fontsize=16)

    # Separate colormaps
    g_colors = get_filtered_colors('summer_r', num_colors=max(len(g_all), 1))
    r_colors = get_filtered_colors('spring_r', num_colors=max(len(r_all), 1))

    # --- gdata overlay (summer) ---
    for i, arr in enumerate(g_all):
        if arr.ndim == 2 and arr.shape[1] == 2:
            axg.plot(arr[:, 0], arr[:, 1],
                     label=labels[i],
                     color=g_colors[i % len(g_colors)],
                     linewidth=1)
    axg.set_title("Green Channel (gdata)")
    axg.set_xlabel("Time")
    axg.set_ylabel("Intensity")
    axg.legend(fontsize='small', loc='upper right')

    # --- rdata overlay (cool) ---
    for i, arr in enumerate(r_all):
        if arr.ndim == 2 and arr.shape[1] == 2:
            axr.plot(arr[:, 0], arr[:, 1],
                     label=labels[i] if i < len(labels) else f"Sample {i+1}",
                     color=r_colors[i % len(r_colors)],
                     linewidth=1)
    axr.set_title("Red Channel (rdata)")
    axr.set_xlabel("Time")
    axr.set_ylabel("Intensity")
    axr.legend(fontsize='small', loc='upper right')

    plt.tight_layout(rect=[0, 0, 1, 0.97])
    if save:
        outpath = os.path.join(topdir, "UIB1B_g_r_overlay_sorted.png")
        fig.savefig(outpath, dpi=300)
        print(f" Saved sorted overlay figure: {outpath}")

    plt.show()
    plt.close(fig)

def plot_g_r_ridgelines(files, topdir, base_offset=100, autoscale=True, padding_factor=0.15):
    """
    Creates two separate ridgeline plots (not subplots):
      - One for all gdata (summer_r colormap)
      - One for all rdata (spring_r colormap)

    Sorting is based on the substring after the last '-' and before '.D'.
    Offset can be fixed (base_offset) or automatically adjusted to just above
    each dataset's maximum y-intensity if autoscale=True.

    padding_factor defines how much taller each offset is than the previous max (default 5%).
    """
    # --- Extract label between last '-' and '.D' ---
    def extract_label(path):
        base = os.path.basename(path)
        m = re.search(r'-([^-\\/.]+)\.D$', base, re.IGNORECASE)
        return m.group(1).lower() if m else base.lower()

    # --- Safely load text data, skipping header if needed ---
    def safe_loadtxt(path):
        if not os.path.exists(path):
            return None
        try:
            with open(path, 'r', encoding='utf-8') as f:
                first = f.readline().strip().lower()
            skip = 1 if first.startswith("time") else 0
            return np.loadtxt(path, delimiter="\t", skiprows=skip)
        except Exception as e:
            print(f"Failed to read {path}: {e}")
            return None

    # --- Collect and sort files ---
    labels = {d: extract_label(d) for d in files}
    files_sorted = sorted(files, key=lambda p: labels[p])

    g_data, r_data, g_labels, r_labels = [], [], [], []

    for d in files_sorted:
        txtdir = os.path.join(d, "Text files")
        gpath = os.path.join(txtdir, "UIB1B_gdata.txt")
        rpath = os.path.join(txtdir, "UIB1B_rdata.txt")
        label = labels[d]

        g_arr = safe_loadtxt(gpath)
        r_arr = safe_loadtxt(rpath)

        if g_arr is not None and g_arr.size > 0:
            g_data.append(g_arr)
            g_labels.append(label)
        if r_arr is not None and r_arr.size > 0:
            r_data.append(r_arr)
            r_labels.append(label)

    # --- Internal helper for one ridgeline group ---
    def plot_ridgeline(data_list, label_list, cmap_name, title, save_prefix):
        if not data_list:
            print(f"No data found for {title}")
            return

        colors = get_filtered_colors(cmap_name, num_colors=len(data_list))
        plt.figure(figsize=(12, 12))

        y_cursor = 0
        for i, (arr, label) in enumerate(zip(data_list, label_list)):
            if arr is None or arr.size == 0:
                continue
            x, y = arr[:, 0], arr[:, 1]
            color = colors[i % len(colors)]

            # Compute dynamic offset
            if autoscale:
                y_shifted = y + y_cursor
                ymax = np.max(y)
                y_cursor += ymax * (1 + padding_factor)
            else:
                y_shifted = y + y_cursor
                y_cursor += base_offset

            #plt.fill_between(x, y_shifted, y_cursor - (np.max(y) if autoscale else base_offset),alpha=0.4, color=color)
            plt.plot(x, y_shifted, color=color, linewidth=0.7)
            plt.text(0, y_shifted[0] + 0.05 * (np.max(y)), label,
                     fontsize=8, fontweight='bold', ha='left', color=color)

        plt.xlabel("Time (min)", fontweight='bold')
        plt.ylabel("")
        plt.title(title + (" (autoscaled)" if autoscale else ""), fontweight='bold')
        plt.yticks([])
        ax = plt.gca()
        for spine in ['top', 'right', 'left']:
            ax.spines[spine].set_visible(False)
        plt.subplots_adjust(hspace=-0.3)

        # Save output
        out_png = os.path.join(topdir, f"{save_prefix}_ridgeline.png")
        out_pdf = os.path.join(topdir, f"{save_prefix}_ridgeline.pdf")
        plt.savefig(out_png, dpi=300, bbox_inches='tight')
        plt.savefig(out_pdf, dpi=300, bbox_inches='tight')
        print(f" Saved ridgeline plots:\n  {out_png}\n  {out_pdf}")

        plt.show()
        plt.close()

    # --- Generate both ridgeline plots ---
    if g_data:
        plot_ridgeline(
            g_data, g_labels,
            cmap_name='winter_r',
            title="UIB1B gdata",
            save_prefix="UIB1B_gdata"
        )

    if r_data:
        plot_ridgeline(
            r_data, r_labels,
            cmap_name='spring_r',
            title="UIB1B rdata",
            save_prefix="UIB1B_rdata"
        )



def main():
    rc_params_setup()
    extract_all_files(files)
    get_just_UIB1B_files(files)
    read_data(files, read_files)
    convert_acq_to_utf8(files)
    parse_instrument_method(files, debug=True)
    plot_instrument_timetables(files)
    process_UIB_data(files, numcompressed=50, gap=5, negcolor="magenta", ramp=0)
    plot_overlay(files, topdir, read_files=read_files, figsize=(16, 7), save=True)
    plot_UIB_overlay(files, topdir, figsize=(14, 8), save=True)
    plot_g_r_ridgelines(files, topdir, base_offset=100, autoscale=True)


if __name__ == "__main__":
    main()



