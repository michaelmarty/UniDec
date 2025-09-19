import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rcParams
from matplotlib.patches import Rectangle
import unidec.tools as ud

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 14
rcParams['font.sans-serif'] = "Arial"


def get_uddir(path):
    dirnew = os.path.splitext(path)[0] + "_unidecfiles"
    return dirnew


def get_mass_data(path):
    header = os.path.splitext(path)[0]
    uddir = get_uddir(path)
    mfile = os.path.join(uddir, header + "_mass.txt")
    mdata = np.loadtxt(mfile)
    return mdata


def make_plotter_file(dir, pfile="plotter.csv"):
    files = os.listdir()
    udfiles = []
    for f in files:
        if os.path.splitext(f)[1] == ".txt":
            uddir = get_uddir(f)
            if os.path.isdir(uddir):
                udfiles.append(f)
    udfiles = np.array(udfiles)
    if not os.path.isfile(pfile):
        np.savetxt(pfile, udfiles, delimiter=",", fmt="%s")


def read_plotter_file(dir, pfile="plotter.csv"):
    fdata = np.loadtxt(os.path.join(dir, pfile), delimiter=",", dtype=str)
    print(fdata)
    return fdata

def write_data_to_plotter(dir, extracts, pfile="plotter.csv"):
    fdata = np.loadtxt(os.path.join(dir, pfile), delimiter=",", dtype=str)
    newfdata = []
    for i, f in enumerate(fdata):
        f = np.array(f)
        x = np.array(extracts[i], dtype=np.str)
        newf = np.concatenate([f, x])
        print(newf)
        newfdata.append(newf)
    np.savetxt(pfile, newfdata, delimiter=",", fmt="%s")


def add_rect(ax, xstart, ystart, xwidth, ywidth, alpha=0.5, facecolor="k", edgecolor='k'):
    ax.add_patch(Rectangle((xstart, ystart), xwidth, ywidth, alpha=alpha, facecolor=facecolor, edgecolor=edgecolor,
                  fill=True))

def plot_fdata(fdata, mmass=11810):
    rows = fdata[:, 1].astype(int)
    cols = fdata[:, 2].astype(int)
    nrows = len(np.unique(rows))
    ncols = len(np.unique(cols))

    fig, ax = plt.subplots(nrows, ncols, sharex='all', sharey='all', figsize=[12, 5]) #Adjust the Figure size here (WxH)
    ydataall=[]
    for i in range(0, len(fdata)):
        f = fdata[i, 0]
        r = rows[i]
        c = cols[i]

        mdata = get_mass_data(f)
        xdata = mdata[:, 0] / 1000.
        ydata = mdata[:, 1] / np.amax(mdata[:, 1]) * 100.

        ax[r, c].plot(xdata, ydata, color="k")
        yvals = []
        if mmass is not None:
            mass = mmass/1000.
            masslist = mass * np.arange(1, 7)
            colors = ["r", "g", "b", "purple", "y", "m"]
            for j, mass in enumerate(masslist):
                window=2
                xstart = mass - window
                xwidth = window * 0.2
                ystart = 0

                yval = ud.data_extract(np.transpose([xdata, ydata]), mass, extract_method=2, window=window)
                ywidth = yval #+10
                yvals.append(yval)
                add_rect(ax[r,c], xstart, ystart, xwidth, ywidth, facecolor=colors[j], alpha=0.4)
        yvals = np.array(yvals)
        yvals /= np.amax(yvals)/100.
        ax[r, c].spines['top'].set_visible(False)
        ax[r, c].spines['right'].set_visible(False)

        ax[r, c].set_xlim(0, 80)  # Set X limits here
        ax[r, c].set_ylim(0, 100)

        ax[r, c].get_yaxis().set_ticks([0, 50, 100])
        ax[r, c].get_yaxis().set_ticklabels(["0", '%', "100"])

        if r == np.amax(rows):
            ax[r, c].set_xlabel("Mass (kDa)")

        ydataall.append(yvals)
    plt.tight_layout()
    plt.savefig("plotter.pdf", transparent="True")
    plt.show()

    return ydataall



monomermass=11810

dir = "Z:\\Julia\\M2 protein\\wildtype\\NEW PURIFICATION\\210104\\detergentscreendata\\3"


if not os.path.isdir(dir):
    print("Not a directory")
    print(dir)
    exit()

os.chdir(dir)

make_plotter_file(dir)
fdata = read_plotter_file(dir)
if len(fdata.shape) == 2:
    yvals = plot_fdata(fdata)
    print(fdata)
    write_data_to_plotter(dir, yvals)
else:
    print("Need to update plotter.csv file to specify rows and columns")
    opencommand = "start \"\" "
    os.system(opencommand + "\"" + dir + "\"")

