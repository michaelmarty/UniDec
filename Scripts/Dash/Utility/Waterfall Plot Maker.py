import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib import colormaps
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from matplotlib import rcParams

rcParams['ps.useafm'] = True
rcParams['ps.fonttype'] = 42
rcParams['pdf.fonttype'] = 42
rcParams['lines.linewidth'] = 1
rcParams['errorbar.capsize'] = 3
rcParams['patch.force_edgecolor'] = True
rcParams['patch.facecolor'] = 'b'
rcParams['lines.markersize'] = 7
rcParams['font.size'] = 12
rcParams['font.weight'] = 500
rcParams['font.sans-serif'] = "Arial"
rcParams['axes.titleweight'] = 500
rcParams['axes.titlesize'] = 8
rcParams['axes.linewidth'] = 3
rcParams['xtick.major.width'] = 1
rcParams['ytick.major.width'] = 1

def read_data(file_path):
    """Read x and y data from a text file."""
    data = np.loadtxt(file_path)
    x = data[:, 0]
    y = data[:, 1]
    return x, y

def format_mass_axis(x, pos):
    """Format the mass axis to use 'K' for multiples of 1000."""
    if x >= 1000:
        return f'{int(x / 1000)}K'
    return f'{int(x)}'

def create_waterfall_plot(data_list, output_file):
    """Create a waterfall plot from a list of (x, y) data tuples."""
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    viridis = colormaps.get_cmap('viridis').resampled(len(data_list))

    for i, (x, y, label) in enumerate(data_list):
        ax.plot(x, y, zs=i, zdir='y', label=label, color=viridis(i))

        # Create vertices for the filled area
        verts = [(x[j], i, y[j]) for j in range(len(x))] + [(x[-1], i, 0), (x[0], i, 0)]
        poly = Poly3DCollection([verts], color=viridis(i), alpha=0.5)
        ax.add_collection3d(poly)

    ax.set_xlabel('Mass (Da)')
    ax.set_ylabel('Timepoint (min)')
    ax.set_zlabel('Intensity')

    # Set y-ticks and y-tick labels to show only alternate values
    yticks = range(len(data_list))
    yticklabels = [label if i % 2 == 0 else '' for i, (_, _, label) in enumerate(data_list)]
    ax.set_yticks(yticks)
    ax.set_yticklabels(yticklabels)

    ax.xaxis.set_major_formatter(FuncFormatter(format_mass_axis))
    ax.set_xlim(0, 1000000)


    # Remove background grids and shading
    ax.grid(False)
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))

    plt.savefig(output_file, format='pdf')
    plt.show()

if __name__ == "__main__":
    directory = "Z:\\Group Share\\Dash\\Abby Data\\Fig Filters\\Rep1_Text Files"
    output_file = os.path.join(directory, 'waterfall_plot.pdf') # Output file for the plot

    data_list = []
    for file_name in sorted(os.listdir(directory)):
        if file_name.endswith('.txt'):
            file_path = os.path.join(directory, file_name)
            x, y = read_data(file_path)
            label = file_name.split('.')[0]  # Use the file name (without extension) as the label
            data_list.append((x, y, label))

    create_waterfall_plot(data_list, output_file)