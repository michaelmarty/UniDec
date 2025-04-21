import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
from unidec.modules.isotopetools import *
from unidec.IsoDec.plots import *
from unidec.IsoDec.match import *
from matplotlib import rcParams


# Set backend to Agg
mpl.use('WxAgg')

rcParams['ps.useafm']=True
rcParams['ps.fonttype']=42
rcParams['pdf.fonttype']=42
rcParams['font.family'] = 'Arial'

#region generate mock data
mock_data = []

noise_factor = 0.2

dist1_factor = 0.3
dist1_original = fast_calc_averagine_isotope_dist(5800, 6)
dist1_original[:,1] *= 100
dist1_noise = np.random.normal(1,noise_factor*2,len(dist1_original))
dist1_original[:, 1] = dist1_original[:, 1] * dist1_factor
dist1_noise_added = np.zeros((len(dist1_original), 2))
for i in range(len(dist1_original)):
    dist1_noise_added[i, 0] = dist1_original[i, 0]
    dist1_noise_added[i, 1] = dist1_original[i, 1] * dist1_noise[i]
    mock_data.append(dist1_noise_added[i])

dist2_factor = 0.5
dist2_original = fast_calc_averagine_isotope_dist(972, 1)
dist2_original[:,1] *= 100
dist2_noise = np.random.normal(1,noise_factor*2,len(dist2_original))
dist2_original[:, 1] = dist2_original[:, 1] * dist2_factor
dist2_noise_added = np.zeros((len(dist2_original), 2))
for i in range(len(dist2_original)):
    dist2_noise_added[i, 0] = dist2_original[i, 0]
    dist2_noise_added[i, 1] = dist2_original[i, 1] * dist2_noise[i]
    mock_data.append(dist2_noise_added[i])

dist3_factor = 0.7
dist3_original = fast_calc_averagine_isotope_dist(3879.4, 4)
dist3_original[:,1] *= 100
dist3_noise = np.random.normal(1,noise_factor,len(dist3_original))
dist3_original[:, 1] = dist3_original[:, 1] * dist3_factor
dist3_noise_added = np.zeros((len(dist3_original), 2))
for i in range(len(dist3_original)):
    dist3_noise_added[i, 0] = dist3_original[i, 0]
    dist3_noise_added[i, 1] = dist3_original[i, 1] * dist3_noise[i]
    mock_data.append(dist3_noise_added[i])

##endregion generate mock data

xtick_size = 14
ytick_size = 14
xlabel_size = 14
ylabel_size = 14

ticks = [50,100]
tick_labels = ["%","100"]

ticks1 = [0,50,100]
tick_labels1 = [0,"%",100]
mz_range = [966,978]

save_figures = True
show_figures = True

mz_window = [-1.05,2.05]

#make the first cluster selection plot
if True:
    first_peak_index = np.argmax(dist3_noise_added[:,1])
    first_peak_intensity = dist3_noise_added[first_peak_index,1]
    first_peak_mz = dist3_noise_added[first_peak_index,0]

    #Create shaded regions on either side of the peak.
    min_mz = first_peak_mz + mz_window[0]
    max_mz = first_peak_mz + mz_window[1]

    r1_width = abs(mz_range[0] - min_mz)
    r1_height = 100
    #need a rectangle from the minimum of the mz_range to the min_mz
    r1 = patches.Rectangle((mz_range[0], 0), r1_width, r1_height, color='b', alpha=0.25)

    r2_width = abs(max_mz - mz_range[1])
    r2_height = 100
    #need a rectangle from the max_mz to the maximum of the mz_range
    r2 = patches.Rectangle((max_mz, 0), r2_width, r2_height, color='b', alpha=0.25)


    plt.figure(figsize=(3,3))
    cplot(np.array(mock_data), color='k')
    plt.plot(first_peak_mz, first_peak_intensity + 5, 'b*', markersize=20)
    plt.text(first_peak_mz + 0.4, first_peak_intensity - 5, s="+4", fontsize=xlabel_size, color='k')
    plt.gca().add_patch(r1)
    plt.gca().add_patch(r2)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.hlines(0, mz_range[0], mz_range[1], color='black', linewidth=1)
    plt.ylim(0,100)
    plt.yticks(ticks, tick_labels, fontsize=ytick_size)
    plt.xlim(mz_range)
    plt.xticks(fontsize=xtick_size)
    plt.xlabel("m/z", fontsize=xlabel_size)
    plt.tight_layout()
    if save_figures:
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\first_cluster_selection.png")
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\first_cluster_selection.pdf")
    if show_figures:
        plt.show()

#make the first shift optimization and scoring plot
if True:
    ##region Annotations
    first_peak_index = np.argmax(dist3_noise_added[:, 1])
    first_peak_intensity = dist3_noise_added[first_peak_index, 1]
    first_peak_mz = dist3_noise_added[first_peak_index, 0]

    # Create shaded regions on either side of the peak.
    min_mz = first_peak_mz + mz_window[0]
    max_mz = first_peak_mz + mz_window[1]

    r1_width = abs(mz_range[0] - min_mz)
    r1_height = 200
    # need a rectangle from the minimum of the mz_range to the min_mz
    r1 = patches.Rectangle((mz_range[0], -100), r1_width, r1_height, color='b', alpha=0.25)

    r2_width = abs(max_mz - mz_range[1])
    r2_height = 200
    # need a rectangle from the max_mz to the maximum of the mz_range
    r2 = patches.Rectangle((max_mz, -100), r2_width, r2_height, color='b', alpha=0.25)

    #get the min mz in the distribution
    min_dist3_mz = np.min(dist3_original[:,0])
    #get the max mz in the distribution
    max_dist3_mz = np.max(dist3_original[:,0])
    arrow_offset = 0.2
    arrow_length = 0.5
    ##endregion Annotations

    ##region Mock Scoring
    #Get the noise-added intensities
    dist3_noise_added_intensities = dist3_noise_added[:,1]
    #Get the original intensities
    dist3_original_intensities = dist3_original[:,1]
    css = calculate_cosinesimilarity(dist3_noise_added_intensities, dist3_original_intensities,shift=0,max_shift=0)
    print("Round 1, Cosine Similarity Score: ", css)

    #Calculate areacovered
    #Get the area within the mz_window of the mock data
    mock_data_within_window = []
    for i in range(len(mock_data)):
        if mock_data[i][0] >= min_mz and mock_data[i][0] <= max_mz:
            mock_data_within_window.append(mock_data[i])
    exp_area = np.sum(np.fromiter(([x[1] for x in mock_data_within_window]), dtype=float))
    matched_area = np.sum(np.fromiter((x[1] for x in dist3_noise_added if x[0] >= min_mz and x[0] <= max_mz), dtype=float))
    print("Round 1, Area Covered: ", matched_area / exp_area)
    print("Round 1, matched peaks: ", len(dist3_noise_added))
    ##endregion Mock Scoring

    ##region Plot
    plt.figure(figsize=(3,3))
    cplot(np.array(mock_data), color='k')
    cplot(np.array(dist3_original), color='b', factor=-1)
    plt.gca().add_patch(r1)
    plt.gca().add_patch(r2)
    plt.arrow(min_dist3_mz - arrow_offset, -10, -arrow_length, 0,head_width=5, head_length=0.1, linewidth=1, fc='k', ec='k')
    plt.arrow(max_dist3_mz, -10, arrow_length, 0, head_width=5, head_length=0.1, linewidth=1, fc='k', ec='k')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.hlines(0, mz_range[0], mz_range[1], color='black', linewidth=1)
    plt.ylim(-100, 100)
    plt.yticks(ticks1, tick_labels1, fontsize=ytick_size)
    plt.xlim(mz_range)
    plt.xticks(fontsize=xtick_size)
    plt.xlabel("m/z", fontsize=xlabel_size)
    plt.tight_layout()
    if save_figures:
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\first_shiftopt_scoring.png")
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\first_shiftopt_scoring.pdf")
    if show_figures:
        plt.show()
    ##endregion Plot

#make the first knockdown plot
if True:
    ##region Plot
    plt.figure(figsize=(3, 3))
    cplot(np.array(mock_data), color='k')
    cplot(np.array(dist3_noise_added), color='r')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.hlines(0, mz_range[0], mz_range[1], color='black', linewidth=1)
    plt.ylim(0, 100)
    plt.yticks(ticks, tick_labels, fontsize=ytick_size)
    plt.xlim(mz_range)
    plt.xticks(fontsize=xtick_size)
    plt.xlabel("m/z", fontsize=xlabel_size)
    plt.tight_layout()
    if save_figures:
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\first_knockdown.png")
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\first_knockdown.pdf")
    if show_figures:
        plt.show()
    ##endregion Plot

#make the second cluster selection plot
if True:
    mock_data_2 = []
    for i in range(len(dist1_noise_added)):
        mock_data_2.append(dist1_noise_added[i])
    for i in range(len(dist2_noise_added)):
        mock_data_2.append(dist2_noise_added[i])

    dist1_peak = np.argmax(dist1_noise_added[:,1])
    dist1_peak_mz = dist1_noise_added[dist1_peak,0]
    dist1_peak_intensity = dist1_noise_added[dist1_peak,1]

    dist2_peak = np.argmax(dist2_noise_added[:,1])
    dist2_peak_mz = dist2_noise_added[dist2_peak,0]
    dist2_peak_intensity = dist2_noise_added[dist2_peak,1]

    #Create shaded regions on either side of either peak.
    min_mz1 = dist1_peak_mz + mz_window[0]
    max_mz1 = dist1_peak_mz + mz_window[1]

    min_mz2 = dist2_peak_mz + mz_window[0]
    max_mz2 = dist2_peak_mz + mz_window[1]

    r1_width = abs(mz_range[0] - min_mz1)
    r1_height = 100

    r1 = patches.Rectangle((mz_range[0], 0), r1_width, r1_height, color='b', alpha=0.25)

    r2_width = min_mz2 - max_mz1
    r2_height = 100

    r2 = patches.Rectangle((max_mz1, 0), r2_width, r2_height, color='b', alpha=0.25)

    r3_width = abs(max_mz2 - mz_range[1])
    r3_height = 100

    r3 = patches.Rectangle((max_mz2, 0), r3_width, r3_height, color='b', alpha=0.25)

    plt.figure(figsize=(3, 3))
    cplot(np.array(mock_data_2), color='k')
    plt.plot(dist1_peak_mz, dist1_peak_intensity + 5, 'b*', markersize=20)
    plt.text(dist1_peak_mz + 0.4, dist1_peak_intensity - 5, s="+6", fontsize=xlabel_size, color='k')
    plt.plot(dist2_peak_mz, dist2_peak_intensity + 5, 'b*', markersize=20)
    plt.text(dist2_peak_mz + 0.4, dist2_peak_intensity - 5, s="+1", fontsize=xlabel_size, color='k')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().add_patch(r1)
    plt.gca().add_patch(r2)
    plt.gca().add_patch(r3)
    plt.hlines(0, mz_range[0], mz_range[1], color='black', linewidth=1)
    plt.ylim(0, 100)
    plt.yticks(ticks, tick_labels, fontsize=ytick_size)
    plt.xlim(mz_range)
    plt.xticks(fontsize=xtick_size)
    plt.xlabel("m/z", fontsize=xlabel_size)
    plt.tight_layout()
    if save_figures:
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\second_cluster_selection.png")
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\second_cluster_selection.pdf")
    if show_figures:
        plt.show()

#make the second shift optimization and scoring plot
if True:
    ##region Annotations
    dist1_peak = np.argmax(dist1_noise_added[:,1])
    dist1_peak_mz = dist1_noise_added[dist1_peak,0]
    dist1_peak_intensity = dist1_noise_added[dist1_peak,1]

    dist2_peak = np.argmax(dist2_noise_added[:,1])
    dist2_peak_mz = dist2_noise_added[dist2_peak,0]
    dist2_peak_intensity = dist2_noise_added[dist2_peak,1]

    #Create shaded regions on either side of either peak.
    min_mz1 = dist1_peak_mz + mz_window[0]
    max_mz1 = dist1_peak_mz + mz_window[1]

    min_mz2 = dist2_peak_mz + mz_window[0]
    max_mz2 = dist2_peak_mz + mz_window[1]

    r1_width = abs(mz_range[0] - min_mz1)
    r1_height = 200

    r1 = patches.Rectangle((mz_range[0], -100), r1_width, r1_height, color='b', alpha=0.25)

    r2_width = min_mz2 - max_mz1
    r2_height = 200

    r2 = patches.Rectangle((max_mz1, -100), r2_width, r2_height, color='b', alpha=0.25)

    r3_width = abs(max_mz2 - mz_range[1])
    r3_height = 200

    r3 = patches.Rectangle((max_mz2, -100), r3_width, r3_height, color='b', alpha=0.25)

    arrow_offset = 0.2
    arrow_length = 0.5
    ##endregion Annotations

    ##region Mock Scoring
    #Get the noise-added intensities
    dist1_noise_added_intensities = dist1_noise_added[:,1]
    #Get the original intensities
    dist1_original_intensities = dist1_original[:,1]
    css1 = calculate_cosinesimilarity(dist1_noise_added_intensities, dist1_original_intensities,shift=0,max_shift=0)
    print("Round 2, Cosine Similarity Score 1: ", css1)

    dist2_noise_added_intensities = dist2_noise_added[:,1]
    dist2_original_intensities = dist2_original[:,1]
    css2 = calculate_cosinesimilarity(dist2_noise_added_intensities, dist2_original_intensities,shift=0,max_shift=0)
    print("Round 2, Cosine Similarity Score 2: ", css2)

    #Calculate areacovered
    within_region1 = []
    for i in range(len(mock_data_2)):
        if mock_data_2[i][0] >= min_mz1 and mock_data_2[i][0] <= max_mz1:
            within_region1.append(mock_data_2[i])
    exp_area1 = np.sum(np.fromiter(([x[1] for x in within_region1]), dtype=float))
    matched_area1 = np.sum(np.fromiter((x[1] for x in dist1_noise_added if x[0] >= min_mz1 and x[0] <= max_mz1), dtype=float))
    print("Round 2, Area Covered 1: ", matched_area1 / exp_area1)
    print("Round 2, matched peaks 1: ", len(dist1_noise_added))

    within_region2 = []
    for i in range(len(mock_data_2)):
        if mock_data_2[i][0] >= min_mz2 and mock_data_2[i][0] <= max_mz2:
            within_region2.append(mock_data_2[i])
    exp_area2 = np.sum(np.fromiter(([x[1] for x in within_region2]), dtype=float))
    matched_area2 = np.sum(np.fromiter((x[1] for x in dist2_noise_added if x[0] >= min_mz2 and x[0] <= max_mz2), dtype=float))
    print("Round 2, Area Covered 2: ", matched_area2 / exp_area2)
    print("Round 2, matched peaks 2: ", len(dist2_noise_added))
    ##endregion Mock Scoring

    ##region Plot
    #only plot the portions of each original distribution within the mz_window
    dist1_original_within_window = []
    for i in range(len(dist1_original)):
        if dist1_original[i][0] >= min_mz1 and dist1_original[i][0] <= max_mz1:
            dist1_original_within_window.append(dist1_original[i])

    dist2_original_within_window = []
    for i in range(len(dist2_original)):
        if dist2_original[i][0] >= min_mz2 and dist2_original[i][0] <= max_mz2:
            dist2_original_within_window.append(dist2_original[i])

    min_dist1_mz = np.min(dist1_original[:,0])
    max_dist1_mz = np.max(dist1_original[:,0])

    min_dist2_mz = np.min(dist2_original[:,0])
    max_dist2_mz = np.max(dist2_original[:,0])


    plt.figure(figsize=(3, 3))
    cplot(np.array(mock_data_2), color='k')
    cplot(np.array(dist1_original_within_window), color='b', factor=-1)
    cplot(np.array(dist2_original_within_window), color='b', factor=-1)
    plt.arrow(min_dist1_mz - arrow_offset, -10, -arrow_length, 0,head_width=5, head_length=0.1, linewidth=1, fc='k', ec='k')
    plt.arrow(max_dist1_mz, -10, arrow_length, 0, head_width=5, head_length=0.1, linewidth=1, fc='k', ec='k')
    plt.arrow(min_dist2_mz - arrow_offset, -10, -arrow_length, 0, head_width=5, head_length=0.1, linewidth=1, fc='k', ec='k')
    plt.arrow(max_dist2_mz - 0.75, -10, arrow_length, 0, head_width=5, head_length=0.1, linewidth=1, fc='k', ec='k')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().add_patch(r1)
    plt.gca().add_patch(r2)
    plt.gca().add_patch(r3)
    plt.hlines(0, mz_range[0], mz_range[1], color='black', linewidth=1)
    plt.ylim(-100, 100)
    plt.yticks(ticks1, tick_labels1, fontsize=ytick_size)
    plt.xlim(mz_range)
    plt.xticks(fontsize=xtick_size)
    plt.xlabel("m/z", fontsize=xlabel_size)
    plt.tight_layout()
    if save_figures:
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\second_shiftop_scoring.png")
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\second_shiftop_scoring.pdf")
    if show_figures:
        plt.show()
    ##endregion Plot

#make the second knockdown plot
if True:
    ##region Plot
    plt.figure(figsize=(3, 3))
    cplot(np.array(mock_data_2), color='k')
    cplot(np.array(dist1_noise_added), color='r')
    cplot(np.array(dist2_noise_added), color='r')
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.hlines(0, mz_range[0], mz_range[1], color='black', linewidth=1)
    plt.ylim(0, 100)
    plt.yticks(ticks, tick_labels, fontsize=ytick_size)
    plt.xlim(mz_range)
    plt.xticks(fontsize=xtick_size)
    plt.xlabel("m/z", fontsize=xlabel_size)
    plt.tight_layout()
    if save_figures:
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\second_knockdown.png")
        plt.savefig("Z:\\Group Share\\JGP\\IsoDec_figs\\second_knockdown.pdf")
    if show_figures:
        plt.show()
    ##endregion Plot