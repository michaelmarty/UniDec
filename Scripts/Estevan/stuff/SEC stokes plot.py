import matplotlib.pyplot as plt
from stokes_plotter import plot_stokes_with_error

cal_D = [1.7, 3.5, 5.62, 10.52, 17]          # known Stokes diameters (nm)
cal_RV = [20.79, 19.05, 18.09, 16.82, 13.81]    # retention volumes

sample_RV_means = [16.99, 13.50, 12.03, 11.40, 12.40, 11.34]
sample_RV_stds  = [0.09, 0.04, 0.06, 0.12, 0.04, 0.01]

plot_stokes_with_error(
    cal_RV=cal_RV,
    cal_D=cal_D,
    sample_RV_means=sample_RV_means,
    sample_RV_stds=sample_RV_stds,
    sample_labels=["MSP1D1", "MSP1D1-sfGFP", "MSP1D1-sfCherryShort", "MSP1D1-mCherry", "MSP1D1-mBlueberry", "MSP1D1-HaloTag"],
    sample_colors=["black", "green", "purple", "pink", "blue", "grey"],
)
plt.show()
