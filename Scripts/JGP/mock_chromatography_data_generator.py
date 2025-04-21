import numpy as np
import matplotlib.pyplot as plt
import unidec.tools as ud
import unidec.IsoDec.datatools as dt
import matplotlib as mpl
import numpy as np


def gaussian(x, A, mu, sigma):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))


if __name__ == "__main__":
    # Set backend to Agg
    mpl.use('WxAgg')

    #need to set a run length and sampling rate
    run_length = 20 #minutes
    run_length_seconds = run_length * 60
    sampling_rate = 20 #Hz

    #generate a time array
    time = np.linspace(0, run_length_seconds, run_length_seconds * sampling_rate)
    time_minutes = time / 60
    #specify the location of the center of the peaks
    centers = np.array([7.8,12.4,15.23]) #minutes
    #specify the peak heights
    extinction_coefficients = np.array([15687, 35568, 55369])
    desired_areas = np.array([20,10,10])
    desired_areas = desired_areas * extinction_coefficients
    #specify the peak widths
    widths = np.array([0.2,0.3,0.4]) #minutes
    #generate the peak heights
    heights = desired_areas / (widths * np.sqrt(2 * np.pi))
    #generate the data
    intensity = np.zeros_like(time)
    real_areas = np.zeros_like(centers)
    for i in range(len(centers)):
        curr_intensity = gaussian(time, heights[i], centers[i]*60, widths[i]*60)
        #Calculate the area under the curve
        area = np.trapz(curr_intensity, time)
        real_areas[i] = area
        intensity += curr_intensity


    print(real_areas)
    scaled_areas = real_areas / extinction_coefficients
    print(scaled_areas)
    normalized_areas = scaled_areas / np.max(scaled_areas)
    print(normalized_areas)
    plt.plot(time_minutes, intensity)
    plt.xlabel("Time (minutes)")
    plt.ylabel("A280")
    plt.show()