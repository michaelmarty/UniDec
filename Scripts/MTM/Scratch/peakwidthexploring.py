import numpy as np
import matplotlib.pyplot as plt
import unidec.modules.fitting as fit
import matplotlib
from unidec.modules.fwhmtools import calc_sigma
from unidec.tools import peakdetect

matplotlib.use("WxAgg")

def area_ratio(data, template, threshold=0.1):
    b1 = template > np.amax(template) * threshold
    area_data = np.sum(data[b1])
    area_template = np.sum(template[b1])
    return area_template / area_data

def apply_ratio(template, ratio):
    return template ** (ratio**2)

def template_decon(data, templates, threshold=0.1):
    ratios = np.array([area_ratio(data[:,1], t, threshold=threshold) for t in templates])
    # b1 = ratios > 1.5
    # ratios[b1] = 1.5
    adjusted_templates = [apply_ratio(t, r) for t, r in zip(templates, ratios)]
    return np.array(adjusted_templates)

def enforce_minmax(xdata, template, minmax):
    sigma = calc_sigma(np.column_stack((xdata, template)))
    if sigma < minmax[0]:
        ratio = minmax[0] / sigma
        template = apply_ratio(template, 1/ratio)
    elif sigma > minmax[1]:
        ratio = minmax[1] / sigma
        template = apply_ratio(template, 1/ratio)
    return template

def richardson_lucy(data, minmax=[0.5, 3]):
    templates = np.array([fit.ndis_std(data[:,0], x, 2, 1) for x in data[:,0]])
    deltas = data[:,1] #/ np.sum(templates, axis=1)
    plt.plot(data[:, 0], deltas, label=f'Iteration {0}', color=plt.cm.viridis(0), alpha=0.5)
    maxit = 40
    for i in range(1, maxit):
        # Adjust templates
        sim_data = np.dot(templates.T, deltas)
        templates = template_decon(np.column_stack((data[:, 0], sim_data)), templates, threshold=0.05)
        templates = np.array([enforce_minmax(data[:, 0], t, minmax) for t in templates])

        # Adjust intensities of deltas
        sim_data = np.dot(templates.T, deltas)
        correction = data[:,1] / sim_data
        deltas = deltas * correction # np.dot(templates.T, correction)
        # deltas = deltas * np.dot(templates.T, correction)
        plt.plot(data[:,0], deltas, label=f'Iteration {i+1}', color=plt.cm.viridis(i / maxit), alpha=0.5)


    sim_data = np.dot(templates.T, deltas)

    # Find peaks in deltas
    peaks = peakdetect(np.column_stack((data[:,0], deltas)), threshold=0.2)
    print(peaks)

    peakindexes = [np.argmin(np.abs(data[:,0]-p[0])) for p in peaks]
    for pi in peakindexes:
        print("Peak at:", data[pi,0], "Value:", deltas[pi], "Sigma:", calc_sigma(np.column_stack((data[:,0], templates[pi]))))

    sigmas = [calc_sigma(np.column_stack((data[:,0], t))) for t in templates]

    plt.plot(data[:,0], sigmas)

    return deltas, sim_data

np.random.seed(0)

x = np.linspace(-10, 20, 100)

y = fit.ndis_std(x, 0, 1, 1)
y += fit.ndis_std(x, 10, 0.5, 0.5)
y = y + np.random.normal(0, 0.01, y.shape)
y[y < 0] = 0

f1 = fit.ndis_std(x, 0, 2, 1)
#
ratio = area_ratio(y, f1)
print(ratio)
#
# f2 = apply_ratio(f1, ratio)

deltas, fit = richardson_lucy(np.column_stack((x, y)))

plt.plot(x, y)
# plt.plot(x, f1, 'r--', label='Fitted Peak')
# plt.plot(x, f2, 'g--', label='Adjusted Peak')
# plt.plot(x, deltas, 'b-', label='Deconvolved Data')
plt.plot(x, fit, 'k--', label='Reconstructed Data')
# plt.legend()
plt.show()