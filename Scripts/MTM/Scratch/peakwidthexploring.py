import numpy as np
import matplotlib.pyplot as plt
import unidec.modules.fitting as fit
import matplotlib
from unidec.modules.fwhmtools import calc_sigma
from unidec.tools import peakdetect
from unidec.tools import safedivide

matplotlib.use("WxAgg")

def autocorr_ratio(deltas):
    s0 = np.sum(deltas**2)
    s1 = np.sum(deltas[:-1] * deltas[1:])
    if s0 == 0:
        return 1.0
    return s1 / s0

def area_ratio(data, template, threshold=0.1):
    b1 = template > np.amax(template) * threshold
    normfactor = np.amax(data[b1]) / np.amax(template[b1])
    area_data = np.sum(data[b1])
    area_template = np.sum(template[b1])
    if area_data == 0:
        return 1.0
    return area_template / area_data * normfactor

def apply_ratio(template, ratio):
    outt = template ** (ratio ** 2)
    # Normalize
    outt = outt / np.amax(outt)
    return outt

def delta_ratio(data, template, threshold=0.1):
    b1 = template > np.amax(template) * threshold
    delta_max = np.amax(data[b1])
    delta_sum = np.sum(data[b1])
    if delta_max == 0:
        return 1.0
    return delta_sum / delta_max

def template_decon(ydata, templates, threshold=0.1):
    ratios = np.array([area_ratio(ydata, t, threshold=threshold) for i, t in enumerate(templates)])
    ratios = np.dot(templates.T, ratios) / np.sum(templates, axis=0) # Convolve ratios
    return np.array([apply_ratio(t, r) for t, r in zip(templates, ratios)])

def apply_delta_ratios(deltas, templates, threshold=0.1):
    ratios = np.array([delta_ratio(deltas, t, threshold=threshold) for i, t in enumerate(templates)])
    # smooth ratios
    ratios = np.dot(templates.T, ratios) / np.sum(templates, axis=0)
    dmax = np.amax(deltas)
    out = deltas ** (ratios ** 0.1)
    out = out / np.amax(out) * dmax
    return out

def enforce_minmax(xdata, template, minmax):
    sigma = calc_sigma(np.column_stack((xdata, template)))
    if sigma < minmax[0]:
        template = apply_ratio(template, sigma/minmax[0])
    elif sigma > minmax[1]:
        template = apply_ratio(template, sigma/minmax[1])
    return template

def local_softmax(deltas, templates, beta=1):
    maxval = np.amax(deltas)
    if maxval * beta > 1000:
        return deltas  # Avoid overflow
    exps = np.exp(beta * deltas)
    localsums = np.dot(templates.T, exps)
    sim_data = np.dot(templates.T, deltas)
    return safedivide(exps * sim_data , localsums)

def richardson_lucy(data, minmax=[0.5, 3], beta=50, threshold=0.1, plot=False):
    # Create initial templates
    templates = np.array([fit.ndis_std(data[:,0], x, 2, 1) for x in data[:,0]])
    deltas = data[:,1] / np.sum(templates, axis=1)
    if plot:
        plt.plot(data[:, 0], deltas, label=f'Iteration {0}', color=plt.cm.viridis(0), alpha=0.5)

    # Adjust templates
    templates = template_decon(data[:, 1], templates, threshold=threshold)
    if minmax is not None:
        templates = np.array([enforce_minmax(data[:, 0], t, minmax) for t in templates])

    # Run Iteration
    maxit = 50
    lastdeltas = deltas.copy()
    for i in range(1, maxit):
        # Adjust intensities of deltas
        sim_data = np.dot(templates.T, deltas)
        correction = safedivide(data[:,1] , sim_data)
        correction = safedivide(np.dot(templates.T, correction) , np.sum(templates, axis=0))
        deltas = deltas * correction

        # maxpos = np.argmax(data[:,1])
        # deltaratios = np.array([delta_ratio(deltas, t, threshold=threshold) for i, t in enumerate(templates)])
        # print(deltaratios[maxpos])
        # deltas = apply_delta_ratios(deltas, templates, threshold=threshold)

        if np.amax(deltas) > np.amax(data[:,1]):
            deltas = deltas / np.amax(deltas) * np.amax(data[:,1])

        if beta > 0:
            # Apply softmax to deltas to stabilize
            acr = autocorr_ratio(deltas)
            deltas = local_softmax(deltas, templates, beta=(beta*acr)**2+beta)

        # Check convergence
        # diff = np.linalg.norm(deltas - lastdeltas) / np.linalg.norm(lastdeltas)
        # print(f"Iteration {i}, Diff: {diff}")
        # if diff < 1e-3:
        #     break
        # lastdeltas = deltas.copy()
        if plot:
            # maxpos = np.argmax(data[:,1])
            # plotdata = templates[maxpos] * data[maxpos, 1]
            # plotdata = sim_data
            # plotdata = deltas / deltas[maxpos] * data[maxpos, 1]
            plotdata = deltas
            # plotdata = correction
            # plotdata = deltaratios
            plt.plot(data[:, 0], plotdata, label=f'Iteration {i + 1}', color=plt.cm.viridis(i / maxit), alpha=0.5)

    sim_data = np.dot(templates.T, deltas)
    # sim_data /= np.amax(sim_data) * np.amax(data[:, 1])
    # Find peaks in deltas
    # peaks = peakdetect(np.column_stack((data[:,0], deltas)), threshold=0.2)
    # peakindexes = [np.argmin(np.abs(data[:,0]-p[0])) for p in peaks]
    # for pi in peakindexes:
    #     print("Peak at:", data[pi,0], "Value:", deltas[pi], "Sigma:", calc_sigma(np.column_stack((data[:,0], templates[pi]))))

    return deltas, sim_data

def gen_random_spectrum(x, npeaks=5, noiselevel=0.01):
    y = np.zeros_like(x)
    centers = np.random.uniform(np.min(x), np.max(x), npeaks)
    sigmas = np.random.uniform(0.5, 2.0, npeaks)
    amplitudes = np.random.uniform(0.1, 1.0, npeaks)
    datas = []
    for c, s, a in zip(centers, sigmas, amplitudes):
        d = fit.ndis_std(x, c, s, a)
        datas.append(d)
        y += d
        print("Peak at:", c, "Sigma:", s, "Amplitude:", a)
    y += np.random.normal(0, noiselevel, y.shape)
    y[y < 0] = 0
    return y, [centers, sigmas, amplitudes], datas

np.random.seed(0)

beta = 0
minmax = [0.5, 3]
minmax = None

x = np.linspace(-50, 50, 2000)
y, info, datas  = gen_random_spectrum(x, npeaks=1, noiselevel=0.02)

plt.figure(figsize=(20,10))
for i, d in enumerate(datas):
    plt.plot(x, d, 'g-', alpha=0.3)
    # Plot line at center
    c, s, a = info[0][i], info[1][i], info[2][i]
    plt.vlines(c, 0, a, colors='b', linestyles='dashed', alpha=0.5)

deltas, dfit = richardson_lucy(np.column_stack((x, y)), beta=beta, minmax=minmax, plot=True, threshold=0.25)
plt.plot(x, y, color="k", label='Original Data')
plt.plot(x, dfit, 'r--', label='Reconstructed Data')
plt.show()


exit()

y = fit.ndis_std(x, 0, 1, 1)
y += fit.ndis_std(x, 10, 0.5, 0.5)
y = y + np.random.normal(0, 0.01, y.shape)
y[y < 0] = 0
plt.figure(figsize=(20,10))
plt.subplot(121)
deltas, dfit = richardson_lucy(np.column_stack((x, y)), beta=beta, minmax=minmax, plot=True)

plt.plot(x, y, color="k", label='Original Data')
plt.plot(x, dfit, 'r--', label='Reconstructed Data')


plt.subplot(122)
x = np.linspace(-10, 20, 1000)
y = fit.ndis_std(x, 0, 1, 1)
y += fit.ndis_std(x, 10, 0.5, 0.5)
y = y + np.random.normal(0, 0.01, y.shape)
y[y < 0] = 0
deltas, dfit = richardson_lucy(np.column_stack((x, y)), beta=beta, minmax=minmax, plot=True)

plt.plot(x, y, color="k", label='Original Data')
plt.plot(x, dfit, 'r--', label='Reconstructed Data')
plt.show()