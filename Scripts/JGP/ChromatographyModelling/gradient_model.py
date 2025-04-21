'''
k' =  (tR-t0)/t0
tR = retention time (seconds)
t0 = dead time
k' = retention factor

rearrange...
tR = k't0 + t0

estimate velocity of the solute in the column
V = L/tR = L/(k't0 + t0)

under the LSS model:
logki = logk0 + b(t/t0)


'''

import numpy as np
import matplotlib.pyplot as plt
import math

def determine_rt(x_vals, t, column_length):
    for i in range(len(x_vals)):
        if x_vals[i] >= column_length:
            return t[i]

def calculate_gaussian_signal(rt, sigma,t, scale=1):
    signal = np.zeros(len(t))
    for i in range(len(signal)):
        signal[i] = scale*np.exp(-((np.square(t[i]-rt))/np.square(2*sigma)))
    return signal

def calculate_signal(gradient, t, L, S, dt, kws, Ds, abs, delay=0):
    #Adding a delay value so that I can have pulsed injections at set intervals
    indices_delayed = int(delay/dt)

    ks = [[0 for x in range(len(t))] for y in range(len(kws))]

    for j in range(len(kws)):
        for i in range(indices_delayed, len(t)):
            ks[j][i] = np.exp(np.log(kws[j]) - S*gradient[i])

    ks_smooth = [[0 for x in range(len(t))] for y in range(len(kws))]
    window_size = 20

    i = 0
    while i < len(t) - window_size + 1:
        for j in range(len(kws)):
            k_j_window = np.zeros(window_size)
            for ind in range(window_size):
                k_j_window[ind] = ks[j][i+ind]
            ks_smooth[j][i] = np.mean(k_j_window)
        i+= 1

    vs = [[0 for x in range(len(t))] for y in range(len(kws))]
    for i in range(len(t)):
        for j in range(len(kws)):
            vs[j][i] = L/(ks_smooth[j][i]*t0 + t0)

    xs = [[0 for x in range(len(t))] for y in range(len(kws))]
    for j in range(len(kws)):
        xs[0][j] = 0

    for i in range(indices_delayed, len(t)):
        for j in range(len(kws)):
            if i == 0:
                continue
            xs[j][i] = xs[j][i-1] + vs[j][i-1]*dt

    sum_signal = np.zeros(len(t))
    for j in range(len(kws)):
        rt = determine_rt(xs[j],t,L)
        sigma = np.sqrt(2*Ds[j]*rt*60)
        signal = calculate_gaussian_signal(rt,sigma,t,abs[j])
        sum_signal += signal

    return sum_signal

#column parameters (set somewhat randomly..)
L = 30 #cm
S = 3
total_grad_length = 30
dt = 0.005
num_t = int(total_grad_length / dt)

t = np.linspace(0,total_grad_length,num_t)
t0 = 1 #dead time, minutes
gradient = np.zeros(num_t)


#this will be used to track and populate various arrays.
current_index = 0
##Let's start with a square wave-style gradient.
step_duration = 0.2
fB_min = 0.2
fB_max = 0.8

step_period = step_duration * 2

num_periods = int(total_grad_length / step_period)
points_per_step = int(step_duration/dt)

for i in range(num_periods):
    for j in range(points_per_step):
        gradient[current_index] = fB_min
        current_index += 1
    for j in range(points_per_step):
        gradient[current_index] = fB_max
        current_index += 1


##Analyte parameters
kw_1 = 20 #these are k values for isocratic separation in pure water (generally extrapolated)
kw_2 = 30

D1 = 1 * (10**-6) #cm2/s
D2 = 1.12 * (10**-6) #centimeters-squared/second


'''
I believe I have a decent pass at modelling a single injection
**It would probably be best to incorporate some sort of lag term with respect
to large changes in the retention of a specific species**

I did this in a somewhat lazy manner by just running a sliding window averaging
of k values to smooth them out..

Now I would like to try to provide a hadamard injection sequence and produce an
encoded chromatogram.
'''

inj_seq = np.array([0,0,0,0,1,0,0,1,0,1,1,0,0,1,1,1,1,1,0,0,0,1,1,0,1,1,1,0,1,0,1])

injection_delay = 0.4 #minutes, must fit with the step duration..

total_signal = np.zeros(num_t)

for i in range(len(inj_seq)):
    if inj_seq[i] == 1:
        signal = calculate_signal(gradient, t, L, S, dt, [kw_1,kw_2], [D1,D2], abs=[1,0.5], delay=i*injection_delay)
        total_signal += signal

plt.plot(t,total_signal)
plt.xlabel("Time (minutes)")
plt.ylabel("Signal")
plt.show()

