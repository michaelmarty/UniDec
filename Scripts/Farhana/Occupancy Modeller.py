import numpy as np

R = 8.31446261815324 / 1000.

def correct_k(k, r):
    return (r/(1+r)) * k + (1/(1+r))

def calc_dg(K):
    return -R * 298.15 * np.log(K)

def inverse_dg(dg):
    return np.exp(-dg/(R*298.15))

ka = 7.5
kma = 5
k = ka/kma

dg = 6
k= inverse_dg(dg)

r = 0.5


dg = calc_dg(correct_k(k, r))
dg2 = calc_dg(k)
print("Raw DG:", dg2)
print("Observed DG:", dg)
print("Difference:", dg-dg2)
