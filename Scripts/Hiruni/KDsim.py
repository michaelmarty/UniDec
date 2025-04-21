import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.special import binom

def min_free(ptot, ltot, kds, n, pguess=None, lguess=None):
    if pguess is None:
        pguess = ptot / 2
    if lguess is None:
        lguess = ltot / 2
    p0 = np.array([pguess, lguess])
    # fit = opt.leastsq(get_free_error, p0, args=(ptot, ltot, kds, n))[0]
    fit = opt.fsolve(get_free_error, p0, args=(ptot, ltot, kds, n))
    pfree = np.clip(fit[0], 0, ptot)
    lfree = np.clip(fit[1], 0, ltot)
    return pfree, lfree


def get_free_error(pO, ptot, ltot, kds, n):
    pfree = pO[0]
    lfree = pO[1]
    intgrid, sumprot, sumlig = make_grid(pfree, lfree, kds, n)
    return [(ptot - sumprot) ** 1, (ltot - sumlig) ** 1]


def make_grid(pfree, lfree, kds, n):
    kdgrid = np.array([kds for i in range(0, n)]).transpose()
    ligtab = np.arange(1, n + 1)

    intgrid = (pfree * (lfree ** ligtab)) / kdgrid ** ligtab

    sumprot = np.sum(intgrid)
    sumlig = np.sum(ligtab * intgrid)
    return intgrid, sumprot, sumlig


def fraction_bound(ptot, ltot, kd):
    a = (ptot + ltot + kd)
    fb = (a - np.sqrt(a ** 2 - 4 * ptot * ltot)) / (2 * ptot)
    return fb


def binding_grid(ptot, ltot, kds, n):
    pfree, lfree = min_free(ptot, ltot, kds, n)
    intgrid, sumprot, sumlig = make_grid(pfree, lfree, kds, n)
    return intgrid, pfree, lfree


def gen_spectrum(ptot, ltot, kds, n):
    intgrid, pfree, lfree = binding_grid(ptot, ltot, kds, n)
    pvec = pfree * np.ones((len(kds), 1))
    intgrid2 = np.concatenate((pvec, intgrid), axis=1)
    intgrid2 /= np.sum(intgrid2)

    print(intgrid2)

    ligtab = np.arange(0, n + 1)
    kdindex = range(0, len(kds))
    dims = tuple([n+1 for i in kdindex])

    indexes = list(np.ndindex(dims))

    intlist = np.array([[intgrid2[j, i[j]] for j in kdindex] for i in indexes])
    intlist2 = np.prod(intlist, axis=1)

    print(intlist2)

    liglist = np.array([[ligtab[i[j]] for j in kdindex] for i in indexes])
    liglist2 = np.sum(liglist, axis=1)

    unis = np.unique(liglist2)
    unilist = np.array([[i, np.sum(intlist2[liglist2 == i])] for i in unis])

    unilist[:, 1] /= np.sum(unilist[:, 1])

    plt.plot(unilist[:, 0], unilist[:, 1])
    plt.show()

    return


print(gen_spectrum(2, 65, [75, 200], 4))

exit()
x = np.linspace(0.1, 100, 50)
kd = 10
y = fraction_bound(x, 1, kd)

plt.plot(x, y)
plt.xlabel('Ligand concentration')
plt.ylabel('Fraction bound')
plt.ylim(0, 1)
plt.show()
