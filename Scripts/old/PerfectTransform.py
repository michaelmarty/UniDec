import numpy as np
import matplotlib.pyplot as plt
from UniDec import unidec


def perfect_transform(mz, ztab, I):
    mtab = []
    for i, z in enumerate(ztab):
        masses = mz * z - z
        mtab.append(masses)
    fm = np.ravel(mtab)
    sorts = np.argsort(fm)

    fi = np.ravel(I)
    sfm = fm[sorts]
    sfi = fi[sorts]

    print(sfm, sfi)

    plt.plot(sfm, sfi, marker='o')
    plt.xlim(66000, 67000)
    plt.show()

    pass

path = "C:\\Data\\AlgorithmTesting\\BSA.txt"
eng = unidec.UniDec()
eng.open_file(path)
eng.unidec_imports(everything=True)

X = np.unique(eng.data.mzgrid[:,0])
Y = np.unique(eng.data.mzgrid[:,1])
Z = eng.data.mzgrid[:,2].reshape(len(X), len(Y)).transpose()

perfect_transform(X, Y, Z)

exit()


print(Z)
plt.contourf(X, Y, Z)
plt.show()



