import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glycodec_eng as gd
import matplotlib
matplotlib.use('WxAgg')

def get_matrix(col, df, target, tolerance):
    m0 = df[col].to_numpy()
    matrix = np.add.outer(m0, m0)
    tri = np.tri(matrix.shape[0], matrix.shape[1]).astype(bool)
    b0 = np.abs(matrix - target) < tolerance
    b0 = b0 * tri
    matrix[np.logical_not(b0)] = 0
    return matrix, m0, b0

def get_diff_matrix(col, df, target, tolerance):
    m0 = df[col].to_numpy()
    matrix = np.subtract.outer(m0, m0)
    tri = np.tri(matrix.shape[0], matrix.shape[1]).astype(bool)
    b0 = np.abs(matrix - target) < tolerance
    b0 = b0 * tri
    matrix[np.logical_not(b0)] = 0
    return matrix, m0, b0

def get_nlinked_matrix(df):
    m0 = df["NCore"].to_numpy()
    matrix = np.multiply.outer(m0, m0)
    tri = np.tri(matrix.shape[0], matrix.shape[1]).astype(bool)
    matrix = matrix * tri
    return matrix

enzyme_sequence = ["ABS", "BTG", "GUH"]

filtpath = "filtered_all.csv"
tolerance = 2
dtol = 2

m0target = 4411+18+18
dm1target = -1166
dm2target = -649
dm3target = -813

df = pd.read_csv(filtpath)
print(len(df))

m0matrix, m0, b0 = get_matrix("Mass0", df, m0target, tolerance)
dm1matrix, dm1, b1 = get_matrix("DM1", df, dm1target, dtol)
dm2matrix, dm2, b2 = get_matrix("DM2", df, dm2target, dtol)
dm3matrix, dm3, b3 = get_matrix("DM3", df, dm3target, dtol)
diffmatrix, m4, b4 = get_diff_matrix("Mass0", df, 0, 5)
b5 = get_nlinked_matrix(df)

ball = np.all([b0,b1, b2, b3, b5], axis=0)
xind, yind = np.where(ball)
goodpairs = np.transpose([xind, yind])
print(goodpairs)
print(len(goodpairs))

index = 0
if len(goodpairs)>0:
    for i in range(len(goodpairs)):
        g1 = df.iloc[goodpairs[i,0]]
        g2 = df.iloc[goodpairs[i,1]]

        if i == 0:
            print(g1)
            print(g1["IUPAC Condensed"])

        m0 = g1["Mass0"]+g2["Mass0"]
        m1 = g1["Mass1"]+g2["Mass1"]
        m2 = g1["Mass2"] + g2["Mass2"]
        m3 = g1["Mass3"] + g2["Mass3"]
        masses = [m1, m2, m3]

        dm1 = g1["DM1"] + g2["DM1"]
        dm2 = g1["DM2"] + g2["DM2"]
        dm3 = g1["DM3"] + g2["DM3"]
        dmasses = [dm1, dm2, dm3]

        eng = gd.GlycoDecEngine()
        eng.load_text(g1["IUPAC Condensed"])
        bad = eng.test_forbidden()

        eng2 = gd.GlycoDecEngine()
        eng2.load_text(g2["IUPAC Condensed"])
        bad2 = eng2.test_forbidden()

        if not bad and not bad2:
            index += 1
            fig = plt.figure()
            ax = plt.subplot(4,2,1)
            ax.set_title(str(m0) + " Da")

            eng.plot_struct(ax)


            for i, e in enumerate(enzyme_sequence):
                eng.enzyme(e, overwrite_structure=True)
                ax = plt.subplot(4, 2, 1+2*(i+1))
                eng.plot_struct(ax)
                ax.set_title(str(round(masses[i],3)) + " Da  (" + str(round(dmasses[i],3)) + ")")

            ax2 = plt.subplot(422)

            eng2.plot_struct(ax2)
            for i, e in enumerate(enzyme_sequence):
                eng2.enzyme(e, overwrite_structure=True)
                newax = plt.subplot(4, 2, 2+2*(i+1))
                eng2.plot_struct(newax)
print(index)
plt.show()


