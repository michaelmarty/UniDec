import numpy as np
import pandas as pd
import unidec.tools as ud
import matplotlib.pyplot as plt
import statsmodels.api as sm
import os

dir = "Z:\\Group Share\\Annika\\APCI Data\\MSMS_ToFMRM Tests\\MS_MSMS_ToFMRM"
os.chdir(dir)
df = pd.read_excel("EIC_Intensities.xlsx")

types = ["ToFMRM", "MS", "MSMS"]
for t in types:
    df1 = df[df["Type"] == t]
    # filter for conc < 120
    df1 = df1[df1["Concentration"] < 120]

    # Make linear plots of concentration vs intensity
    fit = np.polyfit(df1["Concentration"], df1["Intensity"], 1)
    fit_fn = np.poly1d(fit)
    # calculate std error of regression
    lr = sm.OLS(df1["Intensity"], df1["Concentration"]).fit()
    stderr = float(lr.bse[0])
    conc_coeff = fit[0]
    conc_coeff = lr.params[0]
    lod = 3.3 * stderr / conc_coeff
    print(t,fit, stderr, conc_coeff, lod)
    #print(lr.summary())
    fit_fn = lr.predict


    plt.plot(df1["Concentration"], df1["Intensity"], "o", df1["Concentration"], fit_fn(df1["Concentration"]), label=t)

plt.legend()

plt.show()
