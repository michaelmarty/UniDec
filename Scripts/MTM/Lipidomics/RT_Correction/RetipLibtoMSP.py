import os
import RetipFormatter
import pandas as pd
import matplotlib.pyplot as plt
os.chdir(r"C:\Python\pyRetip\notebooks\lipidomics")

csvfile = r"Neg_C1_Ref_FromMSP_predicted.csv"
libfile = r"MSDIAL-TandemMassSpectralAtlas-VS69-Neg-C1.msp"

outfile = libfile.replace("C1.msp", "RTP1.msp")

df = pd.read_csv(csvfile)

plt.scatter(df['RT'], df['RTP'])
plt.xlabel("Experimental RT (min)")
plt.ylabel("Predicted RT (min)")
plt.plot([0, max(df['RT'])], [0, max(df['RT'])], 'r--')
plt.show()

RetipFormatter.df_to_msp(df, libfile, outfile, newcol="RTP", verbose=False)