import os
import pandas as pd

toppath = "Z:\\Group Share\\JGP\\Lipidomics_MassIVE_folder\\Raw\\MS1\\WatersRaw\\"
os.chdir(toppath)

filenames = os.listdir(toppath)

df = pd.DataFrame(filenames, columns=["Filename"])
outputfile = 'output.xlsx'
df.to_excel(outputfile, index=False)
