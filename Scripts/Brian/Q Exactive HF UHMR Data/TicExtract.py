import UniDecImporter.ImporterFactory as importer_factory
import os
import numpy as np

topdir = r"Z:\Group Share\BHT\Q Exactive HF Data\RPLC-MS\Agilent LC\20250828"
os.chdir(topdir)

for file in os.listdir(topdir):
    print(file)
    if file.lower().endswith(".raw"):
        print(f"Processing {file}...")
        filepath = os.path.join(topdir, file)
        importer = importer_factory.create_importer(filepath)
        tic = importer.get_tic()
        outname = os.path.splitext(file)[0] + "_TIC.txt"
        np.savetxt(outname, tic, header="Time\tTIC", delimiter="\t")
        print(f"TIC saved to {outname}")