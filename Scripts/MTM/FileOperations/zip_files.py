import os
import shutil
import zipfile

def zip_files(topdir, ext):
    for root, dirs, files in os.walk(topdir):
        for f in files:
            if f.lower().endswith(ext):
                print("Found:", f)
                file = os.path.join(root, f)
                zipname = os.path.join(root, f + ".zip")
                # Compress
                with zipfile.ZipFile(zipname, "w", compression=zipfile.ZIP_DEFLATED) as z:
                    z.write(file, arcname=f)
                print("Zipped:", f)
                # Remove original file
                os.remove(file)



topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Skippy"
topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Larry Walker\\Papers\\Peptide_Paper2\\Bactenecin\\Data\\DMPG\\20200206"
topdir = "Z:\\Group Share\\James Rohrbough"
topdir = "Z:\\mtmarty\\Archive\\PostDoc\\Data\\PaperData\\NDData\\NDData_Raw"
#topdir = "Z:\\Group Share\\Seun\\LX-MS"

zip_files(topdir, ".raw")
#zip_files(topdir, ".hdf5")