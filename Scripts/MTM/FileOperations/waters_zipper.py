import os
import re
import shutil

def zip_raw_dirs(topdir):
    for root, dirs, files in os.walk(topdir):
        for d in dirs:
            if re.search(r".*\.raw", d):
                print("Found:", d)
                zipdir = os.path.join(root, d)
                shutil.make_archive(zipdir, "zip", zipdir)
                print("Zipped:", d)
                # Remove original directory
                shutil.rmtree(zipdir)


def zip_unidecfiles_dirs(topdir):
    for root, dirs, files in os.walk(topdir):
        for d in dirs:
            if re.search(".*_unidecfiles", d):
                print("Found:", d)
                zipdir = os.path.join(root, d)
                shutil.make_archive(zipdir, "zip", zipdir)
                print("Zipped:", d)
                # Remove original directory
                shutil.rmtree(zipdir)

def zip_pattern_dirs(topdir, pattern):
    for root, dirs, files in os.walk(topdir):
        for d in dirs:
            if re.search(pattern, d):
                print("Found:", d)
                zipdir = os.path.join(root, d)
                shutil.make_archive(zipdir, "zip", zipdir)
                print("Zipped:", d)
                # Remove original directory
                shutil.rmtree(zipdir)

topdir = "Z:\\Group Share\\Group\\Papers\\Nanodisc_assembly\\20230214Nanodisc Assembly\\Raw Data"
topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\"
topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Melanie Odenkirk"
topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\"
topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Skippy"
topdir = "Z:\\mtmarty\\Archive\\Papers\\"
topdir = "Z:\\mtmarty\\Data\\Finished Projects"
topdir = "Z:\\Group Share\\Hiruni Jayasekera"  # Example for Hiruni's directory

zip_unidecfiles_dirs(topdir)
zip_raw_dirs(topdir)
# zip_pattern_dirs(topdir, ".*_datatemp.*")
# # zip_pattern_dirs(topdir, ".*_Raw$")
zip_pattern_dirs(topdir, ".*UniDec_Figures_and_Files$")

# zip_pattern_dirs(topdir, ".*\.Data$")
# zip_pattern_dirs(topdir, ".*old$")
# zip_pattern_dirs(topdir, ".*Old$")
# zip_pattern_dirs(topdir, ".*OLD$")
#
# topdir = "Z:\\Group Share\\Group\\Papers"
#
# zip_pattern_dirs(topdir, ".*old$")
# zip_pattern_dirs(topdir, ".*Old$")
# zip_pattern_dirs(topdir, ".*OLD$")
