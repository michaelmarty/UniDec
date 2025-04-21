import shutil
import os

def zip_subdirs(topdir):
    for d in os.listdir(topdir):
        dir = os.path.join(topdir, d)
        if os.path.isdir(dir):
            print("Found:", dir)
            shutil.make_archive(dir, "zip", dir)
            print("Zipped:", dir)
            # Remove original directory
            shutil.rmtree(dir)

def zip_lowest_directory(topdir):
    for root, dirs, files in os.walk(topdir):
        if len(dirs) == 0:
            print("Found:", root)
            shutil.make_archive(root, "zip", root)
            print("Zipped:", root)
            # Remove original directory
            shutil.rmtree(root)

topdir = "Z:\\Group Share\\Group\\Archive\\Undergrad\\October Owens"
topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Madison De Jesus"
topdir = "Z:\\mtmarty\\VMD\\AmtB"
# zip_subdirs(topdir)
# topdir = "Z:\\mtmarty\\VMD\\AmtB Nanodiscs"
# zip_subdirs(topdir)
# topdir = "Z:\\mtmarty\\VMD\\ANT1"
# zip_subdirs(topdir)
topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Melanie Odenkirk\\2022_Synapt_Experiments_Skyline_Processing\\Example Data from Z"
topdir = "Z:\\mtmarty\\Archive\\PostDoc\\LargeBackup\\Current-Done\\VMD\\AqpZ"
topdir = "Z:\\mtmarty\\Archive\\PostDoc\\LargeBackup\\Current-Done\\Python"
topdir = "Z:\\mtmarty\\Archive\\PostDoc\\Data\\OthersData"
topdir = "Z:\\mtmarty\\Archive\\Papers"
topdir = "Z:\\Group Share\\Hiruni Jayasekera"
zip_subdirs(topdir)
# zip_lowest_directory(topdir)