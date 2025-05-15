import os, shutil, re

def find_files(topdir, filename):
    foundfiles = []
    for root, dirs, files in os.walk(topdir):
        if filename in files:
            print("Found:", os.path.join(root, filename))
            foundfiles.append(os.path.join(root, filename))
    return foundfiles

def delete_files(files):
    for f in files:
        try:
            os.remove(f)
            print("Deleted", f)
        except Exception as e:
            print("Could not delete", f)
            print(e)

# Write a function to find all directories with a given name in the top dir
def find_dirs(topdir, dirname):
    founddirs = []
    for root, dirs, files in os.walk(topdir):
        if dirname in dirs:
            founddirs.append(os.path.join(root, dirname))
            print("Found:", os.path.join(root, dirname))
    return founddirs

# Write a function to delete all dirs
def delete_dirs(dirs):
    for d in dirs:
        try:
            shutil.rmtree(d)
            print("Deleted", d)
        except Exception as e:
            print("Could not delete", d)
            print(e)

def run_delete_dirs(topdir, dirname):
    dirs = find_dirs(topdir, dirname)

    # Get user input on whether to delete
    print("Found", len(dirs), "directories to delete.")
    if len(dirs) > 0:
        print("Delete all?")
        response = input("y/n: ")

        if response == "y":
            delete_dirs(dirs)

def run_delete_files(topdir, filename):
    files = find_files(topdir, filename)
    print(files)
    # Get user input on whether to delete
    print("Found", len(files), "files to delete.")
    if len(files) > 0:
        print("Delete all?")
        response = input("y/n: ")

        if response == "y":
            delete_files(files)

def delete_empty_dirs(topdir):
    count = 0
    for root, dirs, files in os.walk(topdir):
        if len(dirs) == 0 and len(files) == 0:
            try:
                os.rmdir(root)
                print("Deleted Empty:", root)
                count += 1
            except Exception as e:
                print("Could not delete", root)
                print(e)
    print("Deleted", count, "empty directories.")

def find_and_zip(topdir, regex):
    count = 0
    # Find dir in topdir that matches regex
    for root, dirs, files in os.walk(topdir):
        for d in dirs:
            if re.search(regex, d):
                print("Found:", d)
                zipdir = os.path.join(root, d)
                shutil.make_archive(zipdir, "zip", zipdir)
                print("Zipped:", d)
                # Remove original directory
                shutil.rmtree(zipdir)
                count+=1

    print("Zipped", count, "directories.")

def remove_file_regex(topdir, regex, remove=True):
    count = 0
    for root, dirs, files in os.walk(topdir):
        for f in files:
            if len(f) == 0:
                continue
            # if regex matches the file name
            if re.search(regex, f) and remove:
                os.remove(os.path.join(root, f))
                print("Removed:", f)
                count += 1
            elif re.search(regex, f):
                print("Matched:", f)
                count += 1

    print("Removed", count, "files.")

def remove_dir_regex(topdir, regex, remove=True):
    count = 0
    for root, dirs, files in os.walk(topdir):
        for f in dirs:
            if len(f) == 0:
                continue
            # if regex matches the file name
            if re.search(regex, f) and remove:
                shutil.rmtree(os.path.join(root, f))
                print("Removed:", f)
                count += 1
            elif re.search(regex, f):
                print("Matched:", f)
                count += 1

    print("Removed", count, "files.")

def remove_duplicate_files(topdir, good=".*\.raw", bad=".*\.mzML", remove=False):
    count = 0
    for root, dirs, files in os.walk(topdir):
        for d in dirs:
            if re.search(good, d):
                # Check if there is a bad file in the directory
                for f in os.listdir(root):
                    if re.search(bad, f):
                        if os.path.splitext(f)[0] == os.path.splitext(d)[0]:
                            print("Found Match:", d, f)
                            if remove:
                                os.remove(os.path.join(root, f))
                                print("Removed:", f)
                                count += 1

    print("Removed", count, "files.")

topdir = "Z:\\Group Share\\Group\\Archive\\"
topdir = "Z:\\Group Share\\Group\\VMD\\AqpZMutants\\WT_AA_TOCL2_2\\namd"
topdir = "Z:\\Group Share\\Group\\VMD"
topdir = "Z:\\Instrument"
topdir = "Z:\\mtmarty\\Archive\\PostDoc\\"
topdir = "Z:\\Group Share\\Hiruni Jayasekera"
topdir = "C:\\Python\\UniDecDev\\UniDec3"
topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Kimber Focke\\More Synapt"

# topdir = "Z:\\Group Share\\Group"
# topdir = "Z:\\Group Share\\Group\\Archive\\Grad Students and Postdocs\\Marius Kostelic"
# topdir = "Z:\\mtmarty\\VMD"

#remove_duplicate_files(topdir, good=".*\.raw", bad=".*\.mzML", remove=True)
#delete_empty_dirs(topdir)
# run_delete_files(topdir, "Thumbs.db")
# run_delete_files(topdir, ".DS_Store")



# run_delete_dirs(topdir, r"blank")
# run_delete_dirs(topdir, "UniDecPastedSpectra")
# run_delete_files(topdir, ".*.mzML$")

remove_dir_regex(topdir, r'^.*blank.*', remove=True)

# remove_file_regex(topdir, r'^.*.trr', remove=True)
# remove_file_regex(topdir, r'^#.*#$', remove=True)
# remove_file_regex(topdir, r'^.*.old$', remove=True)
# remove_file_regex(topdir, r'^.*-DESKTOP-', remove=True)
# remove_file_regex(topdir, r'^.*.BAK$', remove=True)
# find_and_zip(topdir, ".*_koffs_.*")



