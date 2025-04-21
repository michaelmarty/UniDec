import os
import pandas as pd
import datetime
import ftplib
import time

pride_dir_base = "pride/data/archive/"
pride_url = "ftp.pride.ebi.ac.uk"

massive_url = "massive.ucsd.edu"
massive_dir_base = "v02/"

suppported = ["PRIDE"]


def get_url(row):
    id = row["identifier"]
    date = row["announce_date"]
    date = datetime.datetime.strptime(date, "%Y-%m-%d")
    year = date.year
    # Get month with 0 on the front
    month = '{:02d}'.format(date.month)

    repository = row["repository"]
    if repository == "PRIDE":
        url = f"{pride_dir_base}{year}/{month}/{id}/"
        return pride_url, url
    elif repository == "MassIVE":
        url = f"{massive_dir_base}{id}/"
        return massive_url, url
    else:
        print("Unsupported repository:", repository)
        return None, None


url = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/06/PXD045560/"


def list_files(base, directory):
    ftp = ftplib.FTP(base)
    ftp.login("anonymous", "anonymous")
    ftp.cwd(directory)
    files = ftp.nlst()
    print(files)
    ftp.quit()
    return files


def getFile(ftp, filename):
    try:
        ftp.retrbinary("RETR " + filename, open(filename, 'wb').write)
    except Exception as e:
        print("Error", filename)
        print(e)


file = "C:\\Data\\IsoNN\\proteomexchange_search.tsv"
file = "Z:\\Group Share\\JGP\\proteomexchange_search.tsv"
outdir = "Z:\\Group Share\\JGP\\"
df = pd.read_csv(file, sep="\t")
print(df.columns)
downloadsize = 0
maxsize = 1e12

fileext = ".raw"

for i, row in df.iterrows():
    if int(i) < 17:
        continue
    if row["repository"] not in suppported:
        continue
    id = row["identifier"]
    print(i, id)

    base, directory = get_url(row)

    ftp = ftplib.FTP(base)
    ftp.login("anonymous", "anonymous")
    ftp.cwd(directory)
    files = ftp.nlst()
    ftp.quit()
    print(files)

    if len(files) == 0:
        continue
    # If no files end in .raw continue
    if not any([f.endswith(fileext) for f in files]):
        continue

    newdir = os.path.join(outdir, id)
    if not os.path.exists(newdir):
        os.makedirs(newdir)
    os.chdir(newdir)
    print(newdir)

    for f in files:
        if f.endswith(fileext):
            ftp = ftplib.FTP(base)
            ftp.login("anonymous", "anonymous")
            ftp.cwd(directory)

            if os.path.isfile(f):
                continue
            size = ftp.size(f)
            print("Clock Time:", time.ctime())
            print(f, "File Size:", size/1e9, "GB Total Download Size:", downloadsize/1e9, "GB")
            downloadsize += size
            if downloadsize > maxsize:
                print("Max size reached")
                break
            getFile(ftp, f)

            ftp.quit()
    if downloadsize > maxsize:
        break
