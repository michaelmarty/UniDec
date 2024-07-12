import os
from unidec.IsoDec.trainingdata import *
import platform

slurmheader = """#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4GB
#SBATCH --time=8:00:00
#SBATCH --job-name=JOBNAME
#SBATCH --account=mtmarty
#SBATCH --partition=standard
#SBATCH --output=%x-%j.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtmarty@arizona.edu

cd $SLURM_SUBMIT_DIR

module load python
source /home/u17/mtmarty/mypyenv/bin/activate

ISODEC="/home/u17/mtmarty/UniDec/unidec/IsoDec"

python3 $ISODEC -f FILENAME generate > JOBNAME.log

"""


def generate(file):
    print(file)
    # Just the file name without the path or extension
    fileheader = os.path.basename(file).split(".")[0]
    newheader = slurmheader.replace("FILENAME", file)
    newheader = newheader.replace("JOBNAME", fileheader)

    outfile = fileheader + ".slurm"
    with open(outfile, "w") as f:
        f.write(newheader)

    print("Wrote:", outfile)
    return outfile


def generate_all(directory):
    os.chdir(directory)
    outfiles = []
    script = ""
    for file in os.listdir(directory):
        if is_valid(file):
            outfile = generate(file)
            outfiles.append(outfile)
            sline = "sbatch " + outfile + "\n"
            script += sline
    print(script)
    with open("submit.sh", "w") as f:
        f.write(script)


if __name__ == "__main__":
    if platform.system() == "Windows":
        directory = "C:\\Data\\TabbData"
    elif platform.system() == "Linux":
        directory = "/xdisk/mtmarty/mtmarty/training"

    generate_all(directory)
