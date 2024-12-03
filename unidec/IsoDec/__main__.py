import sys
import getopt
import argparse
import os
import platform
sys.path.append("..")

if platform.system() == "Windows":
    udpath = "C:\\Python\\UniDec3"
elif platform.system() == "Linux":
    udpath = "/home/u17/mtmarty/UniDec"

sys.path.append(udpath)


def main(*args, **kwargs):
    print("Running Command Line IsoDec")

    try:
        opts, args = getopt.getopt(sys.argv[1:], "f:o:", ["file=", "out="])
    except getopt.GetoptError as e:
        print("Error in Argv. Likely unknown option: ", sys.argv, e)
        print("Known options: -f, -o")
        return None

    # print("ARGS:", args)
    # print("KWARGS:", kwargs)
    # print("OPTS:", opts)
    infile = None
    outfile = None
    if opts is not None:
        for opt, arg in opts:
            if opt in ("-f", "--file"):
                infile = arg
                print("Opening File:", infile)
            if opt in ("-o", "--out"):
                outfile = arg
                print("Output File:", outfile)

            if os.path.isfile(infile):
                from unidec.IsoDec.runtime import IsoDecRuntime
                eng = IsoDecRuntime()
                eng.process_file(infile)
                return

    if len(args) == 1:
        file = args[0]
        print("Opening File:", file)
        if os.path.isfile(file):
            from unidec.IsoDec.runtime import IsoDecRuntime
            eng = IsoDecRuntime()
            eng.process_file(file)
            return

    # if "train" in args:
    #     from unidec.IsoDec.train import main
    #     main()
    #
    # elif "generate" in args:
    #     print("Generating Pkl Files")
    #     if infile is None:
    #         print("No Input File Specified")
    #         return None
    #     from unidec.IsoDec.trainingdata import process_file
    #     process_file(infile)
    #
    # elif "generate_script" in args:
    #     print("Generating Pkl Files")
    #     directory = os.getcwd()
    #     print("Directory:", directory)
    #     from unidec.IsoDec.generate import generate_all
    #     generate_all(directory)


if __name__ == '__main__':
    main()
