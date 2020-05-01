# -*- coding: utf-8 -*-
"""
Created on Wed Apr 02 15:52:23 2014

@author: Michael.Marty
"""
import subprocess
import time
import string
import numpy as np
import os

def MakeTexReport(fname, config, path, peaks, labels, names, color, figureflags, output):
    rawfilename = config.filename
    header = config.outfname
    error = config.error
    confname = config.confname
    
    textdict = dict(np.transpose([['\u25CB', '\u25BD', '\u25B3', '\u25B7', '\u25A2', '\u2662', '\u2606'],
                                  ['$\\largecircle$', '$\\medtriangledown$', '$\\medtriangleup$', '$\\medtriangleright$',
                                   '$\\largesquare$', '$\\medlozenge$', '$\\medwhitestar$']]))

    f = open(fname, 'w+')
    f.write("\\documentclass{article}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{amsmath}\n")
    # f.write("\\usepackage{amssymb}\n")- Andrew - using other symbol package for shapes
    f.write("\\usepackage{fdsymbol}\n")
    f.write("\\usepackage{paracol}\n")
    f.write("\\usepackage{color,colortbl}\n")
    f.write("\\usepackage[margin=0.5in]{geometry}\n")
    f.write("\\graphicspath{{" + path + "/}}\n")

    f.write("\\begin{document}\n")
    f.write("\\begin{center}\n")
    f.write("\\Large\\textbf{UniDec Deconvolution}\n")
    # Andrew - Remove "Baker Lab" and just include in the report info string
    f.write("\\end{center}\n")
    f.write("\\noindent\n")

    #Andrew - Use this for manual addtion of report information/MS details
    output1 = output.replace(';', '}\\\{')
    output2 = output1.replace('; ', '}\\\{')
    output3 = output2.replace('_', '\\_')
    f.write("\\text{" + output3 + "}\\\\\n")

    # Andrew - Below is what I used to include some info. This was just off of some code already in here. I imagine it could be better.
    # MSDetails.txt must always be included in the path folder
    #conf = open(path + "\MSDetails.txt", 'r')
    #for line in conf:
    #    linestr = line.rstrip()
    #    f.write("\\text{" + linestr + "}\\\\\n")
    #conf.close()
    
    f.write("Time Stamp: " + time.strftime("%a, %d %b %Y, %H:%M %Z", time.localtime()) + "\\\\\n")
    rawfilename = rawfilename.replace('_', '\\_')
    f.write("File Name: " + rawfilename + "\\\\\n")
    # f.write("\\columnratio{0.5}\n")
    f.write("\\begin{paracol}{2}\n")
    scale = 0.5
    f.write("\\begin{minipage}{3.5in}\n")
    if np.any(np.array(figureflags) == 4 ):
        f.write("\\includegraphics[scale=" + str(scale) + "]{{" + header + "_Figure" + str(4) + "}.pdf}\\\\\n")
    f.write("\\end{minipage}\n")
    f.write("\\switchcolumn\n")
    f.write("\\begin{minipage}{3.5in}\n")
    if np.any(np.array(figureflags) == 2):
        f.write("\\includegraphics[scale=" + str(scale) + "]{{" + header + "_Figure" + str(2) + "}.pdf}\\\\\n")
    f.write("\\end{minipage}\n")
    f.write("\\end{paracol}\n")

    f.write("\\begin{center}\n")
    f.write("\\begin{table}[ht]\n")
    f.write("\\centering\n")
    f.write("\\begin{tabular}{c c c c c}\n")
    f.write("Symbol & Mass (Da) & Intensity & Name/Species & Expected (Da) \\\\\n")
    f.write("\\hline\n")
    f.write("\\hline\n")
    for i in range(0, len(peaks)):
        a = 1
        f.write(
            "\\cellcolor[rgb]{" + "%.4f" % color[i][0] * a + "," + "%.4f" % color[i][1] * a + "," + "%.4f" % color[i][
                2] * a + "}\n")
        f.write(textdict[labels[i]] + '&$' + str(peaks[i][0]) + '$&$' + str(round(peaks[i][1],1)) + '$&$' + names[i] + '$&' "-" '\\\\\hline\n')
    f.write("\\hline\n")
    f.write("\\end{tabular}\n")
    f.write("\\end{table}\n")
    f.write("\\end{center}\n")
    f.write("\\end{document}\n")

    f.close()


def PDFTexReport(fname):
    path="C:\\Program Files (x86)\\MiKTeX 2.9\\miktex\\bin\\pdflatex.exe"
    if os.path.isfile(path):
        subprocess.call([path,fname,"-halt-on-error"])
    else:
        subprocess.call(["pdflatex", fname, "-halt-on-error"])
