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
from copy import deepcopy


# import os
def MakeTexReport(fname, config, path, peaks, labels, names, color, figureflags):
    rawfilename = deepcopy(config.filename)
    # rawfilename = rawfilename.replace('\\', '\\\\')
    # rawfilename = rawfilename.replace(' ', ': ')
    rawfilename = rawfilename.replace('\\', '\\textbackslash ')
    rawfilename = rawfilename.replace('_', '\\_')

    header = deepcopy(config.outfname)
    # header = header.replace('\\', '\\\\')
    # header = header.replace('\\', '\\textbackslash ')
    header = header.replace('_', '\\_')
    # header = header.replace(' ', ': ')
    header2 = os.path.split(header)[1]

    path = deepcopy(path)

    error = config.error
    confname = config.confname

    textdict = dict(np.transpose([['\u25CB', '\u25BD', '\u25B3', '\u25B7', '\u25A2', '\u2662', '\u2606'],
                                  ['$\\bigcirc$', '$\\bigtriangledown$', '$\\bigtriangleup$', '$\\triangleright$',
                                   '$\\square$', '$\\lozenge$', '$\\bigstar$']]))

    f = open(fname, 'w+')
    f.write("\\documentclass{article}\n")
    f.write("\\usepackage{graphicx}\n")
    f.write("\\usepackage{amsmath}\n")
    f.write("\\usepackage{amssymb}\n")
    f.write("\\usepackage{paracol}\n")
    f.write("\\usepackage{color,colortbl}\n")
    f.write("\\usepackage[margin=0.75in]{geometry}\n")
    f.write("\\graphicspath{{" + path + "/}}\n")
    f.write("\\begin{document}\n")
    f.write("\\begin{center}\n")
    f.write("\\Large\\textbf{UniDec Report}\n")
    f.write("\\end{center}\n")

    f.write("\\noindent\n")
    f.write("Time: " + time.strftime("%a, %d %b %Y %H:%M:%S +0000", time.gmtime()) + "\\\\\n")
    f.write("File Name: {" + rawfilename + "}\\\\\n")
    f.write("Error: " + str(error) + "\\\\\n")
    # f.write("\\columnratio{0.5}\n")
    f.write("\\begin{paracol}{2}\n")
    scale = 0.48
    f.write("\\begin{minipage}{3.5in}\n")
    if np.any(np.array(figureflags) == 1):
        f.write("\\includegraphics[scale=" + str(scale) + "]{{" + header2 + "_Figure" + str(1) + "}.pdf}\\\\\n")
    f.write("\\end{minipage}\n")
    f.write("\\switchcolumn\n")
    f.write("\\begin{minipage}{3.5in}\n")
    if np.any(np.array(figureflags) == 3):
        f.write("\\includegraphics[scale=" + str(scale) + "]{{" + header2 + "_Figure" + str(3) + "}.pdf}\\\\\n")
    f.write("\\end{minipage}\n")
    f.write("\\end{paracol}\n")

    f.write("\\begin{center}\n")
    scale = 0.85
    if np.any(np.array(figureflags) == 5):
        f.write("\\includegraphics[scale=" + str(scale) + "]{{" + header2 + "_Figure" + str(5) + "}.pdf}\\\\\n")
    f.write("\\end{center}\n")

    f.write("\\newpage\n")
    f.write("\\begin{center}\n")
    f.write("\\Large\\textbf{Peak Information}\n")
    f.write("\\end{center}\n")

    f.write("\\columnratio{0.65}\n")
    f.write("\\begin{paracol}{2}\n")
    f.write("\\begin{minipage}{3.75in}\n")

    scale = 0.50
    for i in [4, 2, 6]:
        # if(i==0 or i==2 or i==4):
        #    f.write("\\includegraphics[scale="+str(scale)+"]{Luc_Figure"+str(i+1)+".pdf}&")
        if np.any(np.array(figureflags) == i):
            # print "Figure ",i
            f.write("\\includegraphics[scale=" + str(scale) + "]{{" + header2 + "_Figure" + str(i) + "}.pdf}\\\\\n")

    f.write("\\end{minipage}\n")

    f.write("\\switchcolumn\n")
    # f.write("\\\\\n")
    # f.write("\\\\\n")

    f.write("\\begin{table}[ht]\n")
    f.write("\\centering\n")
    f.write("\\begin{tabular}{c c c c}\n")
    f.write("Symbol & Mass (Da) & Intensity & Name \\\\\n")
    f.write("\\hline\n")
    f.write("\\hline\n")
    for i in range(0, len(peaks)):
        a = 1
        f.write(
            "\\rowcolor[rgb]{" + "%.4f" % color[i][0] * a + "," + "%.4f" % color[i][1] * a + "," + "%.4f" % color[i][
                2] * a + "}\n")
        f.write(textdict[labels[i]] + '&$' + str(peaks[i][0]) + '$&$' + str(peaks[i][1]) + '$&' + names[i] + '\\\\\n')
    f.write("\\hline\n")
    f.write("\\hline\n")
    f.write("\\end{tabular}\n")
    f.write("\\end{table}\n")

    f.write("\\end{paracol}\n")
    f.write("\\newpage\n")

    f.write("\\begin{center}\n")
    f.write("\\Large\\textbf{Parameters}\n")
    f.write("\\end{center}\n")
    f.write("\\noindent\n")
    conf = open(confname, 'r')
    # conf.readline()
    # conf.readline()
    for line in conf:
        linestr = line.rstrip()
        #linestr = linestr.replace('\\', '\\textbackspace ')
        linestr = linestr.replace('_', '\\_')
        linestr = linestr.replace(' ', ': ')
        if "input" in linestr or "output" in linestr or "kernel" in linestr or "manualfile" in linestr:
            # f.write("\\text{" + linestr + "}\\\\\n")
            # f.write("\\text{\detokenize{" + linestr + "}}\\\\\n")
            pass
        else:
            # print linestr
            f.write("\\text{" + linestr + "}\\\\\n")
    conf.close()

    f.write("\\end{document}\n")

    f.close()


def PDFTexReport(fname):
    path = "C:\\Program Files (x86)\\MiKTeX 2.9\\miktex\\bin\\pdflatex.exe"
    print(path, fname)
    oldpath = os.getcwd()
    newpath = os.path.dirname(fname)
    os.chdir(newpath)
    if os.path.isfile(path):
        call = [path, str(fname),"-halt-on-error"]
        print(call)
        subprocess.call(call)
    else:
        call = ["pdflatex", str(fname), "-halt-on-error"]
        print(call)
        subprocess.call(call)
    os.chdir(oldpath)
