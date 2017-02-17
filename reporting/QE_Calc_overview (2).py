#!/usr/bin/env python3
import os
import sys

cwd = os.getcwd()
calledfrom = os.path.abspath('') # get calling directory path
sys.path.append(cwd)

scriptdir = cwd

print("Called from " + calledfrom)
print("Now importing from " + scriptdir)
# os.chdir(scriptdir)
from classes import qecalc

os.chdir(cwd)
if len(sys.argv) > 1:
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print("Usage: \n Call without command to use current directory, otherwise:\n  QE_Calc_overview [rootpath] [outputfile] [delimiter]")
    else:
        rootpath = sys.argv[1]
        qecalc.makeSummaryFile(sys.argv[1], sys.argv[2], sys.argv[3])
else:
    qecalc.makeSummaryFile()
