#!/usr/bin/env python3
from Main import Parsers
import sys
if len(sys.argv) > 1:
    if sys.argv[1] == '-h' or sys.argv[1] == '--help':
        print("Usage: \n Call without command to use current directory, otherwise:\n  QE_Calc_overview [rootpath] [outputfile] [delimiter]")
    else:
        rootpath = sys.argv[1]
        Parsers.makeSummaryFile(sys.argv[1], sys.argv[2], sys.argv[3])
else:
    Parsers.makeSummaryFile()
