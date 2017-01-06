import sys
sys.path.append(r"D:\Users\Chris\Documents\SivaLab\Python")

import os
from classes import qecalc
path = r'D:\Users\Chris\Documents\SivaLab\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Batch004 - As-Sb Swap\InAs_2,2,5\InAs_2,2,5.in'
testcalc = qecalc.QECalcIn.import_from_file(path)

for namelist in testcalc.namelistdict:
    print("-" + namelist + "-")
    print(type(testcalc.namelistdict[namelist]))
    for key in testcalc.namelistdict[namelist]:
        print (" Key,Value ("+ key + ", " + str(testcalc.namelistdict[namelist][key]) + ")")

pathlist = os.path.split(path)
testcalc.write_to_file(pathlist[0], pathlist[1] + " (2).in")