# This is my silly pre-amble until i can get custom packages to work
import sys
custompackages = [r'D:\Users\Chris\Documents\SivaLab\Python', r'D:\Users\Chris\Documents\Python', r'~/bin/python']
sys.path.extend(custompackages)
# End Preamble

from biscotti.classes import qecalc
from glob import glob
import os
import matplotlib.pyplot as plt

rootpath = r'D:\Users\Chris\Documents\SivaLab\Type II Superlattices\2016 MDA Type 2 SL\Ab-Initio\Superlattice\Batch005 - As-Sb Swap small'
alloutfiles = [y for x in os.walk(rootpath) for y in glob(os.path.join(x[0], '*.out'))]

file = alloutfiles[3]
ctobj = qecalc.CalcTime.get_time(file)

for type in ['CPU', 'WALL', 'calls']:
    ctobj.pie_charts(os.path.split(file)[1] + " " + type, type)

# print (ctobj.startdt, ctobj.enddt)
# for file in alloutfiles:
#     ctobj = qecalc.CalcTime.get_time(file)
#     ctobj.pie_charts(os.path.split(file)[1], type = 'WALL')

plt.show()