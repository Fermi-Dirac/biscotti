# This is my silly pre-amble until i can get custom packages to work
import sys
custompackages = [r'D:\Users\Chris\Documents\SivaLab\Python', r'D:\Users\Chris\Documents\Python', r'~/bin/python']
sys.path.extend(custompackages)
# End Preamble

from biscotti.classes import atoms

import zipfile, os, zlib
from glob import glob

compression = zipfile.ZIP_DEFLATED

rootpath = os.path.abspath('')  # get the calling directory path

print ("root path is " + rootpath)
newroot = os.path.split(rootpath)[0]
batchname = os.path.split(rootpath)[1]
os.chdir(newroot)
print ("Creating archive for batch " + batchname)
zf = zipfile.ZipFile(batchname + ".zip", mode='w')

alloutfiles = [y for x in os.walk(batchname) for y in glob(os.path.join(x[0], '*.out'))]
alloutfiles += [y for x in os.walk(batchname) for y in glob(os.path.join(x[0], 'slurmout.txt'))]
for outfile in alloutfiles:
    try:
        print ("Adding " + outfile)
        zf.write(outfile, compress_type=compression)
    except:
        pass
zf.close()
print ("Done")