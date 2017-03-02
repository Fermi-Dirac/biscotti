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

alloutpaths = []
ext_list = ['out','in','slurmout.txt']
for folder in os.walk(batchname):
    alloutpaths.extend([folder[0] + os.sep + file for file in os.listdir(folder[0]) if any(file.endswith(ext) for ext in ext_list)])

print(str(len(alloutpaths)) + " files found")
for outfile in alloutpaths:
    try:
        print ("Adding " + outfile)
        zf.write(outfile, compress_type=compression)
    except:
        print("Cannot add " + outfile)
        pass
zf.close()
print ("Done")