import os
from glob import glob
from tkinter import filedialog
rootpath = filedialog.askdirectory()

alloutfiles1 = [y for x in os.walk(rootpath) for y in glob(os.path.join(x[0], '*.out'))] #obfuscated code
alloutfiles2 = []
ext_list = ['out']
for folder in os.walk(rootpath):
    alloutfiles2.extend([folder[0] + os.sep + file for file in os.listdir(folder[0]) if any(file.endswith(ext) for ext in ext_list)])

alloutfiles3= []
for folder in os.walk(rootpath):
    for file in os.listdir(folder[0]):
        if any([file.endswith(ext) for ext in ext_list]):
            alloutfiles3.append(folder[0] + os.sep + file)

print("Using glob, we found " + str(len(alloutfiles1)) + " files")
print("Using os.walk we found " + str(len(alloutfiles2)) + " files such as " + alloutfiles2[0])
print("Using os.walk 2 we found " + str(len(alloutfiles3)) + " files such as " + alloutfiles3[0])
