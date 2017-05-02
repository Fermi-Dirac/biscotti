from biscotti.classes import atoms
from tkinter import filedialog
import os
calcpath = filedialog.askopenfilename()
ext = os.path.splitext(calcpath)
print(ext)
if ext == '.in':
    struct =  atoms.AtomicStructure.from_QEinput(calcpath)
else:
    struct = atoms.AtomicStructure.from_QEOutput(calcpath)[-1]
struct.write_vasp(os.path.dirname(calcpath),"POSCAR " + os.path.basename(calcpath) + ".vasp")