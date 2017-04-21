from biscotti.classes import atoms
import os
from tkinter import filedialog
rootpath = filedialog.askopenfilename()
newstruct = atoms.AtomicStructure.from_QEinput(rootpath)
newfolder, newfile = os.path.split(rootpath)
newstruct.write_vasp(newfolder, 'POSCAR ' + newfile + ".vasp")