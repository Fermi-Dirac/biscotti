from biscotti.classes import atoms
import os
from tkinter import filedialog
rootpath = filedialog.askopenfilename()
newstruct = atoms.AtomicStructure.from_QEinput(rootpath)
laststruct = atoms.AtomicStructure.from_QEOutput(rootpath)[-1]
newfolder, newfile = os.path.split(rootpath)
newstruct.write_vasp(newfolder, 'POSCAR ' + newfile + ".vasp")
laststruct.write_vasp(newfolder, 'POSCAR ' + newfile + "_last.vasp")