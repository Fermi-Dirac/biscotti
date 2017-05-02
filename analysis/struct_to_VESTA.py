from biscotti.classes import atoms
from pyqtgraph.Qt import QtGui
import tempfile
import subprocess
from os.path import sep
vesta_dir = r'C:\Program Files (x86)\VESTA-win64\VESTA.exe'
app = QtGui.QApplication([])
calcpath = QtGui.QFileDialog.getOpenFileName(caption='Select calculation output file')
struct = atoms.AtomicStructure.from_QEinput(calcpath)
struct.write_vasp(folder=tempfile.gettempdir(), filename='POSCAR ' + struct.name + ".vasp")
subprocess.call('"' + vesta_dir + '" "'+ tempfile.gettempdir() + sep + 'POSCAR ' + struct.name + '.vasp"')