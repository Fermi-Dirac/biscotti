from biscotti.classes import atoms
from pyqtgraph.Qt import QtGui, QtCore
import pyqtgraph.opengl as gl
import numpy as np

app = QtGui.QApplication([])  # not sure if i can delete this or wtf.
def render_atoms(structure : atoms.AtomicStructure, window=gl.GLViewWidget(), add_bonds=False, radius=1, scale=1, **kwargs):
    meshdata = gl.MeshData.sphere(rows=10, cols=10, radius=radius*scale)
    colors = {'In': (1,0,0,1), 'As' : (0,1,0,1), 'Sb' : (0,0,1,1)}
    maxlist = np.array([0, 0, 0])
    minlist = np.array([0, 0, 0])
    for atom in structure.atomscart:
        color = colors[atom.species] if atom.species in colors else (1,1,1,1)
        atom_render = gl.GLMeshItem(meshdata=meshdata, smooth=True, color=color, shader='shaded', glOptions='opaque')
        x, y, z = atom.position * scale
        for i, new_val in enumerate([x,y,z]):
            if new_val > maxlist[i]:
                maxlist[i] = new_val
            elif new_val < minlist[i]:
                minlist[i] = new_val
        atom_render.translate(x, y, z)
        window.addItem(atom_render)
    center = (maxlist+minlist) / 2
    window.pan(center[0], center[1], center[2])
    return window
