import pyqtgraph.opengl as gl
from biscotti.classes import atoms
from functions import atoms_plot
from pyqtgraph.Qt import QtCore, QtGui

calcpath = QtGui.QFileDialog.getOpenFileName(caption='Select calculation output file')
# calcpath = r'\\cluster\cbuurma\espresso\Superlattice\Batch013 - 3steps\InAsSb_[1,1,3,3],np=32,ndiag=16,nk=1\InAsSb_[1,1,3,3],np=32,ndiag=16,nk=1.in'
struct = atoms.AtomicStructure.from_QEinput(calcpath)

app = QtGui.QApplication([])
window = gl.GLViewWidget()
# window.opts['distance'] = 1
window.show()
window.setWindowTitle('pyqtgraph example: GLScatterPlotItem')
window.setCameraPosition(distance=30)

grid = gl.GLGridItem()
window.addItem(grid)

atoms_plot.render_atoms(struct, window)

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
