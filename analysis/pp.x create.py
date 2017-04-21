import sys, os
from PyQt4 import QtGui, uic
from PyQt4.QtCore import Qt
from biscotti.classes import qecalc
from biscotti.functions.base import  list_folder
import pyqtgraph as pg
pg.setConfigOption('background', 'w')
pg.setConfigOption('foreground', 'k')

import logging
# Logging level by default
logger = logging.getLogger(__name__)
loglevel = logging.DEBUG
logger.setLevel(loglevel)
console_handler = logging.StreamHandler()
console_handler.setLevel(loglevel)
formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)
# End logging.

qtCreatorFile = r"pw.x_analysis.ui"  # Enter file here.

Ui_MainWindow, QtBaseClass = uic.loadUiType(qtCreatorFile)
# Widgets:
# filename_label    Qlabel
# load_button       QPushButton
# param_table       QTabelWidget
# plot_1            PlotWidget
# unit_combo_1      QComboBox
# free_energy_1     QCheckBox
# delta_1           QCheckbox
# plot_2            PlotWidget
# etc. for plot 2

class PP_Create(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self):
        # Setup calls
        QtGui.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)

        # Main data structures
        self.last_rootpath = None

    def read_folder(self):
        rootpath = QtGui.QFileDialog.getExistingDirectory(caption='Select calculation output folder',
                                               directory=self.last_rootpath)
        self.last_rootpath = rootpath
        pwscf_in = qecalc.QECalcIn.import_from_file(pwscf_in)


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    window = PP_Create()
    window.show()
    sys.exit(app.exec_())
